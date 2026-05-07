import argparse
import configparser
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import reads_alignment
import unmapped_analysis
import reads_preprocess
import reference_prepare
import utils
import variant_calling
import report_generator
import summary_generator
import impurities_prefilter
import cleanup
import db_manager
import e2e_verifier


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class Task:
    def __init__(self):
        self.name = ''


def check_reads_file(task):
    if task.ex_r1 != None and task.ex_r2 != None:
        if Path(task.ex_r1).is_file() and Path(task.ex_r2).is_file():
            return 1
        else:
            logger.error('Reads file not found.')
            return -1
    else:
        logger.error('Reads file path can not be empty.')
        return -1


def check_ref_file(task):
    if Path(task.ref).is_file():
        return True
    else:
        logger.info('Reference sequence file not found.')
        return False


def check_deps(task):
    sys_deps = ['wget', 'git', 'apt', 'conda', 'python3', 'gzip', 'makeblastdb']
    if utils.sys_deps_check(sys_deps) == -1:
        logger.critical('System depency check fail.')
        sys.exit(100)
    if utils.conda_deps_check(task.conda_pkgs) == -1:
        logger.critical('Conda pkg depency check fail.')
        sys.exit(100)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ex_r1', help="Read-R1.")
    parser.add_argument(
        '--ex_r2', help="Read-R2.")
    parser.add_argument(
        '--prefix', help="For output prefix.", default='newtask')
    parser.add_argument(
        '--ref', help="Reference FASTA file path.", default=None)
    parser.add_argument(
        '--threads', help="CPU threads.", default=6)
    parser.add_argument(
        '--alns', help="Reads mapper list.", default='bowtie2,bwa')
    parser.add_argument(
        '--global_trimming', help="Global trimming bases for reads.", default=0)
    parser.add_argument(
        '--remove_host', help="Remove specific host genome (human, dog, vero, chicken, rhesus_monkey).", default=None)
    parser.add_argument(
        '--remove_impurities', help="Remove specific impurity sequences FASTA file path.", default=None)
    parser.add_argument(
        '--test', default=None)
    parser.add_argument(
        '--spades_mem', help="The memory (GB) allocated for spades, apply to both ref and unmapped assemble.", default=22)
    parser.add_argument(
        '--spades_mode', default='metaviral')
    parser.add_argument(
        '--unmapped_spades_mode', default='meta')
    parser.add_argument(
        '--unmapped_bbnorm', help="Enable bbnorm.sh normalization before unmapped reads assemble.", default='False')
    parser.add_argument(
        '--unmapped_bbnorm_target', help="Target coverage for bbnorm.sh.", default='30')
    parser.add_argument(
        '--unmapped_bbnorm_min', help="Min coverage for bbnorm.sh.", default='2')
    parser.add_argument(
        '--min_vc_score', default=1)
    parser.add_argument(
        '--vc_threshold', default='0.7')
    parser.add_argument(
        '--blastdb_path', default=None)
    parser.add_argument(
        '--rvdb_anno_path', default=None)
    parser.add_argument(
        '--unmapped_assemble', help="De novo Assemble the unmapped reads via metaSPAdes. ONLY apply to the first ref alignment.", default='True')
    parser.add_argument(
        '--unmapped_blastdb', help="BLASTDB for reference prepare and unmapped reads assemble.", default=None)
    parser.add_argument(
        '--unmapped_blastdb_extra_list', help="Extra custom BLASTDB list for unmapped reads assemble. Use a single space to seperate DB names.", default=None)
    parser.add_argument(
        '--unmapped_len_filter', help="Min. length (bp) filter to hit in unmapped reads assemble BLAST.", default='500')
    parser.add_argument(
        '--unmapped_ident_filter', help="Min. identity (%) filter to hit in unmapped reads assemble BLAST.", default='95')
    parser.add_argument(
        '--preset_path', help="Load VIVA analysis setting from given preset file path.", default=None)
    parser.add_argument(
        '--task_note', help="Task note. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_product_name', help="Sample (product) name. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_product_lot', help="Sample (product) lot. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_sequencing_date', help="Sample sequencing date. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_note', help="Sample note. Anotation purpose only.", default=None)
    parser.add_argument(
        '--auto_cleanup', help="Automatically clean up intermediate files after pipeline finished.", default='True')

    return parser


def get_latest_rvdb_files(blastdb_path):
    """偵測 blastdb_path 中版號最高的 RVDB 檔案"""
    import re
    if not blastdb_path:
        return None, None, None
    
    p = Path(blastdb_path)
    if not p.is_dir():
        return None, None, None
    
    # 搜尋 [U|C]-RVDBvXX.0.fasta (包含 .gz)
    fasta_pattern = re.compile(r'([UC])-RVDBv(\d+)\.0\.fasta(?:\.gz)?$')
    # 搜尋包含 RVDBvXX 且有 annotation 字樣的 .tab 檔案
    anno_pattern = re.compile(r'.*RVDBv(\d+).*[Aa]nnotation.*\.tab')
    
    max_ver = -1
    best_fasta = None
    
    # 優先權：版號大 > U-RVDB > C-RVDB
    # 先找出最高版號
    for f in p.glob('*RVDBv*.0.fasta*'):
        match = fasta_pattern.match(f.name)
        if match:
            ver = int(match.group(2))
            if ver > max_ver:
                max_ver = ver
    
    if max_ver != -1:
        # 在最高版號中挑選 U (優先) 或 C，且必須符合正則表達式（排除 .nhr 等）
        for f in p.glob(f'*RVDBv{max_ver}.0.fasta*'):
            if fasta_pattern.match(f.name):
                if f.name.startswith('U-'):
                    best_fasta = f.name
                    break
                elif f.name.startswith('C-'):
                    best_fasta = f.name
        
        if best_fasta:
            # 檔名處理：如果是 .gz，回傳解壓後的名稱供 blastdbcmd 使用
            clean_fasta_name = best_fasta.replace('.gz', '')
            
            # 尋找對應版號的 annotation
            best_anno = None
            for f in p.glob(f'*RVDBv{max_ver}*[Aa]nnotation*.tab*'):
                best_anno = str(f.absolute())
                break
            
            # 尋找額外的 C-RVDB (若目前是 U)
            best_extra = None
            if best_fasta.startswith('U-'):
                c_rvdb_name = f'C-RVDBv{max_ver}.0.fasta'
                if (p / c_rvdb_name).exists() or (p / (c_rvdb_name + '.gz')).exists():
                    best_extra = c_rvdb_name
                
            return clean_fasta_name, best_anno, best_extra
    
    return None, None, None


def fix_permissions(path):
    """將路徑下的所有檔案權限改回與父目錄一致 (解決 Docker root 權限問題)"""
    import os
    import subprocess
    try:
        # 嘗試從 /app/tasks 或傳入路徑的父目錄獲取宿主機使用者的 UID/GID
        target_path = Path(path)
        if not target_path.exists():
            return
            
        base_dir = '/app/tasks' if os.path.exists('/app/tasks') else str(target_path.parent)
        stat_info = os.stat(base_dir)
        uid = stat_info.st_uid
        gid = stat_info.st_gid
        
        # 只有在目前是 root 的情況下才需要改權限
        if os.getuid() == 0:
            subprocess.run(['chown', '-R', f'{uid}:{gid}', str(target_path)], check=True)
            logger.info(f"Fixed permissions for {target_path} to {uid}:{gid}")
    except Exception as e:
        logger.debug(f"Permission fix skipped for {path}: {e}")


def run_e2e_tests(input_args):
    """執行全量端到端測試並彙總報告"""
    test_types = ['ref', 'multi_ref', 'rvdb']
    all_results = []
    summary_lines = ["# VIVA E2E Test Suite Report\n"]
    
    for t_type in test_types:
        logger.info(f"===== Starting E2E Sub-test: {t_type} =====")
        # 複製參數並替換 --test 值
        sub_args = []
        skip_next = False
        for i, arg in enumerate(input_args):
            if skip_next:
                skip_next = False
                continue
            if arg == '--test':
                sub_args.append('--test')
                sub_args.append(t_type)
                skip_next = True
            else:
                sub_args.append(arg)
        
        if '--test' not in sub_args:
            sub_args.extend(['--test', t_type])
            
        try:
            task_id = main(sub_args)
            if task_id:
                # main 內部現在會自動調用驗證，我們在這裡抓取結果
                # 這裡為了彙總，我們手動再驗一遍獲取 md
                _, md = e2e_verifier.run_verification(
                    task_id, 
                    Path.cwd().joinpath('tasks'),
                    Path.cwd().joinpath('test_data', 'expected_results.json')
                )
                summary_lines.append(md)
            else:
                summary_lines.append(f"### Test Results for {t_type}\n**Overall Status: ⚠️ SKIPPED**\nReason: Test requirement not met (e.g. missing BLASTDB)\n")
            summary_lines.append("\n---\n")
        except Exception as e:
            logger.error(f"E2E Sub-test {t_type} failed: {e}")
            summary_lines.append(f"### Test Results for {t_type}\n**Overall Status: ❌ CRITICAL FAILURE**\nError: {e}\n")
            summary_lines.append("\n---\n")

    report_path = Path.cwd().joinpath('tasks', 'e2e_test_report.md')
    with open(report_path, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    fix_permissions(report_path)
    fix_permissions(Path.cwd().joinpath('tasks'))
    
    logger.info(f"E2E Test Suite finished. Report saved to {report_path}")
    print(f"\nE2E Test Suite Report generated at: {report_path}")


def main(input_args):
    parser = get_parser()
    args, unknown = parser.parse_known_args(input_args)

    task = Task()
    task.conda_pkgs = [
        'conda', 'python', 'perl',
        'fastp', 'samtools', 'bcftools', 'htslib',
        'bowtie2', 'bwa',
        'varscan', 'lofreq',
        'spades', 'blast', 'bbmap'
    ]
    check_deps(task)
    task.path = Path.cwd().joinpath('tasks')
    task.name = args.prefix
    task.task_note = args.task_note
    task.id = ''
    task.with_ref = False
    task.ex_r1 = args.ex_r1
    task.ex_r2 = args.ex_r2
    task.alns = args.alns.split(',')
    task.ref_num = 0
    task.impurities_prefilter_num = 0
    task.total_reads_after_fastp = 0
    task.preset_path = args.preset_path
    task.sample_product_name = args.sample_product_name
    task.sample_product_lot = args.sample_product_lot
    task.sample_sequencing_date = args.sample_sequencing_date
    task.sample_note = args.sample_note
    if task.preset_path == None:
        task.ref = args.ref
        task.threads = str(args.threads)
        task.global_trimming = str(args.global_trimming)
        task.remove_host = args.remove_host
        task.remove_impurities = args.remove_impurities
        task.spades_mem = str(args.spades_mem)
        task.spades_mode = args.spades_mode
        task.vc_threshold = args.vc_threshold
        task.min_vc_score = args.min_vc_score
        task.blastdb_path = args.blastdb_path
        task.rvdb_anno_path = args.rvdb_anno_path
        task.unmapped_assemble = args.unmapped_assemble
        task.unmapped_spades_mode = args.unmapped_spades_mode
        task.unmapped_bbnorm = args.unmapped_bbnorm
        task.unmapped_bbnorm_target = args.unmapped_bbnorm_target
        task.unmapped_bbnorm_min = args.unmapped_bbnorm_min
        task.unmapped_blastdb = args.unmapped_blastdb
        task.unmapped_blastdb_extra_list = args.unmapped_blastdb_extra_list
        task.unmapped_len_filter = args.unmapped_len_filter
        task.unmapped_ident_filter = args.unmapped_ident_filter
        task.auto_cleanup = args.auto_cleanup
    else:
        # Parse all conf. as strings
        config = configparser.ConfigParser(allow_no_value=True)
        config.read(args.preset_path)
        task.ref = config['PRESET']['ref']
        task.threads = str(config['PRESET']['threads'])
        task.global_trimming = str(config['PRESET']['global_trimming'])
        task.remove_host = config['PRESET']['remove_host']
        task.remove_impurities = config['PRESET']['remove_impurities']
        task.spades_mem = str(config['PRESET']['spades_mem'])
        task.spades_mode = config['PRESET']['spades_mode']
        task.vc_threshold = config['PRESET']['vc_threshold']
        task.min_vc_score = config['PRESET']['min_vc_score']
        task.unmapped_assemble = config['PRESET']['unmapped_assemble']
        task.unmapped_spades_mode = config['PRESET']['unmapped_spades_mode']
        task.unmapped_bbnorm = config['PRESET']['unmapped_bbnorm']
        task.unmapped_bbnorm_target = config['PRESET']['unmapped_bbnorm_target']
        task.unmapped_bbnorm_min = config['PRESET']['unmapped_bbnorm_min']
        task.blastdb_path = config['PRESET']['blastdb_path']
        task.rvdb_anno_path = config['PRESET']['rvdb_anno_path']
        task.unmapped_blastdb = config['PRESET']['unmapped_blastdb']
        task.unmapped_blastdb_extra_list = config['PRESET']['unmapped_blastdb_extra_list']
        task.unmapped_len_filter = config['PRESET']['unmapped_len_filter']
        task.unmapped_ident_filter = config['PRESET']['unmapped_ident_filter']
        task.auto_cleanup = config['PRESET']['auto_cleanup']
        task.preset_id = config['VERSION']['preset_id']
        task.preset_version = config['VERSION']['version']
        task.preset_last_rev_date = config['VERSION']['last_rev_date']
        task.preset_author = config['VERSION']['author']
        task.preset_note = config['VERSION']['note']

    if args.test == 'e2e':
        run_e2e_tests(input_args)
        return

    if args.test != None:
        task.name = 'test_run'
        task.ex_r1 = Path.cwd().joinpath('test_data', 'AdV_R1.fastq.gz')
        task.ex_r2 = Path.cwd().joinpath('test_data', 'AdV_R2.fastq.gz')
        if args.test == 'ref':
            task.name = 'test_ref'
            task.ref = Path.cwd().joinpath('test_data', 'AC_000008.1.fasta')
            task.remove_impurities = Path.cwd().joinpath('test_data', 'impure_test.fasta')
        elif args.test == 'multi_ref':
            task.name = 'test_multi_ref'
            task.ref = Path.cwd().joinpath('test_data', 'adv_multi_ref.fasta')
        elif args.test == 'denovo':
            task.name = 'test_denovo'
            task.remove_host = 'human'
            task.ref = None
            
            # De novo mode requires a BLASTDB to pick a reference.
            import glob
            search_paths = ['/app/blastdb', os.path.expanduser('~/bioapp/blastdb')]
            search_paths.extend(glob.glob('/home/*/bioapp/blastdb'))
            
            found_db = False
            for sp in search_paths:
                fasta, anno, extra = get_latest_rvdb_files(sp)
                if fasta:
                    task.blastdb_path = sp
                    task.unmapped_blastdb = fasta
                    task.rvdb_anno_path = anno
                    logger.info(f"Detected latest RVDB version for denovo test at {sp}: {fasta}")
                    found_db = True
                    break
            
            if not found_db:
                logger.warning("No BLASTDB found for denovo test. Skipping this test.")
                return None
                
        elif args.test == 'rvdb':
            task.name = 'test_rvdb'
            task.ref = Path.cwd().joinpath('test_data', 'AC_000008.1.fasta')
            
            # 自動偵測最新版 RVDB
            import glob
            search_paths = ['/app/blastdb', os.path.expanduser('~/bioapp/blastdb')]
            search_paths.extend(glob.glob('/home/*/bioapp/blastdb'))
            
            found_db = False
            for sp in search_paths:
                fasta, anno, extra = get_latest_rvdb_files(sp)
                if fasta:
                    task.blastdb_path = sp
                    task.unmapped_blastdb = fasta
                    task.rvdb_anno_path = anno
                    task.unmapped_blastdb_extra_list = f"{extra} core_nt" if extra else "core_nt"
                    logger.info(f"Detected latest RVDB version at {sp}: {fasta}")
                    found_db = True
                    break
            
            if not found_db:
                logger.warning("No BLASTDB found for rvdb test. Skipping this test.")
                return None


    if task.unmapped_blastdb != None:
        logger.info('Checking BlastDB.')
        if utils.setup_blastdb(task.blastdb_path, task.unmapped_blastdb) == -1:
            logger.error('BlastDB setup error. Exiting pipeline.')
            sys.exit()
        if task.unmapped_blastdb_extra_list != None:
            for db in task.unmapped_blastdb_extra_list.split():
                if utils.setup_blastdb(task.blastdb_path, db) == -1:
                    logger.error('Extra BlastDB setup error. Exiting pipeline.')
                    sys.exit()
        task.unmapped_assemble = 'True'
    
    if task.rvdb_anno_path != None:
        if not Path(task.rvdb_anno_path).is_file():
            logger.error('RVDB annotation file not found. Exiting pipeline.')
            sys.exit()

    logger.info('Checking reference.')
    if task.ref != None:
        if check_ref_file(task):
            task.with_ref = True
        else:
            logger.error('Input reference not found. Exiting pipeline.')
            sys.exit()
    else:
        logger.info('Input reference not provided. Will go de novo')

    if task.with_ref == False:
        if task.unmapped_blastdb != None:
            logger.info('Checking BlastDB.')
            if utils.setup_blastdb(task.blastdb_path, task.unmapped_blastdb) == -1:
                logger.error('BlastDB setup error. Exiting pipeline.')
                sys.exit()
        else:
            logger.critical('BLASTDB is required for de novo mode. Please provide --unmapped_blastdb.')
            sys.exit()

    if task.remove_host != None:
        if utils.setup_genomes(task.remove_host) == -1:
            logger.error('Host genome not found. Exiting pipeline.')
            sys.exit()
    
    if task.remove_impurities != None:
        if not Path(task.remove_impurities).is_file():
            logger.error('Impurities source file not found. Exiting pipeline.')
            sys.exit()

    logger.info('Checking reads files.')
    if check_reads_file(task) != -1:
        task.id = "%s_%s" % (task.name, time.strftime(
            "%Y%m%d%H%M%S", time.localtime()))
        logger.info('Creating new task %s.' % task.id)
        Path.mkdir(task.path.joinpath(task.id), parents=True)
        logger.info('Starting pipeline.')
        utils.write_log_file(
            task.path.joinpath(task.id),
            'Starting pipeline.'
        )


        db = db_manager.VIVADatabase()
        start_date = time.strftime("%Y-%m-%d %H:%M", time.localtime())
        db.create_task(
            task_id=task.id, 
            task_name=task.name, 
            start_date=start_date, 
            preset_id=getattr(task, 'preset_id', None),
            task_note=task.task_note,
            product=task.sample_product_name,
            lot=task.sample_product_lot,
            seq_date=task.sample_sequencing_date
        )

        try:
            # main pipeline
            reads_preprocess.run(task)
            reference_prepare.run(task)
            impurities_prefilter.run(task)
            reads_alignment.run(task)
            unmapped_analysis.run(task)
            variant_calling.run(task)
            logger.info('Pipeline finished.')
            utils.write_log_file(
                task.path.joinpath(task.id),
                'Pipeline finished.'
            )

            # report generator
            summary_generator.run(task)
            report_generator.run(task)
            
            # DB finalizer triggers in summary_generator, we just mark completed here after safe wrap.
            db.update_task_status(task.id, 'Completed')

            # Auto cleanup
            if task.auto_cleanup == 'True':
                logger.info('Starting auto cleanup.')
                utils.write_log_file(
                    task.path.joinpath(task.id),
                    'Starting auto cleanup.'
                )
                cleanup.cleanup_task(task.path.joinpath(task.id), force=True)
                logger.info('Auto cleanup finished.')
                utils.write_log_file(
                    task.path.joinpath(task.id),
                    'Auto cleanup finished.'
                )
            
        except Exception as e:
            logger.error(f'Pipeline error: {e}')
            db.update_task_status(task.id, 'Failed', error_log=str(e))
            raise e

        # 測試模式下自動執行驗證
        if args.test != None:
            logger.info("Test mode detected, performing auto-verification.")
            expected_json = Path.cwd().joinpath('test_data', 'expected_results.json')
            report, md = e2e_verifier.run_verification(task.id, task.path, expected_json)
            # 將驗證結果寫入任務目錄
            with open(task.path.joinpath(task.id, 'verification_report.md'), 'w') as f:
                f.write(md)
            if not report['passed']:
                logger.error("Verification failed!")
            else:
                logger.info("Verification passed.")

        # 權限修復
        fix_permissions(task.path.joinpath(task.id))
        if task.blastdb_path:
            fix_permissions(task.blastdb_path)

        return task.id

    else:
        logger.error('Reads not found. Exiting pipeline.')
        sys.exit()


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])

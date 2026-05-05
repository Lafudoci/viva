import logging
import subprocess
import sys
import os
import shutil
from decimal import Decimal
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_de_novo(task):
    logger.info('Runing De Novo assembly of unmapped reads.')
    unmapped_assembly_cwd = task.path.joinpath(task.id, 'unmapped_analysis')
    Path.mkdir(unmapped_assembly_cwd, parents=True, exist_ok=True)
    bwa_aligner_cwd = task.path.joinpath(task.id, 'alignment', 'bwa')
    output_folder_name = '%s_unmapped_spades_%s'%(task.id, task.unmapped_spades_mode)
    # ONLY apply to the first bwa ref alignment.
    r1 = str(task.path.joinpath(bwa_aligner_cwd, task.id+'_ref_1_unmapped_R1.fastq.gz'))
    r2 = str(task.path.joinpath(bwa_aligner_cwd, task.id+'_ref_1_unmapped_R2.fastq.gz'))
    
    if getattr(task, 'unmapped_bbnorm', 'False') == 'True':
        logger.info('Running bbnorm.sh normalization.')
        norm_r1 = str(task.path.joinpath(unmapped_assembly_cwd, task.id+'_ref_1_unmapped_norm_R1.fastq.gz'))
        norm_r2 = str(task.path.joinpath(unmapped_assembly_cwd, task.id+'_ref_1_unmapped_norm_R2.fastq.gz'))
        target = getattr(task, 'unmapped_bbnorm_target', '30')
        min_cov = getattr(task, 'unmapped_bbnorm_min', '2')
        bbnorm_cmd = [
            'bbnorm.sh',
            'in=' + r1,
            'in2=' + r2,
            'out=' + norm_r1,
            'out2=' + norm_r2,
            'target=' + str(target),
            'min=' + str(min_cov)
        ]
        logger.info('CMD: '+' '.join(bbnorm_cmd))
        utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(bbnorm_cmd))
        cmd_run = subprocess.run(bbnorm_cmd, cwd=unmapped_assembly_cwd, capture_output=True)
        print(cmd_run.stdout.decode(encoding='utf-8'))
        print(cmd_run.stderr.decode(encoding='utf-8'))
        
        r1 = norm_r1
        r2 = norm_r2

    assemble_cmd = [
        'spades.py',
        '-t', task.threads,
        '-m', task.spades_mem,
        '--' + task.unmapped_spades_mode,
        '--phred-offset', '33',
        '-1', r1,
        '-2', r2,
        '-o', output_folder_name
    ]
    logger.info('CMD: '+' '.join(assemble_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(assemble_cmd))
    cmd_run = subprocess.run(assemble_cmd, cwd=unmapped_assembly_cwd, capture_output=True)
    stdout_msg = cmd_run.stdout.decode(encoding='utf-8')
    stderr_msg = cmd_run.stderr.decode(encoding='utf-8')
    print(stdout_msg)
    print(stderr_msg)

    if cmd_run.returncode != 0:
        # Retry Phase 1: Half threads (to handle Race Conditions and moderate OOM)
        retry_threads = str(max(1, int(task.threads) // 2))
        logger.warning(f'SPAdes failed (Exit code: {cmd_run.returncode}). Retry Phase 1: reducing threads to {retry_threads}...')
        
        # Preserve log before deletion
        failed_log = Path(unmapped_assembly_cwd, output_folder_name, 'spades.log')
        if failed_log.is_file():
            shutil.copy(failed_log, Path(unmapped_assembly_cwd, f'{output_folder_name}_failed_1.log'))
        
        shutil.rmtree(Path(unmapped_assembly_cwd, output_folder_name), ignore_errors=True)
        
        retry_cmd = assemble_cmd[:]
        try:
            t_idx = retry_cmd.index('-t')
            retry_cmd[t_idx+1] = retry_threads
        except (ValueError, IndexError):
            pass
        
        logger.info('Retry Phase 1 CMD: '+' '.join(retry_cmd))
        utils.write_log_file(task.path.joinpath(task.id), 'Retry Phase 1 CMD: '+' '.join(retry_cmd))
        cmd_run = subprocess.run(retry_cmd, cwd=unmapped_assembly_cwd, capture_output=True)
        stdout_msg = cmd_run.stdout.decode(encoding='utf-8')
        stderr_msg = cmd_run.stderr.decode(encoding='utf-8')
        print(stdout_msg)
        print(stderr_msg)
        
        # Final Retry Phase 2: --only-assembler (to bypass all Error Correction bugs/limits)
        if cmd_run.returncode != 0:
            logger.warning(f'SPAdes failed again (Exit code: {cmd_run.returncode}). Final Retry Phase 2: using --only-assembler...')
            
            # Preserve log before second deletion
            failed_log = Path(unmapped_assembly_cwd, output_folder_name, 'spades.log')
            if failed_log.is_file():
                shutil.copy(failed_log, Path(unmapped_assembly_cwd, f'{output_folder_name}_failed_2.log'))
            
            shutil.rmtree(Path(unmapped_assembly_cwd, output_folder_name), ignore_errors=True)
            
            final_cmd = assemble_cmd + ['--only-assembler']
            logger.info('Retry Phase 2 CMD: '+' '.join(final_cmd))
            utils.write_log_file(task.path.joinpath(task.id), 'Retry Phase 2 CMD: '+' '.join(final_cmd))
            cmd_run = subprocess.run(final_cmd, cwd=unmapped_assembly_cwd, capture_output=True)
            print(cmd_run.stdout.decode(encoding='utf-8'))
            print(cmd_run.stderr.decode(encoding='utf-8'))

    # check if assemble result exists
    contigs_path = Path(unmapped_assembly_cwd, output_folder_name, 'contigs.fasta')
    if contigs_path.is_file() != True:
        logger.error('Assemble file not found.')
        return -1


def blast_assembled(task):
    logger.info('BLASTing unmapped reads assembled.')
    assembled_cwd = task.path.joinpath(task.id, 'unmapped_analysis', '%s_unmapped_spades_%s'%(task.id, task.unmapped_spades_mode))
    qeury_filename = 'contigs.fasta'
    if task.unmapped_blastdb_extra_list != None:
        blastdbs = [task.unmapped_blastdb] + task.unmapped_blastdb_extra_list.split()
    else:
        blastdbs = [task.unmapped_blastdb]
    
    highly_match_result_dict = {}
    m_env = os.environ.copy()
    m_env['BLASTDB'] = task.blastdb_path
    for db in blastdbs:
        blast_result_filename = '%s_spades_%s_%s.tsv'%(task.id, task.unmapped_spades_mode, db)
        blast_cmd = [
            'blastn',
            '-db',
            db,
            '-query',
            qeury_filename,
            '-out',
            blast_result_filename,
            '-outfmt',
            '6 qseqid sacc pident qlen length evalue stitle bitscore qcovs qstart qend sstart send',
            '-num_threads',
            task.threads,
            '-max_target_seqs',
            '100',
            '-max_hsps',
            '1',
            '-evalue',
            '1e-6',
        ]
        logger.info('CMD: '+' '.join(blast_cmd))
        utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(blast_cmd))
        cmd_run = subprocess.run(blast_cmd, cwd=assembled_cwd, capture_output=True, env=m_env)
        print(cmd_run.stdout.decode(encoding='utf-8'))
        print(cmd_run.stderr.decode(encoding='utf-8'))
    
        # filter highly matched hits and add annotation from RVDB
        logger.info('Filter highly matched hits')
        blast_result_path = assembled_cwd.joinpath(blast_result_filename)
        annotation_dict = utils.load_rvdb_anno_tab(task.rvdb_anno_path)
        filtered_hits_list = []
        hit = []
        with open(blast_result_path, 'r') as f:
            for line in f.readlines():
                hit = line.strip().split('\t')
                if blast_hits_significant_filter(task, hit):
                    filtered_hit = blast_hits_string_formater(db, hit)
                    filtered_hit = blast_hits_anno_finder(db, filtered_hit, annotation_dict)
                    filtered_hits_list.append(filtered_hit)

        # Get the date information for the current BLAST database
        db_date = "Unknown"
        try:
            info_cmd = f"blastdbcmd -db {db} -info | awk '/^Date:/ {{print $2, $3, $4}}'"

            cmd_result = subprocess.run(
                info_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                env=m_env
            )
            db_date = cmd_result.stdout.strip()
            if not db_date:
                db_date = "Date not found"
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to get info for BLAST DB '{db}'. Error: {e.stderr.strip()}")
            db_date = "Error fetching date"
        except Exception as e:
            logger.error(f"An unexpected error occurred while fetching DB info for '{db}': {e}")
            db_date = "Error"
    
        highly_match_result_dict[db] = {
            'BLASTdb_name': db,
            'BLASTdb_date': db_date,
            'highly_matched_result': blast_hits_max1_bitscore_filter(task, filtered_hits_list)
        }

    return highly_match_result_dict

def build_unmapped_json(task, result_dict):
    unmapped_analysis_json_path = task.path.joinpath(task.id, 'unmapped_analysis', 'unmapped_analysis.json')
    utils.build_json_file(unmapped_analysis_json_path, result_dict)


def blast_hits_max1_bitscore_filter(task, hits_list):
    filtered_dict = {}
    filtered_list = []
    for hit in hits_list:
        if not hit['qseqid'] in filtered_dict:
            filtered_dict[hit['qseqid']] = {
                'current_best_score': Decimal('0'),
                'hit': {}
            }
        if Decimal(hit['bitscore']) > Decimal(filtered_dict[hit['qseqid']]['current_best_score']):
            filtered_dict[hit['qseqid']]['hit'] = hit.copy()
            filtered_dict[hit['qseqid']]['current_best_score'] = str(hit['bitscore'])
    for v in filtered_dict.values():
        filtered_list.append(v['hit'])
    return filtered_list


def blast_hits_significant_filter(task, hit):
    if Decimal(hit[2]) < Decimal(task.unmapped_ident_filter):
        return False
    if Decimal(hit[4]) < Decimal(task.unmapped_len_filter):
        return False
    else:
        return True


def blast_hits_string_formater(db, hit):
    if db.startswith(("U-RVDB","C-RVDB")):
        hit_split_list = hit[6].split('|')
        clean_sacc = hit_split_list[2]
        clean_stitle = hit_split_list[3]
        if len(hit_split_list) >= 5:
            clean_stitle_org = hit_split_list[4]
        else:
            clean_stitle_org = ''
    else:
        clean_sacc = hit[1].replace('|', ' ')
        clean_stitle = hit[6].replace('|', ' ')
        clean_stitle_org = ''
    return {
        'qseqid': hit[0],
        'sacc': hit[1],
        'pident': hit[2],
        'qlen': hit[3],
        'length': hit[4],
        'evalue': hit[5],
        'stitle': hit[6],
        'bitscore': hit[7],
        'qcovs': hit[8],
        'qstart': hit[9],
        'qend': hit[10],
        'sstart': hit[11],
        'send': hit[12],
        'clean_sacc': clean_sacc,
        'clean_stitle': clean_stitle,
        'clean_stitle_org': clean_stitle_org
    }


def blast_hits_anno_finder(db, hit, annotation):
    hit['anno'] = ''
    if db.startswith(("U-RVDB","C-RVDB")):
        clean_sacc = hit.get('clean_sacc')
        sstart = hit.get('sstart')
        send = hit.get('send')
        if clean_sacc and sstart is not None and send is not None:
            if clean_sacc in annotation:
                logger.info('Hit matches RVDB annotation.')
                for anno_entry in annotation[clean_sacc]:
                    anno_start = anno_entry['start']
                    anno_end = anno_entry['end']
                    category = anno_entry['category']
                    # swap alignment direction to forward if it's not
                    if sstart > send:
                        sstart, send = send, sstart
                    # Check if alignment overlaps annotation range
                    overlap = not(int(send) < int(anno_start) or int(anno_end) < int(sstart))
                    if overlap:
                        hit['anno'] = category
        else:
            logger.warning('Hit result is not completed, skipping.')
    return hit


def run(task):
    if task.unmapped_assemble == 'True':
        contigs = run_de_novo(task)
        if contigs != -1:
            if task.unmapped_blastdb != None:
                unmapped_results_list = blast_assembled(task)
                build_unmapped_json(task, unmapped_results_list)
            else:
                logger.warning('unmapped_blastdb not set, skipping blast.')
                unmapped_results_list = {}
                build_unmapped_json(task, unmapped_results_list)
        else:
            logger.warning('Contigs not found, skipping blast.')
            unmapped_results_list = {}
            build_unmapped_json(task, unmapped_results_list)
    else:
        logger.warning('unmapped_assemble not set, skipping assemble and blast.')
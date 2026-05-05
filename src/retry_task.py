"""
retry_task.py
職責：從既有 task 目錄恢復執行參數，重新執行缺失的分析步驟，並固定重新產製報告。
"""
import configparser
import json
import logging
import re
import sys
from pathlib import Path

import db_manager
import impurities_prefilter
import reads_alignment
import reference_prepare
import report_generator
import summary_generator
import unmapped_analysis
import utils
import variant_calling

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


# ─────────────────────────────────────────────────────────────────────────────
# 參數恢復輔助函式（三層向前相容策略）
# ─────────────────────────────────────────────────────────────────────────────

def _load_log_lines(log_file):
    """讀取 log.txt 全部行。"""
    with open(log_file, 'r', encoding='utf-8') as f:
        return f.readlines()


def _parse_params_from_log(log_lines):
    """層一：從 [PARAMS] JSON 行取得完整參數 dict。"""
    for line in log_lines:
        parts = line.strip().split('\t', 1)
        if len(parts) == 2 and parts[1].startswith('[PARAMS] '):
            try:
                params = json.loads(parts[1][len('[PARAMS] '):])
                logger.info('Retry：成功從 log.txt 讀取 [PARAMS]。')
                return params
            except json.JSONDecodeError as e:
                logger.error('Retry：解析 [PARAMS] JSON 失敗：%s' % e)
                sys.exit(1)
    return None


def _parse_params_from_summary(summary_path):
    """層二：從 summary.json 恢復可取得的基本參數。"""
    if not summary_path.is_file():
        logger.info('Retry：找不到 summary.json，跳過層二。')
        return {}
    logger.info('Retry：找不到 [PARAMS]，嘗試從 summary.json 恢復參數。')
    try:
        s = utils.load_json_file(summary_path)
        rm = s.get('reads_meta', {})
        ref_meta = s.get('ref_meta_dict', {})
        sample_meta = rm.get('sample_meta', {})
        reads_file_meta = rm.get('reads_file_meta', {})
        unmapped_dbs = [k for k in s.get('unmapped_analysis', {}).keys() if k != 'N/A']
        params = {
            'ex_r1':  reads_file_meta.get('file_name', {}).get('r1'),
            'ex_r2':  reads_file_meta.get('file_name', {}).get('r2'),
            'ref':    ref_meta.get('origin_file_path') if ref_meta.get('ref_from_user') == 'Yes' else None,
            'spades_mode': ref_meta.get('spades_mode', 'metaviral'),
            'task_note': s.get('task_note'),
            'sample_product_name':    sample_meta.get('sample_product_name'),
            'sample_product_lot':     sample_meta.get('sample_product_lot'),
            'sample_sequencing_date': sample_meta.get('sample_sequencing_date'),
            'sample_note':            sample_meta.get('sample_note'),
            'unmapped_blastdb': unmapped_dbs[0] if unmapped_dbs else None,
        }
        logger.info('Retry：從 summary.json 恢復基本參數完成。')
        return params
    except Exception as e:
        logger.warning('Retry：讀取 summary.json 失敗（%s），繼續嘗試 CMD 解析。' % e)
        return {}


def _parse_params_from_cmd(log_lines, params):
    """
    層三：從 log.txt CMD 行正則解析，補充 params 中仍缺少的欄位。
    僅在對應 key 尚未被層一/層二填入時才寫入。
    """
    # 蒐集所有 CMD 行文字
    cmd_lines = [
        l.strip().split('\t', 1)[1][5:]
        for l in log_lines
        if '\t' in l and l.strip().split('\t', 1)[1].startswith('CMD: ')
    ]
    all_cmds = ' '.join(cmd_lines)

    def _find(keyword):
        for c in cmd_lines:
            if keyword in c:
                return c
        return ''

    # threads
    if not params.get('threads'):
        m = re.search(r'spades\.py\b.*?\s-t\s+(\d+)', all_cmds) or \
            re.search(r'fastp\b.*?\s-w\s+(\d+)', all_cmds)
        params['threads'] = m.group(1) if m else '6'

    # global_trimming
    if not params.get('global_trimming'):
        m = re.search(r'fastp\b.*?\s-f\s+(\d+)', _find('fastp'))
        params['global_trimming'] = m.group(1) if m else '0'

    # ex_r1 / ex_r2（從 fastp -i/-I）
    if not params.get('ex_r1'):
        m = re.search(r'fastp\b.*?\s-i\s+(\S+)', _find('fastp'))
        params['ex_r1'] = m.group(1) if m else None
    if not params.get('ex_r2'):
        m = re.search(r'fastp\b.*?\s-I\s+(\S+)', _find('fastp'))
        params['ex_r2'] = m.group(1) if m else None

    # remove_host（從 bowtie2 -x /app/genomes/<name> 反推）
    if not params.get('remove_host'):
        m = re.search(r'-x\s+\S+/genomes/(\S+)', _find('host_mapped.sam'))
        if m:
            genome_map = {
                'GRCh38.p14': 'human', 'dog10k': 'dog',
                'vero': 'vero', 'grcg6a': 'chicken', 'mmul_10': 'rhesus_monkey'
            }
            params['remove_host'] = genome_map.get(m.group(1), m.group(1))

    # alns（檢查 ref alignment CMD）
    if not params.get('alns'):
        detected = []
        if any('bowtie2' in c and 'ref_1' in c and 'sorted.bam' not in c for c in cmd_lines):
            detected.append('bowtie2')
        if any('bwa mem' in c and 'ref_' in c for c in cmd_lines):
            detected.append('bwa')
        params['alns'] = ','.join(detected) if detected else 'bowtie2,bwa'

    # spades_mem
    if not params.get('spades_mem'):
        m = re.search(r'spades\.py\b.*?\s-m\s+(\d+)', _find('spades.py'))
        params['spades_mem'] = m.group(1) if m else '22'

    # unmapped_spades_mode（從輸出目錄名稱或 -- 參數）
    if not params.get('unmapped_spades_mode'):
        uc = _find('unmapped_spades')
        m = re.search(r'unmapped_spades_(\w+)', uc) or \
            re.search(r'spades\.py\b.*?--(metaviral|meta|rnaviral|corona)\b', uc)
        params['unmapped_spades_mode'] = m.group(1) if m else 'meta'

    # unmapped_blastdb
    if not params.get('unmapped_blastdb'):
        m = re.search(r'blastn\b.*?\s-db\s+(\S+)', _find('blastn'))
        params['unmapped_blastdb'] = m.group(1) if m else None

    # unmapped_assemble
    if not params.get('unmapped_assemble'):
        params['unmapped_assemble'] = 'True' if any('unmapped_spades' in c for c in cmd_lines) else 'False'

    # 無法從 log 推算的參數，填預設值並警告
    _defaults = {
        'vc_threshold': '0.7', 'min_vc_score': '1',
        'unmapped_bbnorm': 'False', 'unmapped_bbnorm_target': '30', 'unmapped_bbnorm_min': '2',
        'unmapped_len_filter': '500', 'unmapped_ident_filter': '95',
        'blastdb_path': None, 'rvdb_anno_path': None,
        'unmapped_blastdb_extra_list': None, 'remove_impurities': None,
        'preset_path': None,
    }
    # 使用 is None 判斷，避免 'False'、'0' 等 falsy 值被預設值錯誤覆寫
    missing_keys = [k for k, v in _defaults.items() if params.get(k) is None]
    for k in missing_keys:
        params[k] = _defaults[k]

    if missing_keys:
        logger.warning(
            'Retry（向前相容模式）：%s 等欄位無法從 log 推算，使用預設值。'
            '若有需要請手動確認報告結果。' % ', '.join(missing_keys)
        )
    return params


def recover_params(task_dir):
    """
    依三層優先順序恢復參數：
    1. log.txt [PARAMS] JSON
    2. summary.json
    3. log.txt CMD 行正則解析
    回傳最終的 params dict。
    """
    log_file = task_dir.joinpath('log.txt')
    if not log_file.is_file():
        logger.error('Retry 失敗：找不到 log.txt：%s' % log_file)
        sys.exit(1)

    log_lines = _load_log_lines(log_file)

    # 層一：從 [PARAMS] JSON 完整恢復
    params = _parse_params_from_log(log_lines)

    if params is None:
        # 層二（僅在層一失敗時）：從 summary.json 恢復基本欄位
        params = _parse_params_from_summary(task_dir.joinpath(task_dir.name + '_summary.json'))
        # 層三（僅在層一失敗時）：從 CMD 行正則解析補充缺漏欄位，並發出相容警告
        params = _parse_params_from_cmd(log_lines, params)
    else:
        logger.info('Retry：以 [PARAMS] JSON 完整恢復參數，跳過向前相容解析。')

    return params


# ─────────────────────────────────────────────────────────────────────────────
# Task 物件重建
# ─────────────────────────────────────────────────────────────────────────────

def _build_task_from_params(task_dir, params):
    """從 params dict 建立並回傳 task 物件。需從 new_task 匯入 Task 類別。"""
    from new_task import Task

    task = Task()
    task.id   = task_dir.name
    task.path = task_dir.parent
    task.name = '_'.join(task.id.split('_')[:-1])

    task.ex_r1 = params.get('ex_r1')
    task.ex_r2 = params.get('ex_r2')
    task.ref   = params.get('ref')
    task.threads              = params.get('threads', '6')
    task.alns                 = params.get('alns', 'bowtie2,bwa').split(',')
    task.global_trimming      = params.get('global_trimming', '0')
    task.remove_host          = params.get('remove_host')
    task.remove_impurities    = params.get('remove_impurities')
    task.spades_mem           = params.get('spades_mem', '22')
    task.spades_mode          = params.get('spades_mode', 'metaviral')
    task.unmapped_spades_mode = params.get('unmapped_spades_mode', 'meta')
    task.unmapped_bbnorm        = params.get('unmapped_bbnorm', 'False')
    task.unmapped_bbnorm_target = params.get('unmapped_bbnorm_target', '30')
    task.unmapped_bbnorm_min    = params.get('unmapped_bbnorm_min', '2')
    task.vc_threshold           = params.get('vc_threshold', '0.7')
    task.min_vc_score           = params.get('min_vc_score', '1')
    task.blastdb_path           = params.get('blastdb_path')
    task.rvdb_anno_path         = params.get('rvdb_anno_path')
    task.unmapped_assemble      = params.get('unmapped_assemble', 'True')
    task.unmapped_blastdb       = params.get('unmapped_blastdb')
    task.unmapped_blastdb_extra_list = params.get('unmapped_blastdb_extra_list')
    task.unmapped_len_filter    = params.get('unmapped_len_filter', '500')
    task.unmapped_ident_filter  = params.get('unmapped_ident_filter', '95')
    task.preset_path            = params.get('preset_path')
    task.task_note              = params.get('task_note')
    task.sample_product_name    = params.get('sample_product_name')
    task.sample_product_lot     = params.get('sample_product_lot')
    task.sample_sequencing_date = params.get('sample_sequencing_date')
    task.sample_note            = params.get('sample_note')

    # 狀態欄位預設值
    task.ref_num = 0
    task.impurities_prefilter_num = 0
    task.total_reads_after_fastp  = 0
    task.with_ref = task.ref not in (None, 'None')

    # 恢復 preset 版本資訊（若 preset 仍可存取）
    if task.preset_path and task.preset_path not in (None, 'None'):
        try:
            config = configparser.ConfigParser(allow_no_value=True)
            config.read(task.preset_path)
            task.preset_id           = config['VERSION']['preset_id']
            task.preset_version      = config['VERSION']['version']
            task.preset_last_rev_date = config['VERSION']['last_rev_date']
            task.preset_author       = config['VERSION']['author']
            task.preset_note         = config['VERSION']['note']
        except Exception:
            task.preset_path = None
    else:
        task.preset_path = None

    task.conda_pkgs = [
        'conda', 'python', 'perl',
        'fastp', 'samtools', 'bcftools', 'htslib',
        'bowtie2', 'bwa',
        'varscan', 'lofreq',
        'spades', 'blast', 'bbmap'
    ]
    return task


# ─────────────────────────────────────────────────────────────────────────────
# 主入口
# ─────────────────────────────────────────────────────────────────────────────

def run(task_dir_path):
    """
    retry 模式主入口。
    於原 task 目錄內重新執行缺失的分析步驟，固定重新產製報告。
    """
    task_dir = Path(task_dir_path).resolve()
    if not task_dir.is_dir():
        logger.error('Retry 失敗：目錄不存在：%s' % task_dir)
        sys.exit(1)

    # ── 恢復參數並建立 task 物件 ──────────────────────────────────────────
    params = recover_params(task_dir)
    task   = _build_task_from_params(task_dir, params)

    # ── 從 fastp.json 恢復 total_reads_after_fastp ───────────────────────
    fastp_json = task.path.joinpath(task.id, 'reads', 'fastp.json')
    if fastp_json.is_file():
        logger.info('[RETRY] reads_preprocess：已找到 fastp.json，跳過此步驟。')
        try:
            task.total_reads_after_fastp = summary_generator.fastp_parser(task)['after_total_reads']
        except Exception as e:
            logger.warning('[RETRY] 讀取 fastp.json 失敗：%s' % e)
    else:
        logger.error('[RETRY] reads_preprocess：未找到 fastp.json，'
                     '請先執行完整的 reads_preprocess 步驟。')
        sys.exit(1)

    # ── 更新資料庫狀態 ───────────────────────────────────────────────────
    logger.info('Retry：開始對 task %s 執行 retry。' % task.id)
    utils.write_log_file(task.path.joinpath(task.id), 'Retry 開始。')

    db = db_manager.VIVADatabase()
    try:
        db.update_task_status(task.id, 'Retrying')
    except Exception as e:
        logger.warning('[RETRY] 更新資料庫狀態失敗（可忽略）：%s' % e)

    # ── 執行各分析步驟（含跳過邏輯） ─────────────────────────────────────
    try:
        reference_prepare.run(task,    is_retry=True)
        impurities_prefilter.run(task, is_retry=True)
        reads_alignment.run(task,      is_retry=True)
        unmapped_analysis.run(task,    is_retry=True)
        variant_calling.run(task,      is_retry=True)

        logger.info('[RETRY] 分析步驟完成。')
        utils.write_log_file(task.path.joinpath(task.id), 'Retry 分析步驟完成。')

        # 固定重新產製報告
        logger.info('[RETRY] 重新產製 summary 與 report。')
        summary_generator.run(task)
        report_generator.run(task)

        try:
            db.update_task_status(task.id, 'Completed')
        except Exception as e:
            logger.warning('[RETRY] 更新資料庫完成狀態失敗（可忽略）：%s' % e)
        logger.info('Retry 完成。')

    except Exception as e:
        logger.error('Retry pipeline 錯誤：%s' % e)
        try:
            db.update_task_status(task.id, 'Failed', error_log=str(e))
        except Exception:
            pass
        raise e

    return task.id

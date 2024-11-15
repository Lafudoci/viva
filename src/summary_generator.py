import json
import logging
import subprocess
from datetime import datetime, timezone, timedelta
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def report_summary(task):
    logger.info('Generating summary.')
    s = {}
    log_abs = log_parser(task)
    s['task_name'] = task.id.split('_')[:-1][0]
    s['task_id'] = task.id
    s['task_note'] = task.task_note
    start_t = datetime.utcfromtimestamp(int(log_abs['start_timestamp']))+ timedelta(hours=8)
    s['start_date'] = start_t.strftime('%Y-%m-%d %H:%M:%S UTC+8')
    finish_t = datetime.utcfromtimestamp(int(log_abs['finish_timestamp']))+ timedelta(hours=8)
    s['finish_date'] = finish_t.strftime('%Y-%m-%d %H:%M:%S UTC+8')
    if task.preset_path != None:
        s['preset'] = {
            'preset_id': task.preset_id,
            'version': task.preset_version,
            'last_rev_date': task.preset_last_rev_date,
            'author': task.preset_author,
            'note': task.preset_note
        }
    s['log_dict'] = log_abs['log_dict']
    if task.sample_product_name in (None,''):
        task.sample_product_name = 'N/A'
    if task.sample_product_lot in (None,''):
        task.sample_product_lot = 'N/A'
    if task.sample_sequencing_date in (None,''):
        task.sample_sequencing_date = 'N/A'
    if task.sample_note in (None,''):
        task.sample_note = 'N/A'
    sample_meta = {
        'sample_product_name': task.sample_product_name,
        'sample_product_lot': task.sample_product_lot,
        'sample_sequencing_date': task.sample_sequencing_date,
        'sample_note': task.sample_note
    }
    s['reads_meta'] = {
        'sample_meta': sample_meta,
        'reads_file_meta': single_meta_parser(task, 'reads', 'reads_meta.json')
    }
    s['fastp_abs'] = fastp_parser(task)
    if task.remove_host != None:
        dehost_meta = single_meta_parser(task, 'reads', 'dehost_meta.json')
        s['remove_genome'] = {'genome': dehost_meta['genome'], 'mapped_reads': dehost_meta['mapped_reads'], 'remove_percentage': dehost_meta['remove_percentage']}
    else:
        s['remove_genome'] = {'genome': 'N/A', 'mapped_reads': 'N/A', 'remove_percentage': 'N/A'}
    s['ref_meta_dict'] = single_meta_parser(task, 'reference', task.id + '_ref.json').copy()
    if task.with_ref == False:
        best_hit = single_meta_parser(task, 'assembly', 'best_hit.json')
        s['ref_source_contig_name'] = best_hit['qseqid']
        s['ref_source_contig_pident'] = best_hit['pident']
    else:
        s['ref_source_contig_name'] = 'N/A'
        s['ref_source_contig_pident'] = 'N/A'
    if task.remove_impurities != None:
        s['impurit_filter_meta'] = single_meta_parser(task, 'impurities_prefilter', 'impurities_prefilter_meta.json')
        s['impurit_filter_results'] = single_meta_parser(task, 'impurities_prefilter', 'impurities_remove.json')
    else:
        s['impurit_filter_meta'] = {}
        s['impurit_filter_results'] = {'1':{
            'bt2':{'mapped_reads':'N/A','remove_percentage':'N/A'},
            'bwa':{'mapped_reads':'N/A','remove_percentage':'N/A'}
            }
        }
    s['aln'] = single_meta_parser(task, 'alignment', 'flagstat.json')
    s['cov'] = single_meta_parser(task, 'alignment', 'coverage_stat.json')
    s['vc'] = vc_parser(task)
    s['draft_meta'] = single_meta_parser(task, 'draft_genome', task.id + '_draft_summary.json').copy()
    if task.unmapped_blastdb != None:
        s['unmapped_analysis'] = single_meta_parser(task, 'unmapped_analysis', 'unmapped_analysis.json')
    else:
        s['unmapped_analysis'] = {'spades_mode':task.unmapped_spades_mode, 'BLASTdb_name':'N/A','highly_matched_result':[]}
    s['version'] = tool_version_caller(task)
    utils.build_json_file(
        task.path.joinpath(task.id, task.id + '_summary.json'),
        s
    )


def single_meta_parser(task, folder_name, file_name):
    meta_json_path = task.path.joinpath(
        task.id, folder_name, file_name
    )
    meta_dict = utils.load_json_file(meta_json_path)
    return meta_dict


def fastp_parser(task):
    fastp_json_path = task.path.joinpath(task.id, 'reads', 'fastp.json')
    fastp_dict = utils.load_json_file(fastp_json_path)
    before_total_reads = fastp_dict['summary']['before_filtering']['total_reads']
    before_total_bases = fastp_dict['summary']['before_filtering']['total_bases']
    before_total_q30 = fastp_dict['summary']['before_filtering']['q30_rate']
    before_r1_length = fastp_dict['summary']['before_filtering']['read1_mean_length']
    before_r2_length = fastp_dict['summary']['before_filtering']['read2_mean_length']
    after_total_reads = fastp_dict['summary']['after_filtering']['total_reads']
    after_total_bases = fastp_dict['summary']['after_filtering']['total_bases']
    after_total_q30 = fastp_dict['summary']['after_filtering']['q30_rate']
    after_r1_length = fastp_dict['summary']['after_filtering']['read1_mean_length']
    after_r2_length = fastp_dict['summary']['after_filtering']['read2_mean_length']
    duplication_rate = fastp_dict['duplication']['rate']
    fastp_abs = {
        'before_total_reads': before_total_reads,
        'before_total_bases': before_total_bases,
        'before_total_q30': before_total_q30,
        'before_r1_length': before_r1_length,
        'before_r2_length': before_r2_length,
        'after_total_reads': after_total_reads,
        'after_total_bases': after_total_bases,
        'after_total_q30': after_total_q30,
        'after_r1_length': after_r1_length,
        'after_r2_length': after_r2_length,
        'duplication_rate': duplication_rate
    }
    return fastp_abs


def log_parser(task):
    log_path = task.path.joinpath(task.id)
    log_list = utils.load_log_file(log_path)
    log_dict = {}
    start_stamp = log_list[0].split('\t')[0]
    end_stamp = log_list[-1].split('\t')[0]
    i = 1
    for log in log_list:
        log_dict[i] = {}
        log_dict[i]['time_stamp'] = int(log.split('\t')[0])
        log_dict[i]['string'] = log.split('\t')[1].strip()
        i += 1
    # duration calculation
    last_time_stamp = start_stamp
    copy_dict = log_dict.copy()
    for order, entry in copy_dict.items():
        log_dict[order]['duration'] = 0
        previous_log_duration = int(entry['time_stamp']) - int(last_time_stamp)
        if order > 1:
            log_dict[order-1]['duration'] = previous_log_duration
        last_time_stamp = int(entry['time_stamp'])
    log_abs = {
        'start_timestamp': start_stamp,
        'finish_timestamp': end_stamp,
        'log_dict': log_dict
    }
    return log_abs


def vc_parser(task):
    vc_json_path = task.path.joinpath(
        task.id, task.id + '_vc_summary.json')
    vc_dict = utils.load_json_file(vc_json_path)
    return vc_dict


def tool_version_caller(task):
    version_dict = utils.conda_pkg_versions(task.conda_pkgs)
    viva_commit_cmd = ['git', 'describe', '--always']
    version_dict['viva'] = subprocess.run(
        viva_commit_cmd,
        capture_output=True).stdout.decode(encoding='utf-8').strip()
    return version_dict
    

def run(task):
    report_summary(task)


if __name__ == "__main__":
    from new_task import Task
    task = Task()
    task.path = Path.cwd()
    task.name = 'TFDA-SARSCOV2-210112'
    task.id = 'TFDA-SARSCOV2-210112_202102010507'
    run(task)
    # print(tool_version_caller())

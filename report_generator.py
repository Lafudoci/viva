import json
import logging
from datetime import datetime
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def report_summary(task_name, base_path):
    s = {}
    log_abs = log_parser(task_name, base_path)
    s['start_date'] = datetime.fromtimestamp(
        int(log_abs['start_timestamp'])).strftime('%Y-%m-%d %H:%M:%S')
    s['finish_date'] = datetime.fromtimestamp(
        int(log_abs['finish_timestamp'])).strftime('%Y-%m-%d %H:%M:%S')
    s['cmd_list'] = log_abs['cmd_list']
    s['reads_meta'] = reads_meta_parser(task_name, base_path)
    s['fastp_abs'] = fastp_parser(task_name, base_path)
    ref_meta = ref_meta_parser(task_name, base_path)
    s['ref_file_name'] = ref_meta['file_name']
    s['ref_fasta_header'] = ref_meta['fasta_header']
    s['ref_from_user'] = ref_meta['user_provide']

    print(s)


def reads_meta_parser(task_name, base_path):
    meta_json_path = base_path.joinpath(
        task_name, 'reads', 'reads_meta.json'
    )
    meta_dict = utils.load_json_file(meta_json_path)
    return meta_dict


def ref_meta_parser(task_name, base_path):
    meta_json_path = base_path.joinpath(
        task_name, 'reference', task_name + '_ref.json'
    )
    meta_dict = utils.load_json_file(meta_json_path)
    return meta_dict


def fastp_parser(task_name, base_path):
    fastp_json_path = base_path.joinpath(task_name, 'reads', 'fastp.json')
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
    duplication_rate = fastp_dict['summary']['duplication']['rate']
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


def log_parser(task_name, base_path):
    log_path = base_path.joinpath(task_name)
    log_list = utils.load_log_file(log_path)
    cmd_list = []
    for cmd in log_list:
        cmd_list.append(cmd.split('\t')[1][5:])
    log_abs = {
        'start_timestamp': log_list[0].split('\t')[0],
        'finish_timestamp': log_list[-1].split('\t')[0],
        'cmd_list': cmd_list
    }
    return log_abs


def run(task_name, base_path):
    report_summary(task_name, base_path)


if __name__ == "__main__":
    # fastp_parser('test', Path.cwd())
    run('test', Path.cwd())

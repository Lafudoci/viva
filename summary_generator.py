import json
import logging
import subprocess
from datetime import datetime, timezone, timedelta
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def report_summary(task_name, base_path):
    logger.info('Generating summary.')
    s = {}
    log_abs = log_parser(task_name, base_path)
    s['task_name'] = task_name.split('_')[:-1][0]
    s['task_id'] = task_name
    start_t = datetime.utcfromtimestamp(int(log_abs['start_timestamp']))+ timedelta(hours=8)
    s['start_date'] = start_t.strftime('%Y-%m-%d %H:%M:%S UTC+8')
    finish_t = datetime.utcfromtimestamp(int(log_abs['finish_timestamp']))+ timedelta(hours=8)
    s['finish_date'] = finish_t.strftime('%Y-%m-%d %H:%M:%S UTC+8')
    s['cmd_list'] = log_abs['cmd_list']
    s['reads_meta'] = reads_meta_parser(task_name, base_path)
    s['fastp_abs'] = fastp_parser(task_name, base_path)
    ref_meta = ref_meta_parser(task_name, base_path)
    s['ref_file_name'] = ref_meta['file_name']
    s['ref_fasta_header'] = ref_meta['fasta_header']
    s['ref_from_user'] = ref_meta['user_provide']
    s['aln'] = aln_meta_parser(task_name, base_path)
    s['vc'] = vc_parser(task_name, base_path)
    s['version'] = tool_version_caller()
    utils.build_json_file(
        base_path.joinpath(task_name, task_name + '_summary.json'),
        s
    )


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


def aln_meta_parser(task_name, base_path):
    meta_json_path = base_path.joinpath(
        task_name, 'alignment', 'flagstat.json'
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


def log_parser(task_name, base_path):
    log_path = base_path.joinpath(task_name)
    log_list = utils.load_log_file(log_path)
    cmd_list = []
    for log in log_list:
        if log.split('\t')[1].startswith('CMD:'):
            cmd_list.append(log.split('\t')[1][5:])
    log_abs = {
        'start_timestamp': log_list[0].split('\t')[0],
        'finish_timestamp': log_list[-1].split('\t')[0],
        'cmd_list': cmd_list
    }
    return log_abs


def vc_parser(task_name, base_path):
    vc_json_path = base_path.joinpath(
        task_name, task_name + '_vc_summary.json')
    vc_dict = utils.load_json_file(vc_json_path)
    return vc_dict


def tool_version_caller():
    version_dict = {}
    fastp_cmd = ['fastp', '--version']
    samtools_cmd = ['samtools', '--version']
    bowtie2_cmd = ['bowtie2', '--version']
    bwa_cmd = ['bwa']
    lofreq_cmd = ['lofreq', 'version']
    varscan2_cmd = ['java', '-jar', '/home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar']
    version_dict['fastp'] = subprocess.run(
        fastp_cmd, capture_output=True).stderr.decode(encoding='utf-8').split(' ')[1].strip()
    version_dict['samtools'] = subprocess.run(
        samtools_cmd, capture_output=True).stdout.decode(encoding='utf-8').split(' ')[1].split('\n')[0]
    version_dict['bowtie2'] = subprocess.run(
        bowtie2_cmd, capture_output=True).stdout.decode(encoding='utf-8').split(' ')[2].split('\n')[0]
    version_dict['bwa'] = subprocess.run(
        bwa_cmd, capture_output=True).stderr.decode(encoding='utf-8').split(' ')[6].split('\n')[0]
    version_dict['lofreq'] = subprocess.run(
        lofreq_cmd, capture_output=True).stdout.decode(encoding='utf-8').split(' ')[1].split('\n')[0]
    version_dict['varscan2'] = subprocess.run(
        varscan2_cmd, capture_output=True).stderr.decode(encoding='utf-8').split(' ')[1].split('\n')[0][1:]
    return version_dict
    

def run(task_name, base_path):
    report_summary(task_name, base_path)


if __name__ == "__main__":
    # fastp_parser('test', Path.cwd())
    run('TFDA-CMV-210125_202101290325', Path.cwd())
    # tool_version_caller()

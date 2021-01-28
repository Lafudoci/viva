import json
import logging
from datetime import datetime
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def build_md_report(task_name, base_path):
    summary_path = base_path.joinpath(task_name, task_name + '_summary.json')
    s = utils.load_json_file(summary_path)
    headline = '# Virus Variant Calling Report'

    meta_t = '## Meta'
    meta_c = '\n'.join([
        'Task start time : %s'% s['start_date'],
        'Task finish time : %s'% s['finish_date']
    ])

    reads_t = '## Input Reads'
    reads_c = '\n'.join([
        '| Reads | File name |',
        '| ----- | --------- |',
        '| R1 | %s |'%(s['reads_meta']['r1']),
        '| R2 | %s |'%(s['reads_meta']['r2']),
        '\n',
        '| Reads | MD5 hash |',
        '| ----- | --------- |',
        '| R1 | %s |'%(s['reads_meta']['r1']),
        '| R2 | %s |'%(s['reads_meta']['r2'])
    ])

    ref_t = '## Input Reference'
    ref_c = '\n'.join([
        'Reference by user : %s'%(s['ref_from_user']),
        'FASTA Header : %s'%(s['ref_fasta_header']),
        'FASTA file name : %s'%(s['ref_file_name']),
    ])
    
    fastp_t = '## Reads Filter Statistics'
    fastp_f_t = '### Filter'
    fastp_f_c = '\n'.join([
        '|   | Before | After |',
        '| - | ------ | ----- |',
        '| Total reads | %s | %s |'%(s['fastp_abs']['before_total_reads'], s['fastp_abs']['after_total_reads']),
        '| Total base  | %s | %s |'%(s['fastp_abs']['before_total_bases'], s['fastp_abs']['after_total_bases']),
        '| Q30 rate    | %s | %s |'%(s['fastp_abs']['before_total_q30'], s['fastp_abs']['after_total_q30']),
        '| Mean length R1 | %s | %s |'%(s['fastp_abs']['before_r1_length'], s['fastp_abs']['after_r1_length']),
        '| Mean length R2 | %s | %s |'%(s['fastp_abs']['before_r2_length'], s['fastp_abs']['after_r2_length'])
    ])
    fastp_d_t = '### Duplication'
    fastp_d_c = 'Duplication rate : %.02f%%' % float(s['fastp_abs']['duplication_rate']*100)

    aln_t = '## Alignment'
    aln_c = '\n'.join([
        '| Aligner | Overall mapped rate |',
        '| ------- | ------------------- |',
        '| Bowtie2 | %s |'%(s['aln']['mapped_rate']['bowtie2']),
        '| BWA MEM | %s |'%(s['aln']['mapped_rate']['bwa'])
    ])
    vc_t = '## Variant Calling'
    vc_l_t = '### LoFreq'
    vc_l_c = '| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |\n| -------- | --- | --- | ------- | --- |\n'
    for pos in s['vc']['lofreq']:
        ref = s['vc']['lofreq'][pos]['REF']
        for alt, aln_result in s['vc']['lofreq'][pos]['SNV'].items():
            bt2_ft = aln_result.get('bowtie2', {}).get('FILTER', '-')
            bt2_af = aln_result.get('bowtie2', {}).get('FREQ', '-')
            bt2_dp = aln_result.get('bowtie2', {}).get('DP', '-')
            bt2_gq = aln_result.get('bowtie2', {}).get('QUAL', '-')
            bwa_ft = aln_result.get('bwa', {}).get('FILTER', '-')
            bwa_af = aln_result.get('bwa', {}).get('FREQ', '-')
            bwa_dp = aln_result.get('bwa', {}).get('DP', '-')
            bwa_gq = aln_result.get('bwa', {}).get('QUAL', '-')
            vc_l_c += '| %s | %s | %s | %s / %s / %s / %s | %s / %s / %s / %s |\n'%(pos, ref, alt, bt2_ft, bt2_af, bt2_dp, bt2_gq, bwa_ft, bwa_af, bwa_dp, bwa_gq)
    vc_v_t = '### Varscan2'
    vc_v_c = '| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |\n| -------- | --- | --- | ------- | --- |\n'
    for pos in s['vc']['varscan']:
        ref = s['vc']['varscan'][pos]['REF']
        for alt, aln_result in s['vc']['varscan'][pos]['SNV'].items():
            bt2_ft = aln_result.get('bowtie2', {}).get('FILTER', '-')
            bt2_af = aln_result.get('bowtie2', {}).get('FREQ', '-')
            bt2_dp = aln_result.get('bowtie2', {}).get('DP', '-')
            bt2_gq = aln_result.get('bowtie2', {}).get('QUAL', '-')
            bwa_ft = aln_result.get('bwa', {}).get('FILTER', '-')
            bwa_af = aln_result.get('bwa', {}).get('FREQ', '-')
            bwa_dp = aln_result.get('bwa', {}).get('DP', '-')
            bwa_gq = aln_result.get('bwa', {}).get('QUAL', '-')
            vc_v_c += '| %s | %s | %s | %s / %s / %s / %s | %s / %s / %s / %s |\n'%(pos, ref, alt, bt2_ft, bt2_af, bt2_dp, bt2_gq, bwa_ft, bwa_af, bwa_dp, bwa_gq)
    genome_t = '## Draft Genome'
    cmd_t = '## Commands'
    cmd_c = "```\n%s\n```\n"%('\n'.join(s['cmd_list']))
    ver_t = '## Versions'
    ver_c = '| Tool | Version |\n| ---- | ------- |\n'
    for tool, version in s['version'].items():
        ver_c += '| %s | %s |\n'%(tool, version)


    md_str = '\n'.join([
        headline,
        meta_t,
        meta_c,
        reads_t,
        reads_c,
        ref_t,
        ref_c,
        fastp_t,
        fastp_f_t,
        fastp_f_c,
        fastp_d_t,
        fastp_d_c,
        aln_t,
        aln_c,
        vc_t,
        vc_l_t,
        vc_l_c,
        vc_v_t,
        vc_v_c,
        genome_t,
        cmd_t,
        cmd_c,
        ver_t,
        ver_c
    ])

    print(md_str)
    utils.build_text_file(base_path.joinpath(task_name, task_name + '_report.md'), md_str)

def run(task_name, base_path):
    build_md_report(task_name, base_path)



if __name__ == "__main__":
    build_md_report('TFDA-CMV-2020dec_202101280652', Path.cwd())
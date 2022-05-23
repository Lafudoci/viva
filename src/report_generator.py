import json
import logging
from datetime import datetime
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def build_md_report(task):
    summary_path = task.path.joinpath(task.id, task.id + '_summary.json')
    s = utils.load_json_file(summary_path)
    headline = '# Virus Integrity and Variant Analyzer Report'

    meta_t = '## Meta'
    meta_c = '\n\n'.join([
        'Task name : %s'% s['task_name'],
        'Task ID : %s'% s['task_id'],
        'Task start time : %s'% s['start_date'],
        'Task finish time : %s'% s['finish_date']
    ])


    reads_t = '## Input Reads'
    reads_c = '\n'.join([
        '| Reads | File name |',
        '| ----- | --------- |',
        '| R1 | %s |'%(s['reads_meta']['file_name']['r1']),
        '| R2 | %s |'%(s['reads_meta']['file_name']['r2']),
        '\n',
        '| Reads | MD5 hash |',
        '| ----- | --------- |',
        '| R1 | %s |'%(s['reads_meta']['md5']['r1']),
        '| R2 | %s |'%(s['reads_meta']['md5']['r2'])
    ])


    ref_t = '## Input Reference'
    ref_c_m = '\n\n'.join([
        'Reference by user : %s'%(s['ref_meta_dict']['ref_from_user']),
        'Origin file path : %s'%(s['ref_meta_dict']['origin_file_path']),
        'Numeber of reference : %s'%(s['ref_meta_dict']['ref_num']),
        'SPAdes assembly mode : %s'%(s['ref_meta_dict']['spades_mode']),
        'Best contig name : %s'%(s['ref_source_contig_name']),
        'Best contig identity : %s %%'%(s['ref_source_contig_pident']),
    ])
    ref_c_t = '\n'.join([
        '| Order | Header | Length |',
        '| ----- | ------ | ------ |'
    ])
    for ref_order in range(1, task.ref_num+1):
        header = s['ref_meta_dict']['seq_meta'][str(ref_order)]['fasta_header_escape']
        length = s['ref_meta_dict']['seq_meta'][str(ref_order)]['seq_length']
        ref_c_t += '\n| %d | %s | %s |'%(
            ref_order,
            header,
            length
            )

    
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
    

    dehost_t = '## Host Genome Removal'
    dehost_i_g = 'Remove genome: %s'%(s['remove_genome']['genome'])
    dehost_i_p = '\nPercentage of removed reads: %s'%(s['remove_genome']['remove_percentage'])

    aln_t = '## Alignment'
    aln_m_t = '### Mapping rate'
    aln_m_c = ''
    for ref_order in range(1, task.ref_num+1):
        aln_m_c += '\n'.join([
            '\n#### Reference #%d :%s'%(ref_order, s['ref_meta_dict']['seq_meta'][str(ref_order)]['fasta_header_escape']),
            '| Aligner | Overall mapped rate |',
            '| ------- | ------------------- |',
            '| Bowtie2 | %s |'%(s['aln']['mapped_rate']['bowtie2'][str(ref_order)]),
            '| BWA MEM | %s |'%(s['aln']['mapped_rate']['bwa'][str(ref_order)])
        ])
    
    aln_c_t = '### Coverage'
    aln_c_c = ''
    for ref_order in range(1, task.ref_num+1):
        aln_c_c += '\n'.join([
            '\n#### Reference #%d :%s'%(ref_order, s['ref_meta_dict']['seq_meta'][str(ref_order)]['fasta_header_escape']),
            '| Aligner | Start base | End base | Covered base | Mean depth |',
            '| ------- | ---------- | -------- | ------------ | ---------- |',
            '| Bowtie2 | %s | %s | %s %% | %s x |'%(s['cov']['bowtie2'][str(ref_order)]['startpos'], s['cov']['bowtie2'][str(ref_order)]['endpos'], s['cov']['bowtie2'][str(ref_order)]['coverage'], s['cov']['bowtie2'][str(ref_order)]['meandepth']),
            '| BWA MEM | %s | %s | %s %% | %s x |'%(s['cov']['bwa'][str(ref_order)]['startpos'], s['cov']['bwa'][str(ref_order)]['endpos'], s['cov']['bwa'][str(ref_order)]['coverage'], s['cov']['bwa'][str(ref_order)]['meandepth'])])
    vc_t = '## Variant Calling'
    vc_c = ''
    for ref_order in range(1, task.ref_num+1):
        vc_r = '### Reference #%d :%s'%(ref_order, s['ref_meta_dict']['seq_meta'][str(ref_order)]['fasta_header_escape'])
        vc_l_c = '#### LoFreq'
        if len(s['vc']['lofreq'][str(ref_order)]) > 0:
            vc_l_t = '| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |\n| -------- | --- | --- | ------- | --- |\n'
            for pos in s['vc']['lofreq'][str(ref_order)]:
                ref = s['vc']['lofreq'][str(ref_order)][pos]['REF']
                for alt, aln_result in s['vc']['lofreq'][str(ref_order)][pos]['SNV'].items():
                    bt2_ft = aln_result.get('bowtie2', {}).get('FILTER', '-')
                    bt2_af = aln_result.get('bowtie2', {}).get('FREQ', '-')
                    bt2_dp = aln_result.get('bowtie2', {}).get('DP', '-')
                    bt2_gq = aln_result.get('bowtie2', {}).get('QUAL', '-')
                    bwa_ft = aln_result.get('bwa', {}).get('FILTER', '-')
                    bwa_af = aln_result.get('bwa', {}).get('FREQ', '-')
                    bwa_dp = aln_result.get('bwa', {}).get('DP', '-')
                    bwa_gq = aln_result.get('bwa', {}).get('QUAL', '-')
                    vc_l_t += '| %s | %s | %s | %s / %s / %s / %s | %s / %s / %s / %s |\n'%(pos, ref, alt, bt2_ft, bt2_af, bt2_dp, bt2_gq, bwa_ft, bwa_af, bwa_dp, bwa_gq)
        else:
            vc_l_t = 'Lofreq did not report any SNV or indel.\n'
        vc_v_c = '#### Varscan2'
        if len(s['vc']['varscan'][str(ref_order)]) > 0:
            vc_v_t = '| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |\n| -------- | --- | --- | ------- | --- |\n'
            for pos in s['vc']['varscan'][str(ref_order)]:
                ref = s['vc']['varscan'][str(ref_order)][pos]['REF']
                for alt, aln_result in s['vc']['varscan'][str(ref_order)][pos]['SNV'].items():
                    bt2_ft = aln_result.get('bowtie2', {}).get('FILTER', '-')
                    bt2_af = aln_result.get('bowtie2', {}).get('FREQ', '-')
                    bt2_dp = aln_result.get('bowtie2', {}).get('DP', '-')
                    bt2_gq = aln_result.get('bowtie2', {}).get('QUAL', '-')
                    bwa_ft = aln_result.get('bwa', {}).get('FILTER', '-')
                    bwa_af = aln_result.get('bwa', {}).get('FREQ', '-')
                    bwa_dp = aln_result.get('bwa', {}).get('DP', '-')
                    bwa_gq = aln_result.get('bwa', {}).get('QUAL', '-')
                    vc_v_t += '| %s | %s | %s | %s / %s / %s / %s | %s / %s / %s / %s |\n'%(pos, ref, alt, bt2_ft, bt2_af, bt2_dp, bt2_gq, bwa_ft, bwa_af, bwa_dp, bwa_gq)
        else:
            vc_v_t = 'Varscan did not report any SNV or indel.\n'
        vc_c += '\n'.join([vc_r, vc_l_c, vc_l_t, vc_v_c, vc_v_t])
    genome_t = '## Draft Genome'
    genome_c = ''
    for ref_order in range(1, task.ref_num+1):
        genome_ref = '\n### Reference #%d :%s'%(ref_order, s['ref_meta_dict']['seq_meta'][str(ref_order)]['fasta_header_escape'])
        genome_pth = 'FASTA was saved to : %s' % s['draft_meta'][str(ref_order)]['file_path']
        genome_snv = ('Apllied SNV : %s' % s['draft_meta'][str(ref_order)]['snv_list']).replace('[]', 'None').replace("'", '').replace('[', '').replace(']', '')
        genome_cf = ('Conflict calling : %s' % s['draft_meta'][str(ref_order)]['conflicts']).replace('[]', 'None').replace("'", '').replace('[', '').replace(']', '')
        genome_mm = ('Mismatch calling : %s' % s['draft_meta'][str(ref_order)]['error']).replace('[]', 'None').replace("'", '').replace('[', '').replace(']', '')
        genome_c += '\n\n'.join([genome_ref, genome_pth, genome_snv, genome_cf, genome_mm])

    unmapped_t = '## Unmapped reads analysis'
    unmapped_c = ''
    if s['unmapped_analysis'] != {}:
        unmapped_as = 'Assembly mode: %s'%s['unmapped_analysis']['spades_mode']
        unmapped_db = 'BLAST database: %s'%s['unmapped_analysis']['BLASTdb_name']
        unmapped_hits = ''
        if s['unmapped_analysis']['highly_matched_result'] != []:
            unmapped_hits = 'BLAST Hits:\n\n| Hits | Acc. | Orgnism | Ident(%) | Query len.(bp) | Align len.(bp) | E-value |\n| ---- | ---- | ------- | -------- | ---------- | ---------- | ------- |\n'
            hit_order = 1
            for hit in s['unmapped_analysis']['highly_matched_result']:
                unmapped_hits += '| %d | %s | %s | %s | %s | %s | %s |\n'%(hit_order, hit['clean_sacc'], hit['clean_stitle_org'], hit['pident'], hit['qlen'], hit['length'], hit['evalue'])
                hit_order += 1
            unmapped_c += '\n\n'.join([unmapped_as, unmapped_db, unmapped_hits])
        else:
            unmapped_hits = 'N/A'
            unmapped_c += '\n\n'.join([unmapped_as, unmapped_db, unmapped_hits])
    else:
        unmapped_c = "N/A"

    cmd_t = '## Commands'
    cmd_c = '| Duration(s) | Executed command |\n| ---- | ------- |\n'
    for order, entry in s['log_dict'].items():
        if entry['string'].startswith('CMD:'):
            cmd_c += '| %d | %s |\n'%(entry['duration'], entry['string'][5:])
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
        ref_c_m,
        ref_c_t,
        fastp_t,
        fastp_f_t,
        fastp_f_c,
        fastp_d_t,
        fastp_d_c,
        dehost_t,
        dehost_i_g,
        dehost_i_p,
        aln_t,
        aln_m_t,
        aln_m_c,
        aln_c_t,
        aln_c_c,
        vc_t,
        vc_c,
        genome_t,
        genome_c,
        unmapped_t,
        unmapped_c,
        cmd_t,
        cmd_c,
        ver_t,
        ver_c
    ])

    print(md_str)
    utils.build_text_file(task.path.joinpath(task.id, task.id + '_report.md'), md_str)

def run(task):
    build_md_report(task)



if __name__ == "__main__":
    from new_task import Task
    task = Task()
    task.path = Path.cwd()
    task.id = 'test_run_202102090359'
    task.name = 'test_run'
    run(task)
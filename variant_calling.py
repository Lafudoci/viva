import logging
import subprocess
from decimal import Decimal
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def variant_calling_lofreq(task_name, base_path, aligner_list):
    logger.info('Starting variant calling by LoFreq.')
    thread_cmd = ['call-parallel', '--pp-threads', '6']
    other_cmd = ['-q', '20', '-Q', '20', '-m', '20']

    for aligner in aligner_list:
        logger.info('Running VC for %s output.' % aligner)
        aln_data_cwd = base_path.joinpath(task_name, 'alignment', aligner)
        aln_input_cmd = [str(aln_data_cwd.joinpath(task_name + '.sorted.bam'))]
        ref_path = aln_data_cwd.joinpath(task_name+'_ref.fasta')
        # index ref
        faidx_cmd = ['lofreq', 'faidx', str(ref_path)]
        logger.info('CMD: '+' '.join(faidx_cmd))
        utils.write_log_file(
            base_path.joinpath(task_name),
            'CMD: '+' '.join(faidx_cmd)
        )
        faidx_run = subprocess.run(
            faidx_cmd, cwd=aln_data_cwd, capture_output=True)
        print(faidx_run.stdout.decode(encoding='utf-8'))
        print(faidx_run.stderr.decode(encoding='utf-8'))
        # vc
        ref_cmd = ['-f', str(ref_path)]
        output_cmd = ['-o', '%s_%s_lofreq.vcf' % (task_name, aligner)]
        vc_cmd = ['lofreq'] + thread_cmd + \
            ref_cmd + output_cmd + \
            other_cmd + aln_input_cmd
        logger.info('CMD: '+' '.join(vc_cmd))
        utils.write_log_file(
            base_path.joinpath(task_name),
            'CMD: '+' '.join(vc_cmd)
        )
        vc_run = subprocess.run(vc_cmd, cwd=aln_data_cwd, capture_output=True)
        print(vc_run.stdout.decode(encoding='utf-8'))
        print(vc_run.stderr.decode(encoding='utf-8'))


def variant_calling_varscan2(task_name, base_path, aligner_list):
    logger.info('Starting variant calling by VarScan2.')
    mpileup_cmd = ['samtools', 'mpileup', '-B']
    mpileup2cns_cmd = [
        'java', '-jar',
        '/home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar',
        'mpileup2cns'
    ]
    output_cmd = ['--output-vcf', '1']
    other_cmd = ['--min-avg-qual', '20', '--P-value', '0.01']

    for aligner in aligner_list:
        logger.info('Running VC for %s output.' % aligner)
        aln_data_cwd = base_path.joinpath(task_name, 'alignment', aligner)
        aln_input_cmd = [str(aln_data_cwd.joinpath(task_name + '.sorted.bam'))]
        ref_path = aln_data_cwd.joinpath(task_name+'_ref.fasta')
        ref_cmd = ['-f', str(ref_path)]
        output_path = str(
            aln_data_cwd.joinpath(
                '%s_%s_varscan.vcf' % (task_name, aligner)
            )
        )
        # Run samtools mpileup and pipe to varscan2
        samtools_cmd = mpileup_cmd + ref_cmd + aln_input_cmd
        logger.info('CMD: '+' '.join(samtools_cmd))
        utils.write_log_file(
            base_path.joinpath(task_name),
            'CMD: '+' '.join(samtools_cmd)
        )
        samtools_run = subprocess.run(
            samtools_cmd,
            cwd=aln_data_cwd,
            capture_output=True
        )
        varscan2_cmd = mpileup2cns_cmd + other_cmd + output_cmd
        logger.info('CMD: '+' '.join(varscan2_cmd))
        utils.write_log_file(
            base_path.joinpath(task_name),
            'CMD: '+' '.join(varscan2_cmd)
        )
        vc_run = subprocess.run(
            varscan2_cmd,
            cwd=aln_data_cwd,
            input=samtools_run.stdout,
            capture_output=True
        )
        utils.build_text_file(
            output_path, vc_run.stdout.decode(encoding='utf-8'))
        print(vc_run.stderr.decode(encoding='utf-8'))


def build_vc_summary_json(task_name, base_path, aligner_list):
    vc_summary_path = base_path.joinpath(
        task_name, task_name+'_vc_summary.json'
    )
    varscan_vc_dict = {}
    lofreq_vc_dict = {}
    for aligner in aligner_list:
        # parse varscan vcf file
        varscan_vcf_file_path = base_path.joinpath(
            task_name, 'alignment', aligner, '%s_%s_varscan.vcf' % (
                task_name, aligner)
        )
        original_varscan_vc_dict = utils.load_vcf_file(varscan_vcf_file_path)
        for pos_result in original_varscan_vc_dict['vc']:
            if pos_result[4] != '.':
                pos = pos_result[1]
                ref = pos_result[3]
                alt = pos_result[4]
                qual = pos_result[9].split(':')[1]
                vc_filter = pos_result[6]
                freq = pos_result[9].split(':')[6]
                sdp = pos_result[9].split(':')[2]
                if varscan_vc_dict.get(pos) == None:
                    varscan_vc_dict[pos] = {'REF': ref, 'SNV': {}}
                if varscan_vc_dict[pos]['SNV'].get(alt) == None:
                    varscan_vc_dict[pos]['SNV'][alt] = {}
                if varscan_vc_dict[pos]['SNV'][alt].get(aligner) == None:
                    varscan_vc_dict[pos]['SNV'][alt][aligner] = {}
                varscan_vc_dict[pos]['SNV'][alt][aligner] = {
                    'FILTER': vc_filter, 'FREQ': freq, 'QUAL': qual, 'DP': sdp
                }
        # parse lofreq vcf file
        lofreq_vcf_file_path = base_path.joinpath(
            task_name, 'alignment', aligner, '%s_%s_lofreq.vcf' % (
                task_name, aligner)
        )
        original_lofreq_vc_dict = utils.load_vcf_file(lofreq_vcf_file_path)
        for pos_result in original_lofreq_vc_dict['vc']:
            af = pos_result[7].split(';')[1].split('=')[1]
            if Decimal(af) >= Decimal('0.1'):
                pos = pos_result[1]
                ref = pos_result[3]
                alt = pos_result[4]
                qual = pos_result[5]
                vc_filter = pos_result[6]
                af_formated = '%.02f%%' % (Decimal(af)*Decimal('100'))
                dp = pos_result[7].split(';')[0].split('=')[1]
                if lofreq_vc_dict.get(pos) == None:
                    lofreq_vc_dict[pos] = {'REF': ref, 'SNV': {}}
                if lofreq_vc_dict[pos]['SNV'].get(alt) == None:
                    lofreq_vc_dict[pos]['SNV'][alt] = {}
                if lofreq_vc_dict[pos]['SNV'][alt].get(aligner) == None:
                    lofreq_vc_dict[pos]['SNV'][alt][aligner] = {}
                lofreq_vc_dict[pos]['SNV'][alt][aligner] = {
                    'FILTER': vc_filter, 'FREQ': af_formated, 'QUAL': qual, 'DP': dp
                }

    all_vc_dict = {'lofreq': lofreq_vc_dict, 'varscan': varscan_vc_dict}
    utils.build_json_file(vc_summary_path, all_vc_dict)


def run(task_name, base_path):
    aligner_list = ['bowtie2', 'bwa']
    variant_calling_lofreq(task_name, base_path, aligner_list)
    variant_calling_varscan2(task_name, base_path, aligner_list)
    build_vc_summary_json(task_name, base_path, aligner_list)


if __name__ == '__main__':
    build_vc_summary_json('test', Path.cwd(), ['bowtie2', 'bwa'])

import logging
import subprocess
from decimal import Decimal
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def variant_calling_lofreq(task):
    logger.info('Starting variant calling by LoFreq.')
    thread_cmd = ['call-parallel', '--pp-threads', str(task.threads)]
    other_cmd = ['--call-indels', '-N', '-B', '-q', '20', '-Q', '20', '-m', '20']

    for aligner in task.alns:
        logger.info('Running VC for %s output.' % aligner)
        aln_data_cwd = task.path.joinpath(task.id, 'alignment', aligner)
        aln_input_name = task.id + '.sorted.bam'
        aln_indelqual_name = task.id + '.indelqual.sorted.bam'
        ref_name = task.id + '_ref.fasta'
        # index ref
        faidx_cmd = ['lofreq', 'faidx', ref_name]
        logger.info('CMD: '+' '.join(faidx_cmd))
        utils.write_log_file(
            task.path.joinpath(task.id),
            'CMD: '+' '.join(faidx_cmd)
        )
        faidx_run = subprocess.run(
            faidx_cmd, cwd=aln_data_cwd, capture_output=True)
        print(faidx_run.stdout.decode(encoding='utf-8'))
        print(faidx_run.stderr.decode(encoding='utf-8'))
        # indelqual
        indelqual_cmd = ['lofreq', 'indelqual', '--dindel', '--ref', ref_name, '--out', aln_indelqual_name, aln_input_name]
        indelqual_run = subprocess.run(
            indelqual_cmd, cwd=aln_data_cwd, capture_output=True)
        print(indelqual_run.stdout.decode(encoding='utf-8'))
        print(indelqual_run.stderr.decode(encoding='utf-8'))
        # index indelqual-ed BAM
        indelqual_index_cmd = ['samtools', 'index', aln_indelqual_name]
        indelqual_index_run = subprocess.run(
            indelqual_index_cmd, cwd=aln_data_cwd, capture_output=True)
        print(indelqual_index_run.stdout.decode(encoding='utf-8'))
        print(indelqual_index_run.stderr.decode(encoding='utf-8'))
        # vc
        ref_cmd = ['-f', ref_name]
        output_cmd = ['-o', '%s_%s_lofreq.vcf' % (task.id, aligner)]
        vc_cmd = ['lofreq'] + thread_cmd + \
            ref_cmd + output_cmd + \
            other_cmd + [aln_indelqual_name]
        logger.info('CMD: '+' '.join(vc_cmd))
        utils.write_log_file(
            task.path.joinpath(task.id),
            'CMD: '+' '.join(vc_cmd)
        )
        vc_run = subprocess.run(vc_cmd, cwd=aln_data_cwd, capture_output=True)
        print(vc_run.stdout.decode(encoding='utf-8'))
        print(vc_run.stderr.decode(encoding='utf-8'))


def variant_calling_varscan2(task):
    logger.info('Starting variant calling by VarScan2.')
    mpileup_cmd = ['samtools', 'mpileup', '-B']
    mpileup2cns_cmd = [
        'java', '-jar',
        '/home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar',
        'mpileup2cns'
    ]
    output_cmd = ['--output-vcf', '1']
    other_cmd = ['--min-avg-qual', '20', '--P-value', '0.01']

    for aligner in task.alns:
        logger.info('Running VC for %s output.' % aligner)
        aln_data_cwd = task.path.joinpath(task.id, 'alignment', aligner)
        aln_input_cmd = [str(aln_data_cwd.joinpath(task.id + '.sorted.bam'))]
        ref_path = aln_data_cwd.joinpath(task.id+'_ref.fasta')
        ref_cmd = ['-f', str(ref_path)]
        output_path = str(
            aln_data_cwd.joinpath(
                '%s_%s_varscan.vcf' % (task.id, aligner)
            )
        )
        # Run samtools mpileup and pipe to varscan2
        samtools_cmd = mpileup_cmd + ref_cmd + aln_input_cmd
        logger.info('CMD: '+' '.join(samtools_cmd))
        utils.write_log_file(
            task.path.joinpath(task.id),
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
            task.path.joinpath(task.id),
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


def build_vc_summary_json(task):
    vc_summary_path = task.path.joinpath(
        task.id, task.id+'_vc_summary.json'
    )
    varscan_vc_dict = {}
    lofreq_vc_dict = {}
    for aligner in task.alns:
        # parse varscan vcf file
        varscan_vcf_file_path = task.path.joinpath(
            task.id, 'alignment', aligner, '%s_%s_varscan.vcf' % (
                task.id, aligner)
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
        lofreq_vcf_file_path = task.path.joinpath(
            task.id, 'alignment', aligner, '%s_%s_lofreq.vcf' % (
                task.id, aligner)
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


def build_draft_genome_seq(task):
    draft_genome_summary = {'conflicts':[], 'snv_list':[], 'error':[], 'file_path':''}
    vc_summary_path = task.path.joinpath(
        task.id, task.id+'_vc_summary.json'
    )
    vc_dict = utils.load_json_file(vc_summary_path)
    dominant_vc = {}
    for caller, vc_table in vc_dict.items():
        for pos, snvs in vc_table.items():
            ref = snvs['REF']
            for snv, alns in snvs['SNV'].items():
                for aligner in alns:
                    if Decimal(alns[aligner]['FREQ'][:-1])/100 > Decimal(task.vc_threshold):
                        if dominant_vc.get(pos) == None:
                            dominant_vc[pos] = {'REF':ref, 'ALT':{}}
                        if dominant_vc[pos]['ALT'].get(snv) == None:
                            dominant_vc[pos]['ALT'].update({snv: {'SCORE':0}})
                        dominant_vc[pos]['ALT'][snv]['SCORE'] += 1
    
    fasta_base_list = []
    imported_ref = task.path.joinpath(task.id, 'reference', task.id + '_ref.fasta')
    ref_fasta_dict = utils.load_fasta_file(imported_ref)
    for base in ref_fasta_dict[task.id+'_ref']:
        fasta_base_list.append(base)

    for pos, vc in dominant_vc.items():
        if len(vc['ALT']) > 2:
            # skip conflict results
            draft_genome_summary['conflicts'].append(pos)
        else:
            # apply snv onto reference sequence
            i = 0
            for base in vc['REF'][1:]:
                if fasta_base_list[int(pos)+i] != base:
                    draft_genome_summary['error'].append(pos)
                    break
                fasta_base_list[int(pos)+i] = ''
                i += 1
            fasta_base_list[int(pos)-1] = list(vc)[0]
            # record apllied snv
            ref_mer = vc['REF']
            alt_mer = list(vc['ALT'].keys())[0]
            print(ref_mer, alt_mer)
            if len(ref_mer) == len(alt_mer):
                draft_genome_summary['snv_list'].append('%s%s%s'%(ref_mer, pos, alt_mer))
    
    draft_fasta_dict = { '%s_draft'%task.id :''.join(fasta_base_list) }

    draft_cwd = task.path.joinpath(task.id, 'draft_genome')
    draft_fasta_path = draft_cwd.joinpath('%s_draft.fasta'%task.id)
    draft_genome_summary['file_path'] = str(draft_fasta_path)
    Path.mkdir(draft_cwd, parents=True, exist_ok=True)
    utils.build_fasta_file(draft_fasta_path, draft_fasta_dict)
    utils.build_json_file(draft_cwd.joinpath('%s_draft_summary.json'%task.id), draft_genome_summary)


def run(task):
    variant_calling_lofreq(task)
    variant_calling_varscan2(task)
    build_vc_summary_json(task)
    build_draft_genome_seq(task)

if __name__ == '__main__':
    import summary_generator
    from new_task import Task
    task = Task()
    task.path = Path.cwd()
    task.id = 'test_run_202102090359'
    task.vc_threshold = '0.7'
    build_draft_genome_seq(task)
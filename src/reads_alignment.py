import logging
import os
import shutil
import subprocess
from decimal import Decimal
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def ref_index(task, aligner):
    for ref_order in range(1, task.ref_num+1):
        logger.info('Building Bowtie2 ref index for ref #%d.'%ref_order)
        aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
        ref_fasta_path = task.path.joinpath(
            task.id, 'reference', '%s_ref_%d.fasta'%(task.id, ref_order))
        Path.mkdir(aligner_cwd, parents=True, exist_ok=True)
        shutil.copy2(ref_fasta_path, aligner_cwd)
        if aligner == 'bowtie2':
            index_cmd = ['bowtie2-build', '--threads', task.threads, '%s_ref_%d.fasta'%(task.id, ref_order), '%s_ref_%d'%(task.id, ref_order)]
        elif aligner == 'bwa':
            index_cmd = ['bwa', 'index', '-p', '%s_ref_%d'%(task.id, ref_order), '%s_ref_%d.fasta'%(task.id, ref_order)]
        logger.info('CMD: '+' '.join(index_cmd))
        utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(index_cmd))
        ref_index_run = subprocess.run(index_cmd, cwd=aligner_cwd, capture_output=True)
        print(ref_index_run.stdout.decode(encoding='utf-8'))
        print(ref_index_run.stderr.decode(encoding='utf-8'))


def align_bowtie2(task):
    for ref_order in range(1, task.ref_num+1):
        logger.info('Running Bowtie2 alignment for ref #%d.'%ref_order)
        aligner_cwd = task.path.joinpath(task.id, 'alignment', 'bowtie2')
        ref_index_path = str(aligner_cwd.joinpath('%s_ref_%d'%(task.id, ref_order)))
        # impurities_prefilter is latter process, so need to be check first, then remove_host
        if task.impurities_prefilter_num > 0:
            filterd_R1 = str(task.path.joinpath(
                task.id, 'reads', '%s_%d_impurity_removed_R1.fastq.gz' % (task.id, task.impurities_prefilter_num)))
            filterd_R2 = str(task.path.joinpath(
                task.id, 'reads', '%s_%d_impurity_removed_R2.fastq.gz' % (task.id, task.impurities_prefilter_num)))
        elif task.remove_host != None:
            filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R1.fastq.gz'))
            filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R2.fastq.gz'))
        else:
            filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz'))
            filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz'))
        reads_cmd = ['-1', filterd_R1, '-2', filterd_R2]
        thread_cmd = ['-p', str(task.threads)]
        output_cmd = ['-S', '%s_ref_%d.sam'%(task.id, ref_order)]
        other_cmd = ['--very-sensitive-local']
        aln_cmd = ['bowtie2', '-x', ref_index_path] + reads_cmd + output_cmd + thread_cmd + other_cmd
        logger.info('CMD: '+' '.join(aln_cmd))
        utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(aln_cmd))
        bt2_run = subprocess.run(aln_cmd, cwd=aligner_cwd, capture_output=True)
        print(bt2_run.stdout.decode(encoding='utf-8'))
        print(bt2_run.stderr.decode(encoding='utf-8'))
        bam_sort_n_index(task, 'bowtie2', ref_order)


def align_bwa(task):
    for ref_order in range(1, task.ref_num+1):
        logger.info('Running BWA alignment for ref #%d.'%ref_order)
        aligner_cwd = task.path.joinpath(task.id, 'alignment', 'bwa')
        ref_index_path = str(aligner_cwd.joinpath('%s_ref_%d'%(task.id, ref_order)))
        # impurities_prefilter is latter process, so need to be check first, then remove_host
        if task.impurities_prefilter_num > 0:
            filterd_R1 = str(task.path.joinpath(
                task.id, 'reads', '%s_%d_impurity_removed_R1.fastq.gz' % (task.id, task.impurities_prefilter_num)))
            filterd_R2 = str(task.path.joinpath(
                task.id, 'reads', '%s_%d_impurity_removed_R2.fastq.gz' % (task.id, task.impurities_prefilter_num)))
        elif task.remove_host != None:
            filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R1.fastq.gz'))
            filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R2.fastq.gz'))
        else:
            filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz'))
            filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz'))
        reads_cmd = [filterd_R1, filterd_R2]
        option_cmd = ['-K 100000000 -Y']
        thread_cmd = ['-t', str(task.threads)]
        output_cmd = ['-o', '%s_ref_%d.sam'%(task.id, ref_order)]
        aln_cmd = ['bwa', 'mem'] + option_cmd + thread_cmd + [ref_index_path] + reads_cmd + output_cmd
        logger.info('CMD: '+' '.join(aln_cmd))
        utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(aln_cmd))
        bt2_run = subprocess.run(aln_cmd, cwd=aligner_cwd, capture_output=True)
        print(bt2_run.stdout.decode(encoding='utf-8'))
        print(bt2_run.stderr.decode(encoding='utf-8'))
        bam_sort_n_index(task, 'bwa', ref_order)


def bam_sort_n_index(task, aligner, ref_order):
    logger.info('Sorting & indexing BAM file for aln #%d.'%ref_order)
    aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
    # sorting
    sorting_cmd = ['samtools', 'sort', '-@', task.threads, '%s_ref_%d.sam'%(task.id, ref_order), '-o', '%s_ref_%d.sorted.bam'%(task.id, ref_order)]
    logger.info('CMD: '+' '.join(sorting_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(sorting_cmd))
    sorting_run = subprocess.run(sorting_cmd, cwd=aligner_cwd, capture_output=True)
    print(sorting_run.stdout.decode(encoding='utf-8'))
    print(sorting_run.stderr.decode(encoding='utf-8'))
    # indexing
    indexing_cmd = ['samtools', 'index', '-@', task.threads, '%s_ref_%d.sorted.bam'%(task.id, ref_order)]
    logger.info('CMD: '+' '.join(indexing_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(indexing_cmd))
    index_run = subprocess.run(indexing_cmd, cwd=aligner_cwd, capture_output=True)
    print(index_run.stdout.decode(encoding='utf-8'))
    print(index_run.stderr.decode(encoding='utf-8'))
    # remove sam file to release disk space
    os.remove(aligner_cwd.joinpath('%s_ref_%d.sam'%(task.id, ref_order)))


def align_flagstat(task, aligners):
    stats_dict = {'mapped_rate':{}, 'mapped_reads':{}}
    for aligner in aligners:
        logger.info('Analysis BAM file from %s' % aligner)
        aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
        stats_dict['mapped_rate'][aligner] = {}
        stats_dict['mapped_reads'][aligner] = {}
        for ref_order in range(1, task.ref_num+1):
            flagstat_cmd = ['samtools', 'flagstat', '-@', task.threads, '%s_ref_%d.sorted.bam'%(task.id, ref_order)]
            logger.info('CMD: '+' '.join(flagstat_cmd))
            utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(flagstat_cmd))
            flagstat_run = subprocess.run(flagstat_cmd, cwd=aligner_cwd, capture_output=True)
            stats_text = flagstat_run.stdout.decode(encoding='utf-8')
            flagstat_file_path = task.path.joinpath(aligner_cwd, 'flagstat_ref_%d.txt'%ref_order)
            utils.build_text_file(flagstat_file_path, stats_text)
            primary_mapped_reads = utils.primary_mapped_from_flagstat(flagstat_file_path)
            mapped_rate = Decimal(primary_mapped_reads)/Decimal(task.total_reads_after_fastp)
            stats_dict['mapped_rate'][aligner][ref_order]= "%f%%" % (mapped_rate*Decimal('100'))
            stats_dict['mapped_reads'][aligner][ref_order]= primary_mapped_reads
    utils.build_json_file(task.path.joinpath(task.id, 'alignment', 'flagstat.json'), stats_dict)


def align_coverage_stat(task, aligners):
    cov_dict = {}
    for aligner in aligners:
        cov_dict[aligner] = {}
        logger.info('Analysis coverage stats from %s BAM files.' % aligner)
        aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
        for ref_order in range(1, task.ref_num+1):
            cov_dict[aligner][ref_order] = {}
            flagstat_cmd = ['samtools', 'coverage', '%s_ref_%d.sorted.bam'%(task.id, ref_order)]
            logger.info('CMD: '+' '.join(flagstat_cmd))
            utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(flagstat_cmd))
            flagstat_run = subprocess.run(flagstat_cmd, cwd=aligner_cwd, capture_output=True)
            stats_text = flagstat_run.stdout.decode(encoding='utf-8')
            titles = stats_text.split('\n')[0].split('\t')
            stats = stats_text.split('\n')[1].split('\t')
            for i in range(len(titles)):
                cov_dict[aligner][ref_order][titles[i]] = stats[i]
    utils.build_json_file(task.path.joinpath(task.id, 'alignment', 'coverage_stat.json'), cov_dict)


def extract_unmapped_reads(task, aligners):
    for aligner in aligners:
        logger.info('Extract unmapped reads from %s BAM files.' % aligner)
        aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
        for ref_order in range(1, task.ref_num+1):
            samtools_option_cmd = ['samtools', 'fastq', '-f 13']
            samtools_fastq_cmd = ['-1', '%s_ref_%d_unmapped_R1.fastq.gz'%(task.id, ref_order), '-2', '%s_ref_%d_unmapped_R2.fastq.gz'%(task.id, ref_order)]
            samtools_run_cmd = samtools_option_cmd + samtools_fastq_cmd + ['%s_ref_%d.sorted.bam'%(task.id, ref_order)]
            subprocess.run(samtools_run_cmd, cwd=aligner_cwd, check=True)


def align_disp(task, aligner):
    if aligner == 'bowtie2':
        align_bowtie2(task)
    elif aligner == 'bwa':
        align_bwa(task)

def run(task):
    aligners = task.alns
    for aligner in aligners:
        ref_index(task, aligner)
        align_disp(task, aligner)
    align_flagstat(task, aligners)
    align_coverage_stat(task, aligners)
    extract_unmapped_reads(task, aligners)

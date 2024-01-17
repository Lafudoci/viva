import hashlib
import logging
import os
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import summary_generator
import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def import_reads(task):
    external_reads_R1 = Path(task.ex_r1)
    external_reads_R2 = Path(task.ex_r2)
    original_reads_path = task.path.joinpath(task.id, 'reads', 'original')
    Path.mkdir(original_reads_path, parents=True, exist_ok=True)
    if str(external_reads_R1).endswith('.fastq.gz') and str(external_reads_R2).endswith('.fastq.gz'):
        shutil.copy(external_reads_R1, original_reads_path.joinpath(
            task.id + '_R1.fastq.gz'))
        shutil.copy(external_reads_R2, original_reads_path.joinpath(
            task.id + '_R2.fastq.gz'))
        md5 = reads_hash_md5(task)
        reads_meta = {
            'file_name': {'r1': str(external_reads_R1), 'r2':str(external_reads_R2)},
            'md5': {'r1': md5[0], 'r2': md5[1]}
        }
        reads_meta_path = task.path.joinpath(task.id, 'reads', 'reads_meta.json')
        
        utils.build_json_file(reads_meta_path, reads_meta)
    else:
        logger.error('Unsupported reads file type.')


def run_fastp(task):
    logger.info('Running fastp to filter reads.')
    original_R1 = task.path.joinpath(
        task.id, 'reads', 'original', task.id + '_R1.fastq.gz')
    original_R2 = task.path.joinpath(
        task.id, 'reads', 'original', task.id + '_R2.fastq.gz')
    filterd_R1 = task.path.joinpath(
        task.id, 'reads', task.id + '_R1.fastq.gz')
    filterd_R2 = task.path.joinpath(
        task.id, 'reads', task.id + '_R2.fastq.gz')
    report_json = task.path.joinpath(task.id, 'reads', 'fastp.json')
    report_html = task.path.joinpath(task.id, 'reads', 'fastp.html')
    fastp_cmd = [
        "fastp",
        "-i", str(original_R1),
        "-I", str(original_R2),
        "-o", str(filterd_R1),
        "-O", str(filterd_R2)
    ]
    reports_cmd = [
        '-j', str(report_json),
        '-h', str(report_html)
    ]
    parameter_cmd = [
        '-f', str(task.global_trimming),
        '-t', str(task.global_trimming),
        '-F', str(task.global_trimming),
        '-T', str(task.global_trimming),
        '-w', str(task.threads)
    ]
    logger.info('CMD: '+' '.join(fastp_cmd + reports_cmd + parameter_cmd))
    utils.write_log_file(
        task.path.joinpath(task.id),
        'CMD: '+' '.join(fastp_cmd + reports_cmd + parameter_cmd)
    )
    fastp_run = subprocess.run(
        fastp_cmd + reports_cmd + parameter_cmd,
        capture_output=True
    )
    print(fastp_run.stderr.decode(encoding='utf-8'))
    # record total reads after fastp
    task.total_reads_after_fastp = summary_generator.fastp_parser(task)['after_total_reads']
    mininal_reads = 100
    if task.total_reads_after_fastp <= mininal_reads:
        logger.critical('Reads too less (%s) to run, minimal is %s'%(task.total_reads_after_fastp, mininal_reads))
        sys.exit(-1)



def reads_hash_md5(task):
    logger.info('Calculating reads MD5 hash.')
    original_reads_path = task.path.joinpath(task.id, 'reads', 'original')
    R1_path = task.path.joinpath(original_reads_path, task.id + '_R1.fastq.gz')
    R2_path = task.path.joinpath(original_reads_path, task.id + '_R2.fastq.gz')
    r1 = hashlib.md5()
    with open(R1_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            r1.update(chunk)
    R1_hash = r1.hexdigest()
    r2 = hashlib.md5()
    with open(R2_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            r2.update(chunk)
    R2_hash = r2.hexdigest()
    return (R1_hash, R2_hash)


def remove_host(task):
    logger.info('Removing host genome.')
    dehost_meta = {'genome': '', 'remove_percentage': ''}
    host_remove_cwd = task.path.joinpath(task.id, 'reads')
    Path.mkdir(host_remove_cwd, parents=True, exist_ok=True)

    if task.remove_host == 'dog':
        dehost_meta['genome'] = 'Dog (Dog10K_Boxer_Tasha, GCF_000002285.5)'
        genome_path = '/app/genomes/' + 'dog10k'
    elif task.remove_host == 'human':
        dehost_meta['genome'] = 'Human (GRCh38.p13, GCF_000001405.39)'
        genome_path = '/app/genomes/' + 'grch38'
    elif task.remove_host == 'vero':
        dehost_meta['genome'] = 'Vero (Vero_WHO_p1.0, GCF_015252025.1)'
        genome_path = '/app/genomes/' +  'vero'
    elif task.remove_host == 'chicken':
        dehost_meta['genome'] = 'Chicken (GRCg6a, GCF_000002315.6)'
        genome_path = '/app/genomes/' + 'grcg6a'
    elif task.remove_host == 'rhesus_monkey':
        dehost_meta['genome'] = 'Rhesus monkey (Mmul_10, GCF_003339765.1)'
        genome_path = '/app/genomes/' + 'mmul_10'
    else:
        dehost_meta['genome'] = 'Custom sequence file (%s)'%task.remove_host
        genome_path = '/app/genomes/' + task.remove_host

    mapped_reads_out = 'host_mapped.sam'
    align_cmd = [
        'bowtie2',
        '-p', str(task.threads),
        '-x', str(genome_path),
        '-1', str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz')),
        '-2', str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz')),
        '-S', str(mapped_reads_out),
        '--very-sensitive-local'
    ]
    logger.info('CMD: '+' '.join(align_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(align_cmd))
    cmd_run = subprocess.run(align_cmd, cwd=host_remove_cwd, capture_output=True)
    # print(cmd_run.stdout.decode(encoding='utf-8'))
    print(cmd_run.stderr.decode(encoding='utf-8'))

    # build meta
    logger.info('Analysis BAM file from host mapped reads')
    # sorting
    sorting_cmd = ['samtools', 'sort', '-@', task.threads, 'host_mapped.sam', '-o', 'host_mapped.sorted.bam']
    logger.info('CMD: '+' '.join(sorting_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(sorting_cmd))
    sorting_run = subprocess.run(sorting_cmd, cwd=host_remove_cwd, capture_output=True)
    print(sorting_run.stdout.decode(encoding='utf-8'))
    print(sorting_run.stderr.decode(encoding='utf-8'))
    # extract unmapped
    unmapped_fastq_r1 = task.id + '_host_removed_R1.fastq.gz'
    unmapped_fastq_r2 = task.id + '_host_removed_R2.fastq.gz'
    samtools_option_cmd = ['samtools', 'fastq', '-f 13']
    samtools_fastq_cmd = ['-1', unmapped_fastq_r1, '-2', unmapped_fastq_r2]
    samtools_run_cmd = samtools_option_cmd + samtools_fastq_cmd + ['host_mapped.sorted.bam']
    subprocess.run(samtools_run_cmd, cwd=host_remove_cwd, check=True)
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(samtools_run_cmd))
    # flagstat
    flagstat_cmd = ['samtools', 'flagstat', '-@', task.threads, 'host_mapped.sorted.bam']
    logger.info('CMD: '+' '.join(flagstat_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(flagstat_cmd))
    flagstat_run = subprocess.run(flagstat_cmd, cwd=host_remove_cwd, capture_output=True)
    stats_text = flagstat_run.stdout.decode(encoding='utf-8')
    flagstat_file_path = task.path.joinpath(host_remove_cwd, 'flagstat.txt')
    utils.build_text_file(flagstat_file_path, stats_text)
    total_reads = task.total_reads_after_fastp
    primary_mapped_reads = utils.primary_mapped_from_flagstat(flagstat_file_path)
    mapped_rate = Decimal(primary_mapped_reads)/Decimal(total_reads)
    dehost_meta['mapped_reads'] = primary_mapped_reads
    dehost_meta['remove_percentage'] = "%f%%" % (mapped_rate*Decimal('100'))
    utils.build_json_file(task.path.joinpath(host_remove_cwd, 'dehost_meta.json'), dehost_meta)
    # remove sam file to release disk space
    os.remove(task.path.joinpath(host_remove_cwd, mapped_reads_out))
    # remove host bam file to release disk space
    os.remove(task.path.joinpath(host_remove_cwd, 'host_mapped.sorted.bam'))

def run(task):
    logger.info('Importing reads.')
    import_reads(task)
    run_fastp(task)
    if task.remove_host != None:
        remove_host(task)

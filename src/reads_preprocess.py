import hashlib
import logging
import os
import shutil
import subprocess
from pathlib import Path

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

    if task.dehost == 'dog':
        dehost_meta['genome'] = 'Dog (Dog10K_Boxer_Tasha, GCF_000002285.5)'
        genome_path = Path('genomes', 'dog', 'dog10k', 'dog10k')
    elif task.dehost == 'human':
        dehost_meta['genome'] = 'Human (GRCh38.p13, GCF_000001405.39)'
        genome_path = '/app/genomes/' + 'grch38'
    elif task.dehost == 'vero':
        dehost_meta['genome'] = 'Vero (Vero_WHO_p1.0, GCF_015252025.1)'
        genome_path = Path('genomes', 'vero', 'vero')
    elif task.dehost == 'chicken':
        dehost_meta['genome'] = 'Chicken (GRCg6a, GCF_000002315.6)'
        genome_path = Path('genomes', 'grcg6a', 'grcg6a')
    elif task.dehost == 'rhesus_monkey':
        dehost_meta['genome'] = 'Rhesus monkey (Mmul_10, GCF_003339765.1)'
        genome_path = Path('genomes', 'mmul_10', 'mmul_10')

    unconc_reads_out = task.id + '_host_removed_R%.fastq.gz'
    mapped_reads_out = 'host_mapped.sam'
    align_cmd = [
        'bowtie2',
        '-p', str(task.threads),
        '-x', str(genome_path),
        '-1', str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz')),
        '-2', str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz')),
        '-S', str(mapped_reads_out),
        '--very-sensitive-local',
        '--un-conc-gz', '%s' % str(unconc_reads_out)
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
    # flagstat
    flagstat_cmd = ['samtools', 'flagstat', '-@', task.threads, 'host_mapped.sorted.bam']
    logger.info('CMD: '+' '.join(flagstat_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(flagstat_cmd))
    flagstat_run = subprocess.run(flagstat_cmd, cwd=host_remove_cwd, capture_output=True)
    stats_text = flagstat_run.stdout.decode(encoding='utf-8')
    stats_list = stats_text.split('\n')
    utils.build_text_file(task.path.joinpath(host_remove_cwd, 'flagstat.txt'), stats_text)
    mapped_rate = stats_list[4].split(' ')[4][1:]
    dehost_meta['remove_percentage'] = mapped_rate
    utils.build_json_file(task.path.joinpath(host_remove_cwd, 'dehost_meta.json'), dehost_meta)
    # remove sam file to release disk space
    os.remove(task.path.joinpath(host_remove_cwd, mapped_reads_out))
    # remove host bam file to release disk space
    os.remove(task.path.joinpath(host_remove_cwd, 'host_mapped.sorted.bam'))

def run(task):
    logger.info('Importing reads.')
    import_reads(task)
    run_fastp(task)
    if task.dehost != None:
        remove_host(task)

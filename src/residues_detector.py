import logging
import os
import subprocess
from decimal import Decimal
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
        genome_path = '/app/genomes/' + 'vero'
    elif task.remove_host == 'chicken':
        dehost_meta['genome'] = 'Chicken (GRCg6a, GCF_000002315.6)'
        genome_path = '/app/genomes/' + 'grcg6a'
    elif task.remove_host == 'rhesus_monkey':
        dehost_meta['genome'] = 'Rhesus monkey (Mmul_10, GCF_003339765.1)'
        genome_path = '/app/genomes/' + 'mmul_10'
    else:
        dehost_meta['genome'] = 'Custom sequence file (%s)' % task.remove_host
        genome_path = '/app/genomes/' + task.remove_host

    unconc_reads_out = task.id + '_host_removed_R%.fastq.gz'
    mapped_reads_out = 'host_mapped.sam'
    align_cmd = [
        'bowtie2',
        '-p', str(task.threads),
        '-x', str(genome_path),
        '-1', str(task.path.joinpath(task.id,
                                     'reads', task.id + '_R1.fastq.gz')),
        '-2', str(task.path.joinpath(task.id,
                                     'reads', task.id + '_R2.fastq.gz')),
        '-S', str(mapped_reads_out),
        '--very-sensitive-local',
        '--un-conc-gz', '%s' % str(unconc_reads_out)
    ]
    logger.info('CMD: '+' '.join(align_cmd))
    utils.write_log_file(task.path.joinpath(task.id),
                         'CMD: '+' '.join(align_cmd))
    cmd_run = subprocess.run(
        align_cmd, cwd=host_remove_cwd, capture_output=True)
    # print(cmd_run.stdout.decode(encoding='utf-8'))
    print(cmd_run.stderr.decode(encoding='utf-8'))

    # build meta
    logger.info('Analysis BAM file from host mapped reads')
    # sorting
    sorting_cmd = ['samtools', 'sort', '-@', task.threads,
                   'host_mapped.sam', '-o', 'host_mapped.sorted.bam']
    logger.info('CMD: '+' '.join(sorting_cmd))
    utils.write_log_file(task.path.joinpath(task.id),
                         'CMD: '+' '.join(sorting_cmd))
    sorting_run = subprocess.run(
        sorting_cmd, cwd=host_remove_cwd, capture_output=True)
    print(sorting_run.stdout.decode(encoding='utf-8'))
    print(sorting_run.stderr.decode(encoding='utf-8'))
    # flagstat
    flagstat_cmd = ['samtools', 'flagstat', '-@',
                    task.threads, 'host_mapped.sorted.bam']
    logger.info('CMD: '+' '.join(flagstat_cmd))
    utils.write_log_file(task.path.joinpath(task.id),
                         'CMD: '+' '.join(flagstat_cmd))
    flagstat_run = subprocess.run(
        flagstat_cmd, cwd=host_remove_cwd, capture_output=True)
    stats_text = flagstat_run.stdout.decode(encoding='utf-8')
    stats_list = stats_text.split('\n')
    utils.build_text_file(task.path.joinpath(
        host_remove_cwd, 'flagstat.txt'), stats_text)
    total_reads = stats_list[0].split(' ')[0]
    mapped_reads = stats_list[4].split(' ')[0]
    mapped_rate = Decimal(mapped_reads)/Decimal(total_reads)
    dehost_meta['remove_percentage'] = "%f%%" % (mapped_rate*Decimal('100'))
    utils.build_json_file(task.path.joinpath(
        host_remove_cwd, 'dehost_meta.json'), dehost_meta)
    # remove sam file to release disk space
    os.remove(task.path.joinpath(host_remove_cwd, mapped_reads_out))
    # remove host bam file to release disk space
    os.remove(task.path.joinpath(host_remove_cwd, 'host_mapped.sorted.bam'))


def run(task):
    logger.info('Detecting residues.')
    if task.remove_host != None:
        remove_host(task)
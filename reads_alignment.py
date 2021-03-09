import logging
import subprocess
import utils
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def ref_index(task, aligner):
    logger.info('Building Bowtie2 ref index.')
    aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
    ref_fasta_path = task.path.joinpath(
        task.id, 'reference', task.id+'_ref.fasta')
    Path.mkdir(aligner_cwd, parents=True, exist_ok=True)
    shutil.copy2(ref_fasta_path, aligner_cwd)
    if aligner == 'bowtie2':
        index_cmd = ['bowtie2-build', '--threads', task.threads, task.id+'_ref.fasta', task.id+'_ref']
    elif aligner == 'bwa':
        index_cmd = ['bwa', 'index', '-p', task.id+'_ref', task.id+'_ref.fasta']
    logger.info('CMD: '+' '.join(index_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(index_cmd))
    ref_index_run = subprocess.run(index_cmd, cwd=aligner_cwd, capture_output=True)
    print(ref_index_run.stdout.decode(encoding='utf-8'))
    print(ref_index_run.stderr.decode(encoding='utf-8'))


def align_bowtie2(task):
    logger.info('Running Bowtie2 alignment.')
    aligner_cwd = task.path.joinpath(task.id, 'alignment', 'bowtie2')
    ref_index_path = str(aligner_cwd.joinpath(task.id+'_ref'))
    if task.dehost != None:
        filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R1.fastq.gz'))
        filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R2.fastq.gz'))
    else:
        filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz'))
        filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz'))
    reads_cmd = ['-1', filterd_R1, '-2', filterd_R2]
    thread_cmd = ['-p', str(task.threads)]
    output_cmd = ['-S', task.id+'.sam']
    other_cmd = ['--very-sensitive-local',
        '--un-conc-gz',
        '%s' % (task.id + '_unmapped_R%.fastq.gz')
    ]
    aln_cmd = ['bowtie2', '-x', ref_index_path] + reads_cmd + output_cmd + thread_cmd + other_cmd
    logger.info('CMD: '+' '.join(aln_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(aln_cmd))
    bt2_run = subprocess.run(aln_cmd, cwd=aligner_cwd, capture_output=True)
    print(bt2_run.stdout.decode(encoding='utf-8'))
    print(bt2_run.stderr.decode(encoding='utf-8'))


def align_bwa(task):
    logger.info('Running BWA alignment.')
    aligner_cwd = task.path.joinpath(task.id, 'alignment', 'bwa')
    ref_index_path = str(aligner_cwd.joinpath(task.id+'_ref'))
    if task.dehost != None:
        filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R1.fastq.gz'))
        filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R2.fastq.gz'))
    else:
        filterd_R1 = str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz'))
        filterd_R2 = str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz'))
    reads_cmd = [filterd_R1, filterd_R2]
    thread_cmd = ['-t', str(task.threads)]
    output_cmd = ['-o', task.id+'.sam']
    aln_cmd = ['bwa', 'mem'] + thread_cmd + [ref_index_path] + reads_cmd + output_cmd
    logger.info('CMD: '+' '.join(aln_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(aln_cmd))
    bt2_run = subprocess.run(aln_cmd, cwd=aligner_cwd, capture_output=True)
    print(bt2_run.stdout.decode(encoding='utf-8'))
    print(bt2_run.stderr.decode(encoding='utf-8'))


def bam_sort_n_index(task, aligner):
    logger.info('Sorting & indexing BAM file.')
    aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
    # sorting
    sorting_cmd = ['samtools', 'sort', '-@', task.threads, task.id+'.sam', '-o', task.id+'.sorted.bam']
    logger.info('CMD: '+' '.join(sorting_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(sorting_cmd))
    sorting_run = subprocess.run(sorting_cmd, cwd=aligner_cwd, capture_output=True)
    print(sorting_run.stdout.decode(encoding='utf-8'))
    print(sorting_run.stderr.decode(encoding='utf-8'))
    # indexing
    indexing_cmd = ['samtools', 'index', '-@', task.threads, task.id+'.sorted.bam']
    logger.info('CMD: '+' '.join(indexing_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(indexing_cmd))
    index_run = subprocess.run(indexing_cmd, cwd=aligner_cwd, capture_output=True)
    print(index_run.stdout.decode(encoding='utf-8'))
    print(index_run.stderr.decode(encoding='utf-8'))


def align_flagstat(task, aligners):
    stats_dict = {'mapped_rate':{}}
    for aligner in aligners:
        logger.info('Analysis BAM file from %s' % aligner)
        aligner_cwd = task.path.joinpath(task.id, 'alignment', aligner)
        flagstat_cmd = ['samtools', 'flagstat', '-@', task.threads, task.id+'.sorted.bam']
        logger.info('CMD: '+' '.join(flagstat_cmd))
        utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(flagstat_cmd))
        flagstat_run = subprocess.run(flagstat_cmd, cwd=aligner_cwd, capture_output=True)
        stats_text = flagstat_run.stdout.decode(encoding='utf-8')
        stats_list = stats_text.split('\n')
        utils.build_text_file(task.path.joinpath(aligner_cwd, 'flagstat.txt'), stats_text)
        mapped_rate = stats_list[4].split(' ')[4][1:]
        stats_dict['mapped_rate'][aligner]= mapped_rate
    utils.build_json_file(task.path.joinpath(task.id, 'alignment', 'flagstat.json'), stats_dict)


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
        bam_sort_n_index(task, aligner)
    align_flagstat(task, aligners)
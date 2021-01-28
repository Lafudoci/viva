import logging
import subprocess
import utils
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def ref_index(task_name, base_path, aligner):
    logger.info('Building Bowtie2 ref index.')
    aligner_cwd = base_path.joinpath(task_name, 'alignment', aligner)
    ref_fasta_path = base_path.joinpath(
        task_name, 'reference', task_name+'_ref.fasta')
    Path.mkdir(aligner_cwd, parents=True, exist_ok=True)
    shutil.copy2(ref_fasta_path, aligner_cwd)
    if aligner == 'bowtie2':
        index_cmd = ['bowtie2-build', task_name+'_ref.fasta', task_name+'_ref']
    elif aligner == 'bwa':
        index_cmd = ['bwa', 'index', '-p', task_name+'_ref', task_name+'_ref.fasta']
    logger.info('CMD: '+' '.join(index_cmd))
    utils.write_log_file(base_path.joinpath(task_name), 'CMD: '+' '.join(index_cmd))
    ref_index_run = subprocess.run(index_cmd, cwd=aligner_cwd, capture_output=True)
    print(ref_index_run.stdout.decode(encoding='utf-8'))
    print(ref_index_run.stderr.decode(encoding='utf-8'))


def align_bowtie2(task_name, base_path):
    logger.info('Running Bowtie2 alignment.')
    aligner_cwd = base_path.joinpath(task_name, 'alignment', 'bowtie2')
    ref_index_path = str(aligner_cwd.joinpath(task_name+'_ref'))
    filterd_R1 = str(base_path.joinpath(task_name, 'reads', task_name + '_R1.fastq.gz'))
    filterd_R2 = str(base_path.joinpath(task_name, 'reads', task_name + '_R2.fastq.gz'))
    reads_cmd = ['-1', filterd_R1, '-2', filterd_R2]
    thread_cmd = ['-p', '6']
    output_cmd = ['-S', task_name+'.sam']
    other_cmd = ['--very-sensitive-local']
    aln_cmd = ['bowtie2', '-x', ref_index_path] + reads_cmd + output_cmd + thread_cmd + other_cmd
    logger.info('CMD: '+' '.join(aln_cmd))
    utils.write_log_file(base_path.joinpath(task_name), 'CMD: '+' '.join(aln_cmd))
    bt2_run = subprocess.run(aln_cmd, cwd=aligner_cwd, capture_output=True)
    print(bt2_run.stdout.decode(encoding='utf-8'))
    print(bt2_run.stderr.decode(encoding='utf-8'))


def align_bwa(task_name, base_path):
    logger.info('Running BWA alignment.')
    aligner_cwd = base_path.joinpath(task_name, 'alignment', 'bwa')
    ref_index_path = str(aligner_cwd.joinpath(task_name+'_ref'))
    filterd_R1 = str(base_path.joinpath(task_name, 'reads', task_name + '_R1.fastq.gz'))
    filterd_R2 = str(base_path.joinpath(task_name, 'reads', task_name + '_R2.fastq.gz'))
    reads_cmd = [filterd_R1, filterd_R2]
    thread_cmd = ['-t', '6']
    output_cmd = ['-o', task_name+'.sam']
    aln_cmd = ['bwa', 'mem'] + thread_cmd + [ref_index_path] + reads_cmd + output_cmd
    logger.info('CMD: '+' '.join(aln_cmd))
    utils.write_log_file(base_path.joinpath(task_name), 'CMD: '+' '.join(aln_cmd))
    bt2_run = subprocess.run(aln_cmd, cwd=aligner_cwd, capture_output=True)
    print(bt2_run.stdout.decode(encoding='utf-8'))
    print(bt2_run.stderr.decode(encoding='utf-8'))


def bam_sort_n_index(task_name, base_path, aligner):
    logger.info('Sorting & indexing BAM file.')
    aligner_cwd = base_path.joinpath(task_name, 'alignment', aligner)
    # sorting
    sorting_cmd = ['samtools', 'sort', task_name+'.sam', '-o', task_name+'.sorted.bam']
    logger.info('CMD: '+' '.join(sorting_cmd))
    utils.write_log_file(base_path.joinpath(task_name), 'CMD: '+' '.join(sorting_cmd))
    sorting_run = subprocess.run(sorting_cmd, cwd=aligner_cwd, capture_output=True)
    print(sorting_run.stdout.decode(encoding='utf-8'))
    print(sorting_run.stderr.decode(encoding='utf-8'))
    # indexing
    indexing_cmd = ['samtools', 'index', task_name+'.sorted.bam']
    logger.info('CMD: '+' '.join(indexing_cmd))
    utils.write_log_file(base_path.joinpath(task_name), 'CMD: '+' '.join(indexing_cmd))
    index_run = subprocess.run(indexing_cmd, cwd=aligner_cwd, capture_output=True)
    print(index_run.stdout.decode(encoding='utf-8'))
    print(index_run.stderr.decode(encoding='utf-8'))


def align_flagstat(task_name, base_path, aligners):
    stats_dict = {}
    for aligner in aligners:
        logger.info('Analysis BAM file from %s' % aligner)
        aligner_cwd = base_path.joinpath(task_name, 'alignment', aligner)
        flagstat_cmd = ['samtools', 'flagstat', task_name+'.sorted.bam']
        logger.info('CMD: '+' '.join(flagstat_cmd))
        utils.write_log_file(base_path.joinpath(task_name), 'CMD: '+' '.join(flagstat_cmd))
        flagstat_run = subprocess.run(flagstat_cmd, cwd=aligner_cwd, capture_output=True)
        stats_text = flagstat_run.stdout.decode(encoding='utf-8')
        stats_list = stats_text.split('\n')
        utils.build_text_file(base_path.joinpath(aligner_cwd, 'flagstat.txt'), stats_text)
        mapping_rate = stats_list[4].split(' ')[4][1:]
        stats_dict['mapping_rate']= {aligner : mapping_rate}
    utils.build_json_file(base_path.joinpath('alignment'), stats_dict)


def align_disp(task_name, base_path, aligner):
    if aligner == 'bowtie2':
        align_bowtie2(task_name, base_path)
    elif aligner == 'bwa':
        align_bwa(task_name, base_path)

def run(task_name, base_path):
    aligners = ['bowtie2', 'bwa']
    for aligner in aligners:
        ref_index(task_name, base_path, aligner)
        align_disp(task_name, base_path, aligner)
        bam_sort_n_index(task_name, base_path, aligner)
    align_flagstat(task_name, base_path, aligners)

if __name__ == "__main__":
    run('test', Path.cwd())
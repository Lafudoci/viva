import utils
import logging
import subprocess
import sys
import shutil
import os
from decimal import Decimal
from pathlib import Path


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def import_impurities_fasta(task):
    logger.info('Importing impurities sequence')
    imported_impurities_path = task.path.joinpath(
        task.id, 'impurities_prefilter')
    Path.mkdir(imported_impurities_path, parents=True, exist_ok=True)
    impurities_fasta_dict = utils.load_fasta_file(task.remove_impurities)
    if impurities_fasta_dict == -1:
        logger.critical('Failed to load FASTA.')
        sys.exit(-1)
    task.impurities_prefilter_num = len(impurities_fasta_dict)
    if task.impurities_prefilter_num == 1:
        pass
    else:
        logger.warning('Impurities file contains muiltple squences.')

    # import impurities
    meta_dict = {
        'seq_meta': {},
        'origin_file_path': str(task.remove_impurities),
        'impurities_num': str(task.impurities_prefilter_num)
    }
    i = 1
    for header, seq in impurities_fasta_dict.items():
        imported_impurities_fasta_dict = {}
        # collect seqs fasta
        imported_impurities_fasta_path = task.path.joinpath(
            task.id,
            'impurities_prefilter',
            '%s_impurities_%d.fasta' % (task.id, i)
        )
        imported_impurities_fasta_dict[header] = seq
        # collect seqs meta
        meta_dict['seq_meta'][i] = {
            'fasta_header': header,
            'fasta_header_escape': header.replace("|", "&#124;"),
            'seq_length': str(len(seq))
        }
        i += 1
        # build impurities fasta
        utils.build_fasta_file(
            imported_impurities_fasta_path,
            imported_impurities_fasta_dict
        )
    # build meta dict
    impurities_prefilter_meta_path = task.path.joinpath(
        task.id,
        'impurities_prefilter',
        'impurities_prefilter_meta.json')
    utils.build_json_file(impurities_prefilter_meta_path, meta_dict)


def build_impurities_index(task):
    for impurities_order in range(1, task.impurities_prefilter_num+1):
        logger.info(
            'Building Bowtie2 index for impurities prefilter #%d.' % impurities_order)
        aligner_cwd = task.path.joinpath(
            task.id, 'impurities_prefilter', "bt2_alignment")
        impurities_fasta_path = task.path.joinpath(
            task.id, 'impurities_prefilter', '%s_impurities_%d.fasta' % (task.id, impurities_order))
        Path.mkdir(aligner_cwd, parents=True, exist_ok=True)
        shutil.copy2(impurities_fasta_path, aligner_cwd)
        index_cmd = ['bowtie2-build', '--threads', task.threads,
                     impurities_fasta_path.name, impurities_fasta_path.stem]
        logger.info('CMD: '+' '.join(index_cmd))
        utils.write_log_file(task.path.joinpath(task.id),
                             'CMD: '+' '.join(index_cmd))
        impurities_index_run = subprocess.run(
            index_cmd, cwd=aligner_cwd, capture_output=True)
        print(impurities_index_run.stdout.decode(encoding='utf-8'))
        print(impurities_index_run.stderr.decode(encoding='utf-8'))


def remove_impurities(task):
    aligner_cwd = task.path.joinpath(
        task.id, 'impurities_prefilter', "bt2_alignment")
    impurities_remove_meta = {}
    for impurities_order in range(1, task.impurities_prefilter_num+1):
        logger.info('Removing impurities #%d' % impurities_order)
        impurities_remove_meta[impurities_order] = {
            'mapped_reads': "", 'remove_percentage': ""}
        unconc_reads_out = task.id + '_%d' % impurities_order + \
            '_impurity_removed_R%.fastq.gz'
        mapped_sam = 'impurity_%d_mapped.sam' % impurities_order
        if impurities_order == 1:
            if task.remove_host != None:
                filterd_R1 = str(task.path.joinpath(
                    task.id, 'reads', task.id + '_host_removed_R1.fastq.gz'))
                filterd_R2 = str(task.path.joinpath(
                    task.id, 'reads', task.id + '_host_removed_R2.fastq.gz'))
            else:
                filterd_R1 = str(task.path.joinpath(
                    task.id, 'reads', task.id + '_R1.fastq.gz'))
                filterd_R2 = str(task.path.joinpath(
                    task.id, 'reads', task.id + '_R2.fastq.gz'))
        else:
            filterd_R1 = str(task.path.joinpath(
                task.id, 'reads', '%s_%d_impurity_removed_R1.fastq.gz' % (task.id, impurities_order-1)))
            filterd_R2 = str(task.path.joinpath(
                task.id, 'reads', '%s_%d_impurity_removed_R2.fastq.gz' % (task.id, impurities_order-1)))
        align_cmd = [
            'bowtie2',
            '-p', str(task.threads),
            '-x', str('%s_impurities_%d' % (task.id, impurities_order)),
            '-1', filterd_R1,
            '-2', filterd_R2,
            '-S', str(mapped_sam),
            '--very-sensitive-local',
            '--un-conc-gz', '%s' % str(unconc_reads_out)
        ]
        logger.info('CMD: '+' '.join(align_cmd))
        utils.write_log_file(task.path.joinpath(task.id),
                             'CMD: '+' '.join(align_cmd))
        cmd_run = subprocess.run(
            align_cmd, cwd=aligner_cwd, capture_output=True)
        print(cmd_run.stderr.decode(encoding='utf-8'))

        # build meta
        logger.info('Analysis BAM file from impurity mapped reads')
        # sorting
        sorted_bam = 'impurity_%d_mapped.sorted.bam' % impurities_order
        sorting_cmd = ['samtools', 'sort', '-@', task.threads,
                       mapped_sam, '-o', sorted_bam]
        logger.info('CMD: '+' '.join(sorting_cmd))
        utils.write_log_file(task.path.joinpath(task.id),
                             'CMD: '+' '.join(sorting_cmd))
        sorting_run = subprocess.run(
            sorting_cmd, cwd=aligner_cwd, capture_output=True)
        print(sorting_run.stdout.decode(encoding='utf-8'))
        print(sorting_run.stderr.decode(encoding='utf-8'))
        # flagstat
        flagstat_cmd = ['samtools', 'flagstat', '-@',
                        task.threads, sorted_bam]
        logger.info('CMD: '+' '.join(flagstat_cmd))
        utils.write_log_file(task.path.joinpath(task.id),
                             'CMD: '+' '.join(flagstat_cmd))
        flagstat_run = subprocess.run(
            flagstat_cmd, cwd=aligner_cwd, capture_output=True)
        stats_text = flagstat_run.stdout.decode(encoding='utf-8')
        stats_list = stats_text.split('\n')
        utils.build_text_file(task.path.joinpath(
            aligner_cwd, 'flagstat.txt'), stats_text)
        total_reads = task.total_reads_after_fastp
        mapped_reads = stats_list[4].split(' ')[0]
        mapped_rate = Decimal(mapped_reads)/Decimal(total_reads)
        impurities_remove_meta[impurities_order]['mapped_reads'] = mapped_reads
        impurities_remove_meta[impurities_order]['remove_percentage'] = "%f%%" % (
            mapped_rate*Decimal('100'))
        utils.build_json_file(task.path.joinpath(
            task.id, 'impurities_prefilter', 'impurities_remove.json'), impurities_remove_meta)
        # remove sam file to release disk space
        os.remove(task.path.joinpath(aligner_cwd, mapped_sam))
        # remove host bam file to release disk space
        os.remove(task.path.joinpath(aligner_cwd, sorted_bam))


def run(task):
    if task.remove_impurities != None:
        logger.info('Running impurities pre-filter.')
        import_impurities_fasta(task)
        build_impurities_index(task)
        remove_impurities(task)
    else:
        logger.info('No impurities pre-filter.')

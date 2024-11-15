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
        # Bowtie2 indexing
        for aligner in ('bt2', 'bwa'):
            logger.info(
                'Building %s index for impurities prefilter #%d.' % (aligner, impurities_order))
            if aligner == 'bt2':
                aligner_cwd = task.path.joinpath(
                    task.id, 'impurities_prefilter', "bt2_alignment")
            else:
                aligner_cwd = task.path.joinpath(
                    task.id, 'impurities_prefilter', "bwa_alignment")
            impurities_fasta_path = task.path.joinpath(
                task.id, 'impurities_prefilter', '%s_impurities_%d.fasta' % (task.id, impurities_order))
            Path.mkdir(aligner_cwd, parents=True, exist_ok=True)
            shutil.copy2(impurities_fasta_path, aligner_cwd)
            if aligner == 'bt2':
                index_cmd = ['bowtie2-build', '--threads', task.threads,
                             impurities_fasta_path.name, impurities_fasta_path.stem]
            else:
                index_cmd = ['bwa', 'index', '-p',
                             impurities_fasta_path.stem, impurities_fasta_path.name]
            logger.info('CMD: '+' '.join(index_cmd))
            utils.write_log_file(task.path.joinpath(task.id),
                                 'CMD: '+' '.join(index_cmd))
            impurities_index_run = subprocess.run(
                index_cmd, cwd=aligner_cwd, capture_output=True)
            print(impurities_index_run.stdout.decode(encoding='utf-8'))
            print(impurities_index_run.stderr.decode(encoding='utf-8'))


def remove_impurities(task):
    impurities_remove_meta = {}
    # loop impurities
    for impurities_order in range(1, task.impurities_prefilter_num+1):
        impurities_remove_meta[impurities_order] = {}
        logger.info('Removing impurities #%d' % impurities_order)
        # loop aligners
        for aligner in ('bt2', 'bwa'):
            if aligner == 'bt2':
                aligner_cwd = task.path.joinpath(
                    task.id, 'impurities_prefilter', "bt2_alignment")
            else:
                aligner_cwd = task.path.joinpath(
                    task.id, 'impurities_prefilter', "bwa_alignment")

            impurities_remove_meta[impurities_order][aligner] = {
                'mapped_reads': "", 'remove_percentage': ""
            }

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
            if aligner == 'bt2':
                align_cmd = [
                    'bowtie2',
                    '-p', str(task.threads),
                    '-x', str('%s_impurities_%d' %
                              (task.id, impurities_order)),
                    '-1', filterd_R1,
                    '-2', filterd_R2,
                    '-S', str(mapped_sam),
                    '--very-sensitive-local'
                ]
            else:
                align_cmd = [
                    'bwa',
                    'mem',
                    '-K 100000000',
                    '-Y',
                    '-t', str(task.threads),
                    str('%s_impurities_%d' % (task.id, impurities_order)),
                    filterd_R1,
                    filterd_R2,
                    '-o', str(mapped_sam)
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
            # extract unmapped with bt2
            if aligner == 'bt2':
                unmapped_fastq_r1 = str(task.path.joinpath(
                    task.id, 'reads', task.id + '_%d' % impurities_order +
                    '_impurity_removed_R1.fastq.gz'))
                unmapped_fastq_r2 = str(task.path.joinpath(
                    task.id, 'reads', task.id + '_%d' % impurities_order +
                    '_impurity_removed_R2.fastq.gz'))
                samtools_option_cmd = ['samtools', 'fastq', '-f 13']
                samtools_fastq_cmd = [
                    '-1', unmapped_fastq_r1, '-2', unmapped_fastq_r2]
                samtools_run_cmd = samtools_option_cmd + \
                    samtools_fastq_cmd + [sorted_bam]
                subprocess.run(samtools_run_cmd, cwd=aligner_cwd, check=True)
                utils.write_log_file(task.path.joinpath(
                    task.id), 'CMD: '+' '.join(samtools_run_cmd))
            # flagstat
            flagstat_cmd = ['samtools', 'flagstat', '-@',
                            task.threads, sorted_bam]
            logger.info('CMD: '+' '.join(flagstat_cmd))
            utils.write_log_file(task.path.joinpath(task.id),
                                 'CMD: '+' '.join(flagstat_cmd))
            flagstat_run = subprocess.run(
                flagstat_cmd, cwd=aligner_cwd, capture_output=True)
            stats_text = flagstat_run.stdout.decode(encoding='utf-8')
            flagstat_file_path = task.path.joinpath(
                aligner_cwd, 'flagstat_impurities_%d.txt' % impurities_order)
            utils.build_text_file(flagstat_file_path, stats_text)

            total_reads = task.total_reads_after_fastp
            primary_mapped_reads = utils.primary_mapped_from_flagstat(
                flagstat_file_path)
            mapped_rate = Decimal(primary_mapped_reads)/Decimal(total_reads)
            impurities_remove_meta[impurities_order][aligner]['mapped_reads'] = primary_mapped_reads
            impurities_remove_meta[impurities_order][aligner]['remove_percentage'] = "%f%%" % (
                mapped_rate*Decimal('100'))
            utils.build_json_file(task.path.joinpath(
                task.id, 'impurities_prefilter', 'impurities_remove.json'), impurities_remove_meta)
            # remove sam file to release disk space
            os.remove(task.path.joinpath(aligner_cwd, mapped_sam))
            # remove host bam file to release disk space
            # os.remove(task.path.joinpath(aligner_cwd, sorted_bam))


def run(task):
    if task.remove_impurities != None:
        logger.info('Running impurities pre-filter.')
        import_impurities_fasta(task)
        build_impurities_index(task)
        remove_impurities(task)
    else:
        logger.info('No impurities pre-filter.')

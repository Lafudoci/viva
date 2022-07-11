import argparse
import configparser
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import reads_alignment
import unmapped_analysis
import reads_preprocess
import reference_prepare
import utils
import variant_calling
import report_generator
import summary_generator


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class Task:
    def __init__(self):
        self.name = ''


def check_reads_file(task):
    if task.ex_r1 != None and task.ex_r2 != None:
        if Path(task.ex_r1).is_file() and Path(task.ex_r2).is_file():
            return 1
        else:
            logger.error('Reads file not found.')
            return -1
    else:
        logger.error('Reads file path can not be empty.')
        return -1


def check_ref_file(task):
    if Path(task.ref).is_file():
        return True
    else:
        logger.info('Reference sequence file not found.')
        return False


def check_deps(task):
    sys_deps = ['wget', 'git', 'apt', 'conda', 'python', 'gzip']
    if utils.deps_check(task.conda_pkgs+sys_deps) == -1:
        logger.critical('Depency check fail.')
        sys.exit(100)


def main(input_args):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ex_r1', help="Read-R1.")
    parser.add_argument(
        '--ex_r2', help="Read-R2.")
    parser.add_argument(
        '--prefix', help="For output prefix.", default='newtask')
    parser.add_argument(
        '--ref', help="Reference FASTA file path.", default=None)
    parser.add_argument(
        '--threads', help="CPU threads.", default=6)
    parser.add_argument(
        '--alns', help="Reads mapper list.", default='bowtie2,bwa')
    parser.add_argument(
        '--global_trimming', help="Global trimming bases for reads.", default=0)
    parser.add_argument(
        '--remove_host', help="Remove specific host genome (human, dog, vero, chicken, rhesus_monkey).", default=None)
    parser.add_argument(
        '--test', default=None)
    parser.add_argument(
        '--spades_mem', help="The memory (GB) allocated for spades, apply to both ref and unmapped assemble.", default=22)
    parser.add_argument(
        '--spades_mode', default='metaviral')
    parser.add_argument(
        '--unmapped_spades_mode', default='meta')
    parser.add_argument(
        '--min_vc_score', default=1)
    parser.add_argument(
        '--vc_threshold', default='0.7')
    parser.add_argument(
        '--unmapped_assemble', help="De novo Assemble the unmapped reads via metaSPAdes. ONLY apply to the first ref alignment.", default=True)
    parser.add_argument(
        '--unmapped_blastdb', help="BLASTDB name in unmapped reads assemble BLAST.", default=None)
    parser.add_argument(
        '--unmapped_len_filter', help="Min. length (bp) filter to hit in unmapped reads assemble BLAST.", default='500')
    parser.add_argument(
        '--unmapped_ident_filter', help="Min. identity (%) filter to hit in unmapped reads assemble BLAST.", default='95')
    parser.add_argument(
        '--preset_path', help="Load VIVA analysis setting from given preset file path.", default=None)
    parser.add_argument(
        '--task_note', help="Task note. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_product_name', help="Sample (product) name. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_product_lot', help="Sample (product) lot. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_sequencing_date', help="Sample sequencing date. Anotation purpose only.", default=None)
    parser.add_argument(
        '--sample_note', help="Sample note. Anotation purpose only.", default=None)
    args, unknown = parser.parse_known_args(input_args)

    task = Task()
    task.conda_pkgs = [
        'fastp', 'samtools', 'bcftools',
        'bowtie2', 'bwa',
        'varscan', 'lofreq',
        'spades.py', 'blastn', 'makeblastdb'
    ]
    check_deps(task)
    task.path = Path.cwd().joinpath('tasks')
    task.name = args.prefix
    task.task_note = args.task_note
    task.id = ''
    task.with_ref = False
    task.ex_r1 = args.ex_r1
    task.ex_r2 = args.ex_r2
    task.alns = args.alns.split(',')
    task.ref_num = 0
    task.preset_path = args.preset_path
    task.sample_product_name = args.sample_product_name
    task.sample_product_lot = args.sample_product_lot
    task.sample_sequencing_date = args.sample_sequencing_date
    task.sample_note = args.sample_note
    if task.preset_path == None:
        task.ref = args.ref
        task.threads = str(args.threads)
        task.global_trimming = str(args.global_trimming)
        task.remove_host = args.remove_host
        task.spades_mem = str(args.spades_mem)
        task.spades_mode = args.spades_mode
        task.vc_threshold = args.vc_threshold
        task.min_vc_score = args.min_vc_score
        task.unmapped_assemble = args.unmapped_assemble
        task.unmapped_spades_mode = args.unmapped_spades_mode
        task.unmapped_blastdb = args.unmapped_blastdb
        task.unmapped_len_filter = args.unmapped_len_filter
        task.unmapped_ident_filter = args.unmapped_ident_filter
    else:
        config = configparser.ConfigParser(allow_no_value=True)
        config.read(args.preset_path)
        task.ref = config['PRESET']['ref']
        task.threads = str(config['PRESET']['threads'])
        task.global_trimming = str(config['PRESET']['global_trimming'])
        task.remove_host = config['PRESET']['remove_host']
        task.spades_mem = str(config['PRESET']['spades_mem'])
        task.spades_mode = config['PRESET']['spades_mode']
        task.vc_threshold = config['PRESET']['vc_threshold']
        task.min_vc_score = config['PRESET']['min_vc_score']
        task.unmapped_assemble = config['PRESET']['unmapped_assemble']
        task.unmapped_spades_mode = config['PRESET']['unmapped_spades_mode']
        task.unmapped_blastdb = config['PRESET']['unmapped_blastdb']
        task.unmapped_len_filter = config['PRESET']['unmapped_len_filter']
        task.unmapped_ident_filter = config['PRESET']['unmapped_ident_filter']
        task.preset_id = config['VERSION']['preset_id']
        task.preset_version = config['VERSION']['version']
        task.preset_last_rev_date = config['VERSION']['last_rev_date']
        task.preset_author = config['VERSION']['author']
        task.preset_note = config['VERSION']['note']

    if args.test != None:
        task.name = 'test_run'
        task.ex_r1 = Path.cwd().joinpath('test_data', 'AdV_R1.fastq.gz')
        task.ex_r2 = Path.cwd().joinpath('test_data', 'AdV_R2.fastq.gz')
        if args.test == 'ref':
            task.ref = Path.cwd().joinpath('test_data', 'AC_000008.1.fasta')
        elif args.test == 'multi_ref':
            task.ref = Path.cwd().joinpath('test_data', 'adv_multi_ref.fasta')
        elif args.test == 'denovo':
            task.remove_host = 'human'
            task.ref = None

    if task.unmapped_blastdb != None:
        task.unmapped_assemble = True

    logger.info('Checking reference.')
    if task.ref != None:
        if check_ref_file(task):
            task.with_ref = True
        else:
            logger.error('Input reference not found. Exiting pipeline.')
            sys.exit()
    else:
        logger.info('Input reference not provided. Will go de novo')

    if task.with_ref == False:
        logger.info('Checking RVDB.')
        if utils.setup_rvdb() == -1:
            logger.error('RVDB setup error. Exiting pipeline.')
            sys.exit()

    if task.remove_host != None:
        if utils.setup_genomes(task.remove_host) == -1:
            logger.error('Host genome not found. Exiting pipeline.')
            sys.exit()

    logger.info('Checking reads files.')
    if check_reads_file(task) != -1:
        task.id = "%s_%s" % (task.name, time.strftime(
            "%Y%m%d%H%M", time.localtime()))
        logger.info('Creating new task %s.' % task.id)
        Path.mkdir(task.path.joinpath(task.id), parents=True)
        logger.info('Starting pipeline.')
        utils.write_log_file(
            task.path.joinpath(task.id),
            'Starting pipeline.'
        )

        # main pipeline
        reads_preprocess.run(task)
        reference_prepare.run(task)
        reads_alignment.run(task)
        unmapped_analysis.run(task)
        variant_calling.run(task)
        logger.info('Pipeline finished.')
        utils.write_log_file(
            task.path.joinpath(task.id),
            'Pipeline finished.'
        )

        # report generator
        summary_generator.run(task)
        report_generator.run(task)

        return task.id

    else:
        logger.error('Reads not found. Exiting pipeline.')
        sys.exit()


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])

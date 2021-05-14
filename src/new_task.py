import argparse
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import reads_alignment
import reads_preprocess
import reference_prepare
import utils
import variant_calling
import report_generator
import summary_generator

parser = argparse.ArgumentParser()
parser.add_argument('--r1', help="Read-R1.")
parser.add_argument('--r2', help="Read-R2.")
parser.add_argument('--prefix', help="For output prefix.", default='newtask')
parser.add_argument('--ref', help="Reference FASTA file path.", default=None)
parser.add_argument('--threads', help="CPU threads.", default=6)
parser.add_argument('--alns', help="Reads mapper list", default='bowtie2,bwa')
parser.add_argument('--trimming', help="Global trimming bases for reads.", default=0)
parser.add_argument('--remove_host', help="Remove specific host genome.", default=None, choices=['human', 'dog', 'vero', 'chicken', 'rhesus_monkey'])
parser.add_argument('--test', default=None)
parser.add_argument('--spades_mem', default=22)
parser.add_argument('--spades_mode', default='metaviral')
parser.add_argument('--min_vc_score', default=1)
args = parser.parse_args()

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class Task:
    def __init__(self):
        self.name = ''


def check_reads_file(task):
    if Path(task.ex_r1).is_file() and Path(task.ex_r2).is_file():
        return 1
    else:
        logger.error('Reads file not found.')
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


def main():
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
    task.id = ''
    task.ref = args.ref
    task.with_ref = False
    task.ex_r1 = args.r1
    task.ex_r2 = args.r2
    task.threads = str(args.threads)
    task.alns = args.alns.split(',')
    task.global_trimming = str(args.trimming)
    task.dehost = args.remove_host
    task.spades_mem = str(args.spades_mem)
    task.spades_mode = args.spades_mode
    task.vc_threshold = '0.7'
    task.ref_num = 0
    task.min_vc_score = args.min_vc_score
    
    if args.test != None:
        task.name = 'test_run'
        task.ex_r1 = Path.cwd().joinpath('test_data','AdV_R1.fastq.gz')
        task.ex_r2 = Path.cwd().joinpath('test_data','AdV_R2.fastq.gz')
        if args.test == 'ref':
            task.ref = Path.cwd().joinpath('test_data', 'AC_000008.1.fasta')
        elif args.test == 'multi_ref':
            task.ref = Path.cwd().joinpath('test_data', 'adv_multi_ref.fasta')
        elif args.test == 'denovo':
            task.dehost = 'human'
            task.ref = None

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
            logger.error('Reads not found. Exiting pipeline.')
            sys.exit()

    if task.dehost != None:
        if utils.setup_genomes(task.dehost) == -1:
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
        variant_calling.run(task)
        logger.info('Pipeline finished.')
        utils.write_log_file(
            task.path.joinpath(task.id),
            'Pipeline finished.'
        )

        # report generator
        summary_generator.run(task)
        report_generator.run(task)
        
    else:
        logger.error('Reads not found. Exiting pipeline.')
        sys.exit()


if __name__ == "__main__":
    main()

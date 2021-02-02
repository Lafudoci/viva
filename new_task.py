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
parser.add_argument('--r1', help="Read-R1.", default='CoV2_R1.fastq.gz')
parser.add_argument('--r2', help="Read-R2.", default='CoV2_R2.fastq.gz')
parser.add_argument('--prefix', help="For output prefix.", default='newtask')
parser.add_argument('--ref', help="Reference FASTA file path.",
                    default='NC_045512.fasta')
parser.add_argument('--threads', help="CPU threads.", default=1)
parser.add_argument('--trimming', help="Global trimming bases for reads.", default=0)
args = parser.parse_args()

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class Task:
    def __init__(self):
        self.name = ''
        self.id = ''
        self.path = ''
        self.with_ref = False
        self.ex_r1 = ''
        self.ex_r2 = ''
        self.ref_path = ''
        self.threads = 1
        self.alns = []
        self.global_trimming = 0


def check_reads_file(task):
    if Path(task.ex_r1).is_file() and Path(task.ex_r2).is_file():
        return 1
    else:
        logger.error('Reads file not found.')
        return -1


def check_ref_file(task):
    if Path(task.ref).is_file():
        return 1
    else:
        logger.info('Reference sequence file not found.')
        return -1


def main():
    task = Task()
    task.path = Path.cwd()
    task.name = args.prefix
    task.id = ''
    task.ref = args.ref
    task.ex_r1 = args.r1
    task.ex_r2 = args.r2
    task.threads = args.threads
    task.alns = ['bowtie2', 'bwa']
    task.global_trimming = args.trimming
    logger.info('Checking input file.')
    if task.ref != None:
        if check_ref_file(task):
            task.with_ref = True
    if check_reads_file(task) != -1:
        task.id = "%s_%s" % (task.name, time.strftime(
            "%Y%m%d%H%M", time.localtime()))
        logger.info('Creating new task %s.' % task.id)
        Path.mkdir(task.path.joinpath(task.id))
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
        summary_generator.run(task)
        report_generator.run(task)
        logger.info('Pipeline finished.')
        utils.write_log_file(
            task.path.joinpath(task.id),
            'Pipeline finished.'
        )
    else:
        logger.error('Reads not found. Exiting pipeline.')
        sys.exit()


if __name__ == "__main__":
    main()

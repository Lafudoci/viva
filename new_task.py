import argparse
import logging
import os
import sys
import subprocess
import time
from pathlib import Path

import reads_preprocess
import reference_prepare
import reads_alignment
import variant_calling

parser = argparse.ArgumentParser()
parser.add_argument('--r1', help="Read-R1.", default='CoV2_R1.fastq.gz')
parser.add_argument('--r2', help="Read-R2.", default='CoV2_R2.fastq.gz')
parser.add_argument('--prefix', help="For output prefix.", default='newtask')
parser.add_argument('--ref', help="Reference FASTA file path.", default='NC_045512.fasta')
args = parser.parse_args()

task_prefix = args.prefix
base_path = Path.cwd()
r1_path = args.r1
r2_path = args.r2
ref_path = args.ref

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def check_reads_file():
    if Path(r1_path).is_file() and Path(r2_path).is_file():
        return 1
    else:
        logger.error('Reads file not found.')
        return -1


def check_ref_file():
    if Path(ref_path).is_file():
        return 1
    else:
        logger.info('Reference sequence file not found.')
        return -1


def main():
    task_with_ref = False
    logger.info('Start pipeline.')
    if ref_path != None:
        if check_ref_file():
            task_with_ref = True
    if check_reads_file() != -1:
        task_name = "%s_%s" % (task_prefix, time.strftime(
            "%Y%m%d%H%M", time.localtime()))
        logger.info('Creating new task %s.' % task_name)
        Path.mkdir(base_path.joinpath(task_name))

        # main pipeline
        reads_preprocess.run(task_name, base_path, Path(r1_path), Path(r2_path))
        reference_prepare.run(task_name, base_path, Path(ref_path), task_with_ref)
        reads_alignment.run(task_name, base_path)
        variant_calling.run(task_name, base_path)
    else:
        logger.error('Exiting pipeline.')
        sys.exit()


if __name__ == "__main__":
    main()

import logging
import shutil
import subprocess
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def import_reads(task_name, base_path, external_reads_R1, external_reads_R2):
    original_reads_path = base_path.joinpath(task_name, 'reads', 'original')
    Path.mkdir(original_reads_path, parents=True, exist_ok=True)
    if str(external_reads_R1).endswith('.fastq.gz') and str(external_reads_R2).endswith('.fastq.gz'):
        shutil.copy(external_reads_R1, original_reads_path.joinpath(
            task_name + '_R1.fastq.gz'))
        shutil.copy(external_reads_R2, original_reads_path.joinpath(
            task_name + '_R2.fastq.gz'))
        reads_meta = {'r1': str(external_reads_R1), 'r2':str(external_reads_R2)}
        reads_meta_path = base_path.joinpath(task_name, 'reads', 'reads_meta.json')
        utils.build_json_file(reads_meta_path, reads_meta)
    else:
        logger.error('Unsupported reads file type.')


def run_fastp(task_name, base_path):
    logger.info('Running fastp to filter reads.')
    original_R1 = base_path.joinpath(
        task_name, 'reads', 'original', task_name + '_R1.fastq.gz')
    original_R2 = base_path.joinpath(
        task_name, 'reads', 'original', task_name + '_R2.fastq.gz')
    filterd_R1 = base_path.joinpath(
        task_name, 'reads', task_name + '_R1.fastq.gz')
    filterd_R2 = base_path.joinpath(
        task_name, 'reads', task_name + '_R2.fastq.gz')
    report_json = base_path.joinpath(task_name, 'reads', 'fastp.json')
    report_html = base_path.joinpath(task_name, 'reads', 'fastp.html')
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
        '-f', '1',
        '-t', '1',
        '-F', '1',
        '-T', '1',
        '-w', '4'       # thread
    ]
    logger.info('CMD: '+' '.join(fastp_cmd + reports_cmd + parameter_cmd))
    utils.write_log_file(
        base_path.joinpath(task_name),
        'CMD: '+' '.join(fastp_cmd + reports_cmd + parameter_cmd)
    )
    fastp_run = subprocess.run(
        fastp_cmd + reports_cmd + parameter_cmd,
        capture_output=True
    )
    print(fastp_run.stderr.decode(encoding='utf-8'))


def export_qc_report():
    pass


def run(task_name, base_path, external_reads_R1, external_reads_R2):
    logger.info('Importing reads.')
    import_reads(task_name, base_path, external_reads_R1, external_reads_R2)
    run_fastp(task_name, base_path)


if __name__ == "__main__":
    base_path = Path.cwd()
    run('test_20210128',
        Path.cwd(),
        Path('C:\\Users\\robot\\Documents\\CoV2_R1.fastq.gz'),
        Path('C:\\Users\\robot\\Documents\\CoV2_R2.fastq.gz')
        )

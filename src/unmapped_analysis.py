import logging
import subprocess
import sys
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_de_novo(task):
    logger.info('Runing De Novo assembly of unmapped reads.')
    unmapped_assembly_cwd = task.path.joinpath(task.id, 'assembly', 'unmapped')
    Path.mkdir(unmapped_assembly_cwd, parents=True, exist_ok=True)
    bwa_aligner_cwd = task.path.joinpath(task.id, 'alignment', 'bwa')
    # ONLY apply to the first bwa ref alignment.
    r1 = str(task.path.joinpath(bwa_aligner_cwd, task.id+'_ref_1_unmapped_R1.fastq.gz'))
    r2 = str(task.path.joinpath(bwa_aligner_cwd, task.id+'_ref_1_unmapped_R2.fastq.gz'))
    
    assemble_cmd = [
        'spades.py',
        '-t', task.threads,
        '-m', task.spades_mem,
        '--' + task.unmapped_spades_mode,
        '-1', r1,
        '-2', r2,
        '-o', '%s_unmapped_spades_%s'%(task.id, task.unmapped_spades_mode)
    ]
    logger.info('CMD: '+' '.join(assemble_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(assemble_cmd))
    cmd_run = subprocess.run(assemble_cmd, cwd=unmapped_assembly_cwd, capture_output=True)
    print(cmd_run.stdout.decode(encoding='utf-8'))
    print(cmd_run.stderr.decode(encoding='utf-8'))


def blast_assembled(task):
    logger.info('BLASTing unmapped reads assembled.')
    assembled_cwd = task.path.joinpath(task.id, 'assembly', 'unmapped', '%s_unmapped_spades_%s'%(task.id, task.unmapped_spades_mode))
    blast_cmd = [
        'blastn',
        '-db',
        task.unmapped_blastdb,
        '-query',
        'contigs.fasta',
        '-out',
        '%s_spades_%s.tsv'%(task.id, task.unmapped_spades_mode),
        '-outfmt',
        '6 qseqid sacc pident qlen length evalue stitle',
        '-num_threads',
        task.threads,
        '-max_target_seqs',
        '5'
    ]
    logger.info('CMD: '+' '.join(blast_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(blast_cmd))
    cmd_run = subprocess.run(blast_cmd, cwd=assembled_cwd, capture_output=True)
    print(cmd_run.stdout.decode(encoding='utf-8'))
    print(cmd_run.stderr.decode(encoding='utf-8'))


def run(task):
    run_de_novo(task)
    if task.unmapped_blastdb != None:
        blast_assembled(task)
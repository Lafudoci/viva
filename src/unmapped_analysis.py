import logging
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_de_novo(task):
    logger.info('Runing De Novo assembly of unmapped reads.')
    unmapped_assembly_cwd = task.path.joinpath(task.id, 'unmapped_analysis')
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
    assembled_cwd = task.path.joinpath(task.id, 'unmapped_analysis', '%s_unmapped_spades_%s'%(task.id, task.unmapped_spades_mode))
    qeury_filename = 'contigs.fasta'
    blast_result_filename = '%s_spades_%s.tsv'%(task.id, task.unmapped_spades_mode)
    blast_cmd = [
        'blastn',
        '-db',
        task.unmapped_blastdb,
        '-query',
        qeury_filename,
        '-out',
        blast_result_filename,
        '-outfmt',
        '6 qseqid sacc pident qlen length evalue stitle',
        '-num_threads',
        task.threads,
        '-max_target_seqs',
        '5',
        '-evalue',
        '1e-15',
    ]
    logger.info('CMD: '+' '.join(blast_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(blast_cmd))
    cmd_run = subprocess.run(blast_cmd, cwd=assembled_cwd, capture_output=True)
    print(cmd_run.stdout.decode(encoding='utf-8'))
    print(cmd_run.stderr.decode(encoding='utf-8'))
    
    # filter highly matched hits
    logger.info('Filter highly matched hits')
    blast_result_path = assembled_cwd.joinpath(blast_result_filename)
    highly_match_result_list = []
    with open(blast_result_path, 'r') as f:
        for line in f.readlines():
            hit = line.strip().split('\t')
            # identity
            if Decimal(hit[2]) <= Decimal('95'):
                break
            if task.unmapped_blastdb == "U-RVDBv21.0.fasta":
                clean_sacc = hit[6].split('|')[2]
                clean_stitle = hit[6].split('|')[3]
                clean_stitle_org = hit[6].split('|')[4]
            else:
                clean_sacc = hit[1]
                clean_stitle = hit[6]
                clean_stitle_org = hit[6]

            highly_match_result_list.append({
                'qseqid': hit[0],
                'sacc': hit[1],
                'pident': hit[2],
                'qlen': hit[3],
                'length': hit[4],
                'evalue': hit[5],
                'stitle': hit[6],
                'clean_sacc': clean_sacc,
                'clean_stitle': clean_stitle,
                'clean_stitle_org': clean_stitle_org
            })

    unmapped_analysis = {
        'spades_mode': task.unmapped_spades_mode,
        'BLASTdb_name': task.unmapped_blastdb,
        'highly_matched_result': highly_match_result_list
    }
    
    unmapped_analysis_json_path = task.path.joinpath(task.id, 'unmapped_analysis', 'unmapped_analysis.json')
    utils.build_json_file(unmapped_analysis_json_path, unmapped_analysis)


def run(task):
    run_de_novo(task)
    if task.unmapped_blastdb != None:
        blast_assembled(task)
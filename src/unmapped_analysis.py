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
        '--phred-offset', '33',
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
        '6 qseqid sacc pident qlen length evalue stitle bitscore',
        '-num_threads',
        task.threads,
        '-max_target_seqs',
        '100',
        '-max_hsps',
        '1',
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
    filtered_hits_list = []
    with open(blast_result_path, 'r') as f:
        for line in f.readlines():
            hit = line.strip().split('\t')
            if blast_hits_significant_filter(task, hit):
                filtered_hits = blast_hits_string_formater(task, hit)
                filtered_hits_list.append(filtered_hits)
    
    highly_match_result_list = blast_hits_max1_bitscore_filter(task, filtered_hits_list)

    unmapped_analysis = {
        'spades_mode': task.unmapped_spades_mode,
        'BLASTdb_name': task.unmapped_blastdb,
        'highly_matched_result': highly_match_result_list
    }
    
    unmapped_analysis_json_path = task.path.joinpath(task.id, 'unmapped_analysis', 'unmapped_analysis.json')
    utils.build_json_file(unmapped_analysis_json_path, unmapped_analysis)


def blast_hits_max1_bitscore_filter(task, hits_list):
    filtered_dict = {}
    filtered_list = []
    for hit in hits_list:
        if not hit['qseqid'] in filtered_dict:
            filtered_dict[hit['qseqid']] = {
                'current_best_score': Decimal('0'),
                'hit': {}
            }
        if Decimal(hit['bitscore']) > Decimal(filtered_dict[hit['qseqid']]['current_best_score']):
            filtered_dict[hit['qseqid']]['hit'] = hit.copy()
            filtered_dict[hit['qseqid']]['current_best_score'] = str(hit['bitscore'])
    for v in filtered_dict.values():
        filtered_list.append(v['hit'])
    return filtered_list


def blast_hits_significant_filter(task, hit):
    if Decimal(hit[2]) < Decimal(task.unmapped_ident_filter):
        return False
    if Decimal(hit[4]) < Decimal(task.unmapped_len_filter):
        return False
    else:
        return True


def blast_hits_string_formater(task, hit):
    if task.unmapped_blastdb == "U-RVDBv21.0.fasta":
        hit_split_list = hit[6].split('|')
        clean_sacc = hit_split_list[2]
        clean_stitle = hit_split_list[3]
        if len(hit_split_list) >= 5:
            clean_stitle_org = hit_split_list[4]
        else:
            clean_stitle_org = ''
    else:
        clean_sacc = hit[1]
        clean_stitle = hit[6]
        clean_stitle_org = ''
    return {
        'qseqid': hit[0],
        'sacc': hit[1],
        'pident': hit[2],
        'qlen': hit[3],
        'length': hit[4],
        'evalue': hit[5],
        'stitle': hit[6],
        'bitscore': hit[7],
        'clean_sacc': clean_sacc,
        'clean_stitle': clean_stitle,
        'clean_stitle_org': clean_stitle_org
    }

def run(task):
    if task.unmapped_assemble == True:
        run_de_novo(task)
        if task.unmapped_blastdb != None:
            blast_assembled(task)
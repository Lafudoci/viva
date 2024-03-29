import logging
import subprocess
import sys
from pathlib import Path

import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_de_novo(task):
    logger.info('Runing De Novo assembly.')
    assembly_cwd = task.path.joinpath(task.id, 'assembly')
    Path.mkdir(assembly_cwd, parents=True, exist_ok=True)
    if task.remove_host != None:
        r1 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R1.fastq.gz'))
        r2 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R2.fastq.gz'))
    else:
        r1 = str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz'))
        r2 = str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz'))
    
    assemble_cmd = [
        'spades.py',
        '-t', task.threads,
        '-m', task.spades_mem,
        '--' + task.spades_mode,
        '--phred-offset', '33',
        '-1', r1,
        '-2', r2,
        '-o', '%s_spades_%s'%(task.id, task.spades_mode)
    ]
    logger.info('CMD: '+' '.join(assemble_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(assemble_cmd))
    cmd_run = subprocess.run(assemble_cmd, cwd=assembly_cwd, capture_output=True)
    print(cmd_run.stdout.decode(encoding='utf-8'))
    print(cmd_run.stderr.decode(encoding='utf-8'))


def blast_assembled(task):
    logger.info('Finding ref seqeunce.')
    assembled_cwd = task.path.joinpath(task.id, 'assembly', '%s_spades_%s'%(task.id, task.spades_mode))
    blast_cmd = [
        'blastn',
        '-db',
        task.unmapped_blastdb,
        '-query',
        'contigs.fasta',
        '-out',
        '%s_spades_%s.tsv'%(task.id, task.spades_mode),
        '-outfmt',
        '6',
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


def extract_virus_refseq(task):
    spades_prefix = '%s_spades_%s'%(task.id, task.spades_mode)
    fmt6_path = task.path.joinpath(task.id, 'assembly', spades_prefix, spades_prefix+'.tsv')
    fmt6_dict = utils.load_blast_fmt6_max1_bitscore(fmt6_path)
    if len(fmt6_dict) > 0:
        best_hit_dict = utils.find_top_score_hits(fmt6_dict)
        utils.build_json_file(task.path.joinpath(task.id, 'assembly', 'best_hit.json'), best_hit_dict)
        refseq_virus_fasta_path = Path('/app/blastdb/%s'%task.unmapped_blastdb)
        fasta_dict = utils.extract_seq_from_fasta(refseq_virus_fasta_path, best_hit_dict['sseqid'])
        task.ref = task.path.joinpath(task.id, 'assembly', '%s.fasta'%best_hit_dict['sseqid'])
        utils.build_fasta_file(task.ref, fasta_dict)
    else:
        logger.critical('De novo result did not provide viable virus reference.')
        sys.exit(-1)

def ref_import(task):
    logger.info('Importing reference sequence')
    imported_ref_path = task.path.joinpath(task.id, 'reference')
    Path.mkdir(imported_ref_path, parents=True, exist_ok=True)
    ref_fasta_dict = utils.load_fasta_file(task.ref)
    if ref_fasta_dict == -1:
        logger.critical('Failed to load FASTA.')
        sys.exit(-1)
    task.ref_num = len(ref_fasta_dict)
    if task.ref_num == 1:
        pass
    else:
        logger.warning('Reference file contains muiltple squences.')
    # import ref    
    meta_dict = {
        'ref_from_user':'',
        'seq_meta':{},
        'spades_mode':'',
        'origin_file_path': str(task.ref),
        'ref_num': str(task.ref_num)
        }
    # set spades meta
    if task.with_ref == False:
        meta_dict['spades_mode'] = task.spades_mode
        meta_dict['ref_from_user'] = 'No'
    else:
        meta_dict['spades_mode'] = 'N/A'
        meta_dict['ref_from_user'] = 'Yes'
    
    i = 1
    for header, seq in ref_fasta_dict.items():
        imported_ref_fasta_dict = {}
        # collect seqs fasta
        imported_ref_fasta_path = task.path.joinpath(
        task.id,
        'reference',
        '%s_ref_%d.fasta'%(task.id, i)
        )
        imported_ref_fasta_dict[header] = seq
        # collect seqs meta
        meta_dict['seq_meta'][i] = {
            'fasta_header': header,
            'fasta_header_escape': header.replace("|", "&#124;"),
            'seq_length': str(len(seq))
        }
        i += 1
        # build ref fasta
        utils.build_fasta_file(
            imported_ref_fasta_path,
            imported_ref_fasta_dict
        )
    # build meta dict
    imported_ref_meta_path = task.path.joinpath(
        task.id,
        'reference',
        task.id+'_ref.json')
    utils.build_json_file(imported_ref_meta_path, meta_dict)


def run(task):
    if task.with_ref:
        ref_import(task)
    else:
        run_de_novo(task)
        blast_assembled(task)
        extract_virus_refseq(task)
        ref_import(task)

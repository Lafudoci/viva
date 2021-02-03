import logging
import subprocess
import utils
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_de_novo(task):
    logger.info('Runing De Novo assembly.')
    assembly_cwd = task.path.joinpath(task.id, 'assembly')
    Path.mkdir(assembly_cwd, parents=True, exist_ok=True)
    if task.dehost == True:
        remove_host(task)
        r1 = str(task.path.joinpath(task.id, 'reads', task.id + '_R1.fastq.gz'))
        r2 = str(task.path.joinpath(task.id, 'reads', task.id + '_R2.fastq.gz'))
    else:
        r1 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R1.fastq.gz'))
        r2 = str(task.path.joinpath(task.id, 'reads', task.id + '_host_removed_R2.fastq.gz'))
    
    assemble_cmd = [
        'spades.py',
        '-t', task.threads,
        '-m', task.spades_mem,
        '--' + task.spades_mode,
        '-1', r1,
        '-2', r2,
        '-o', '%s_spades_%s'%(task.id, task.spades_mode)
    ]
    logger.info('CMD: '+' '.join(assemble_cmd))
    utils.write_log_file(task.path.joinpath(task.id), 'CMD: '+' '.join(assemble_cmd))
    cmd_run = subprocess.run(assemble_cmd, cwd=assembly_cwd, capture_output=True)
    print(cmd_run.stdout.decode(encoding='utf-8'))
    print(cmd_run.stderr.decode(encoding='utf-8'))

def remove_host(task):
    pass


def extract_refseq(task):
    pass


def external_ref_import(task):
    external_ref_path = Path(task.ref)
    logger.info('Importing reference sequence')
    imported_ref_path = task.path.joinpath(task.id, 'reference')
    Path.mkdir(imported_ref_path, parents=True, exist_ok=True)
    ref_fasta_dict = utils.load_fasta_file(external_ref_path)
    if len(ref_fasta_dict) == 1:
        # import ref
        imported_ref_fasta_dict = {}
        imported_ref_fasta_path = task.path.joinpath(
            task.id,
            'reference',
            task.id+'_ref.fasta'
        )
        first_fasta_seq = str(list(ref_fasta_dict.values())[0])
        imported_ref_fasta_dict[task.id + '_ref'] = first_fasta_seq
        utils.build_fasta_file(
            imported_ref_fasta_path,
            imported_ref_fasta_dict
        )
        # build meta file
        imported_ref_meta_path = task.path.joinpath(
            task.id, 'reference', task.id+'_ref.json')
        file_name = str(external_ref_path)
        fasta_header = str(list(ref_fasta_dict)[0])
        meta_dict = {
            'file_name': file_name,
            'fasta_header': fasta_header,
            'user_provide': True
        }
        utils.build_json_file(imported_ref_meta_path, meta_dict)
    else:
        logger.error('Reference sequence file must be single squence.')
    pass


def run(task):
    if task.with_ref:
        external_ref_import(task)
    else:
        run_de_novo(task)
        extract_refseq(task)

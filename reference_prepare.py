import logging
import subprocess
import utils
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def run_de_novo():
    pass


def find_ncbi_refseq():
    pass


def external_ref_import(task_name, base_path, external_ref_path):
    logger.info('Importing reference sequence')
    imported_ref_path = base_path.joinpath(task_name, 'reference')
    Path.mkdir(imported_ref_path, parents=True, exist_ok=True)
    ref_fasta_dict = utils.load_fasta_file(external_ref_path)
    if len(ref_fasta_dict) == 1:
        # import ref
        imported_ref_fasta_dict = {}
        imported_ref_fasta_path = base_path.joinpath(
            task_name, 'reference', task_name+'_ref.fasta')
        imported_ref_fasta_dict[task_name +
                                '_ref'] = str(list(ref_fasta_dict.values())[0])
        utils.build_fasta_file(imported_ref_fasta_path,
                               imported_ref_fasta_dict)
        # build meta file
        imported_ref_meta_path = base_path.joinpath(
            task_name, 'reference', task_name+'_ref.json')
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


def run(task_name, base_path, external_ref_path, task_with_ref):
    if task_with_ref:
        external_ref_import(task_name, base_path, external_ref_path)
    else:
        run_de_novo()
        find_ncbi_refseq()

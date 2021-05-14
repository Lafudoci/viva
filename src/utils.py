import json
import logging
import os
import subprocess
import textwrap
import time
from decimal import Decimal
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def is_valid_fasta(file_path):
    pass


def load_fasta_file(file_path):
    fasta_dict = {}
    with open(file_path, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                header = line[1:].strip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.strip()
    return fasta_dict


def build_fasta_file(file_path, fasta_dict):
    with open(file_path, 'w', encoding='utf-8') as f:
        for header, seq in fasta_dict.items():
            splited_seq_list = textwrap.wrap(seq, 80)
            f.write('>%s\n' % header)
            for line in splited_seq_list:
                f.write(line+'\n')


def load_json_file(file_path):
    with open(file_path, 'r') as f:
        j = json.load(f)
    return j


def build_json_file(file_path, python_dict):
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(python_dict, f)
    pass


def build_text_file(file_path, text):
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(text)


def load_vcf_file(file_path):
    vcf_dict = {'comments': [], 'column_names': [], 'vc': []}
    with open(file_path, 'r') as f:
        for line in f.readlines():
            if line.startswith('##'):
                vcf_dict['comments'].append(line.strip())
            elif line.startswith('#'):
                vcf_dict['column_names'] = line[1:].strip().split('\t')
            else:
                vc = line.strip().split('\t')
                if len(vc) == len(vcf_dict['column_names']):
                    vcf_dict['vc'].append(vc)
                else:
                    logger.error('Parsing VCF error.')
        # print(vcf_dict['comments'][:3])
        # print(vcf_dict['column_names'])
        # print(vcf_dict['vc'][:3])
        return vcf_dict


def write_log_file(log_path, text):
    log_file_path = log_path.joinpath('log.txt')
    with open(log_file_path, 'a') as f:
        f.write('%d\t%s\n' % (int(time.time()), text))


def load_log_file(log_path):
    log_file_path = log_path.joinpath('log.txt')
    with open(log_file_path, 'r') as f:
        log_list = f.readlines()
    return log_list


def load_blast_fmt6_max1_bitscore(file_path):
    fmt6_dict = {}
    with open(file_path, 'r') as f:
        for line in f.readlines():
            hit = line.strip().split('\t')
            if hit[0] in fmt6_dict:
                if Decimal(hit[11]) <= Decimal(fmt6_dict[hit[0]]['bitscore']):
                    continue
            fmt6_dict[hit[0]] = {
                'qseqid': hit[0],
                'sseqid': hit[1],
                'pident': hit[2],
                'length': hit[3],
                'mismatch': hit[4],
                'gapopen': hit[5],
                'qstart': hit[6],
                'qend': hit[7],
                'sstart': hit[8],
                'send': hit[9],
                'evalue': hit[10],
                'bitscore': hit[11]
            }
    return fmt6_dict


def build_fmt6_file(file_path, fmt6_dict):
    with open(file_path, 'w') as f:
        for sseqid, rest in fmt6_dict.items():
            hit_line = '%s'%sseqid
            for col in rest.values():
                hit_line += '\t'+col
            f.write(hit_line+'\n')


def find_top_score_hits(fmt6_dict):
    top_hit = {}
    for hit in fmt6_dict.values():
        if Decimal(hit['bitscore']) > Decimal(top_hit.get('bitscore', 0)):
            top_hit = hit.copy()
    return top_hit


def extract_seq_from_fasta(file_path, target_acc):
    seq_dict = {}
    extracting = False
    target_header = ''
    with open(file_path, 'r') as f:
        for line in f.readlines():
            if target_acc in line:
                extracting = True
                target_header = line.strip()[1:]
                seq_dict[target_header] = ''
                continue
            if extracting:
                if line.startswith('>'):
                    return seq_dict
                else:
                    seq_dict[target_header] += line.strip()


def load_blast_fmt_sciname_max1_bitscore(file_path):
    fmt6_dict = {}
    with open(file_path, 'r') as f:
        for line in f.readlines():
            hit = line.strip().split('\t')
            if hit[0] in fmt6_dict:
                if Decimal(hit[4]) <= Decimal(fmt6_dict[hit[0]]['bitscore']):
                    continue
            fmt6_dict[hit[0]] = {
                'sseqid': hit[1],
                'pident': hit[2],
                'length': hit[3],
                'bitscore': hit[4],
                'sci': hit[5],
                'common': hit[6]
            }
    return fmt6_dict


def deps_check(dep_list):
    try:
        for dep in dep_list:
            logger.info('Dependency check: %s'%dep)
            subprocess.run(['which', dep], check=True)
    except subprocess.CalledProcessError as e:
        logger.error('Dependency check error: %s.'%str(e))
        return -1


def conda_pkg_versions(pkg_list):
    verions_dict = {}
    if deps_check(['conda']) != -1:
        all_pkg_list = subprocess.run(['conda', 'list'], capture_output=True).stdout.decode(encoding='utf-8').split('\n')
        for pkg_string in all_pkg_list:
            print(pkg_string)
            if not pkg_string.startswith('#'):
                if len(pkg_string.split()) == 4:
                    name = pkg_string.split()[0].strip()
                    version = pkg_string.split()[1].strip()
                    if name in pkg_list:
                        verions_dict[name] = version
        print(verions_dict)
        return verions_dict
    else:
        return -1


def setup_rvdb():
    if deps_check(['wget', 'gunzip', 'makeblastdb']) == -1:
        return -1
    if subprocess.run(['blastdbcmd', '-db', 'U-RVDBv21.0.fast', '-info']).returncode == 0:
        return
    try:
        logger.info('Preparing RVDB')
        # logger.info('Downloading RVDB')
        # subprocess.run(
        #     [
        #     'wget',
        #     'https://rvdb.dbi.udel.edu/download/U-RVDBv21.0.fasta.gz',
        #     '-P', '/app/blastdb'
        #     ],
        #     check=True)
        logger.info('Decompressing RVDB')
        subprocess.run(
            [
            'gunzip',
            '--keep',
            'U-RVDBv21.0.fasta.gz',
            ],
            check=True,
            cwd='/app/blastdb')
        logger.info('Building blastdb of RVDB')
        subprocess.run(
            [
            'makeblastdb',
            '-in',
            'U-RVDBv21.0.fasta',
            '-blastdb_version',
            '5',
            '-title',
            'Reference Viral DataBase',
            '-dbtype',
            'nucl'
            ],
            check=True,
            cwd='/app/blastdb')
    except subprocess.CalledProcessError as e:
        logger.error('RVDB setup error: %s.'%str(e))
        return -1


def setup_genomes(host_name):
    if deps_check(['wget', 'gzip', 'bowtie2-inspect', 'bowtie2-build']) == -1:
        return -1
    if subprocess.run(['bowtie2-inspect', '--summary', 'grch38'], cwd='/app/genomes/grch38').returncode == 0:
        return
    try:
        logger.info('Preparing host genome file')
        # logger.info('Downloading genome file')
        # subprocess.run(
        #     [
        #     'wget',
        #     'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz',
        #     '-P', '/app/genome'
        #     ],
        #     check=True)
        # logger.info('Decompressing genome file')
        subprocess.run(
            [
            'gunzip',
            '--keep',
            'GCF_000001405.39_GRCh38.p13_genomic.fna.gz',
            ],
            check=True,
            cwd='/app/genomes/grch38')
        logger.info('Indexing genome file')
        subprocess.run(
            [
            'bowtie2-build',
            '--threads',
            '6',
            'GCF_000001405.39_GRCh38.p13_genomic.fna',
            'grch38'
            ],
            check=True,
            cwd='/app/genome')
    except subprocess.CalledProcessError as e:
        logger.error('RVDB setup error: %s.'%str(e))
        return -1

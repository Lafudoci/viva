import hashlib
import json
import logging
import os
import subprocess
import sys
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
    header = ''
    with open(file_path, 'r', encoding='UTF-8') as f:
        for line in f.readlines():
            if '>' in line:
                header = line.rstrip().split('>')[1]
                if len(header) > 0:
                    fasta_dict[header] = ''
            else:
                # only write base after first ">" symbol
                if header != '':
                    fasta_dict[header] += line.rstrip()
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


def load_rvdb_anno_tab(file_path):
    annotation_data = {}
    try:
        with open(file_path, 'r') as infile:
            for line in infile:
                accession, header_info, start_str, end_str, category = line.strip().split('\t')
                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    logger.warning("Skipping line due to invalid start/end values %s"%accession)
                    continue

                if accession not in annotation_data:
                    annotation_data[accession] = []
                annotation_data[accession].append({
                    'start': start,
                    'end': end,
                    'category': category
                })
    except FileNotFoundError:
        logger.error("Error: Annotation file not found")
        return {}
    return annotation_data

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
            hit_line = '%s' % sseqid
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


def sys_deps_check(dep_list):
    try:
        for dep in dep_list:
            logger.info('Dependency check: %s' % dep)
            subprocess.run(['which', dep], check=True)
    except subprocess.CalledProcessError as e:
        logger.error('Dependency check error: %s.' % str(e))
        return -1


def conda_deps_check(dep_list):
    verions_dict = conda_pkg_versions(dep_list)
    if verions_dict != -1:
        for dep in dep_list:
            logger.info('Dependency check: %s' % dep)
            if dep not in verions_dict:
                logger.error('Dependency check error: %s.' % dep)
                return -1
    else:
        return -1


def conda_pkg_versions(pkg_list):
    verions_dict = {}
    if sys_deps_check(['conda']) != -1:
        all_pkg_list = subprocess.run(['conda', 'list', '-p', sys.prefix], capture_output=True).stdout.decode(
            encoding='utf-8').split('\n')
        for pkg_string in all_pkg_list:
            # print(pkg_string)
            if not pkg_string.startswith('#'):
                if len(pkg_string.split()) >= 3:
                    name = pkg_string.split()[0].strip()
                    version = pkg_string.split()[1].strip()
                    if name in pkg_list:
                        verions_dict[name] = version
        # print(verions_dict)
        return verions_dict
    else:
        return -1


def md5_check(file_path, md5_string):
    try:
        logger.info('Checking md5 hash.')
        hashmd5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hashmd5.update(chunk)
        hashmd5_string = hashmd5.hexdigest()
        if hashmd5_string != md5_string:
            logger.error('md5 hash mismatched\nDownload: %s\n Expected:%s\n')
            return -1
    except subprocess.CalledProcessError as e:
        logger.error('md5 hash error: %s.' % str(e))
        return -1


def setup_blastdb(blastdb_path, blastdb_name):
    auto_rvdb_fasta_ver_list = ['U-RVDBv29.0.fasta', 'C-RVDBv29.0.fasta']
    auto_rvdb_fastagz_md5_dict = {
        'C-RVDBv29.0.fasta.gz': 'deb369751ea32c723f640ee192688e48',
        'U-RVDBv29.0.fasta.gz': 'e6352c74dc691a600e830bceca650c3a'
    }
    try:
        if sys_deps_check(['wget', 'gunzip', 'makeblastdb']) == -1:
            return -1
        if blastdb_path != None:
            # use coustom blastdb
            # copy n modify BLASDB env
            m_env = os.environ.copy()
            m_env['BLASTDB'] = blastdb_path
            if subprocess.run(['blastdbcmd', '-db', blastdb_name, '-info'], env=m_env).returncode == 0:
                logger.info('blastdb %s at %s exists.' %
                            (blastdb_name, blastdb_path))
                return
            else:
                logger.info('blastdb %s at %s not found.' %
                            (blastdb_name, blastdb_path))
                return -1
        else:
            # check built-in app/blastdb
            if subprocess.run(['blastdbcmd', '-db', blastdb_name, '-info']).returncode == 0:
                logger.info('blastdb %s at app/blastdb exists.' %
                            (blastdb_name))
                return
            else:
                logger.info('blastdb %s at app/blastdb not found.' %
                            (blastdb_name))
                # if use rvdb then go setup, else then exit
                if blastdb_name in auto_rvdb_fasta_ver_list:
                    if Path("/app/blastdb_arch/%s.gz" % blastdb_name).is_file():
                        logger.info('blastdb archive gz exists.')
                    else:
                        download_rvdb(blastdb_name)
                    # check md5 hash
                    if md5_check(Path("/app/blastdb_arch/%s.gz" % blastdb_name),
                                 auto_rvdb_fastagz_md5_dict[blastdb_name+'.gz']) == -1:
                        return -1
                    # decompress
                    decompress_rvdb(blastdb_name)
                    # build blastdb
                    Path.mkdir(Path("/app/blastdb"),
                               parents=True, exist_ok=True)
                    logger.info('Building blastdb')
                    subprocess.run(
                        [
                            'makeblastdb',
                            '-in',
                            blastdb_name,
                            '-blastdb_version',
                            '5',
                            '-title',
                            'Reference Viral DataBase (%s)' % blastdb_name,
                            '-dbtype',
                            'nucl'
                        ],
                        check=True,
                        cwd='/app/blastdb')
                else:
                    logger.error('%s not found in app/blastdb' % blastdb_name)
                    return -1
    except subprocess.CalledProcessError as e:
        logger.error('blastdb setup error: %s.' % str(e))
        return -1
def download_rvdb(blastdb_name):
    logger.info('Preparing RVDB')
    rvdb_fasta = blastdb_name
    rvdb_fastagz = blastdb_name + '.gz'
    try:
        logger.info('Downloading RVDB %s' % blastdb_name)
        subprocess.run(
            [
                'wget',
                'https://rvdb.dbi.udel.edu/download/%s' % rvdb_fastagz,
                '-P', '/app/blastdb_arch'
            ],
            check=True)
    except subprocess.CalledProcessError as e:
        logger.error('RVDB setup error: %s.' % str(e))
        return -1


def decompress_rvdb(blastdb_name):
    rvdb_fastagz = blastdb_name + '.gz'
    try:
        logger.info('Decompressing RVDB')
        with open("/app/blastdb/%s" % blastdb_name, "w") as f:
            subprocess.run(
                [
                    'gunzip',
                    '-c',
                    '/app/blastdb_arch/%s' % rvdb_fastagz,
                ],
                check=True,
                cwd='/app/blastdb',
                stdout=f
            )
    except subprocess.CalledProcessError as e:
        logger.error('RVDB setup error: %s.' % str(e))
        return -1



def check_genome_availability(host_file_name, genome_source_dir):
    try:
        if sys_deps_check(['bowtie2-inspect']) == -1:
            return -1
        genome_index_prefix = Path("/app/genomes").joinpath(host_file_name)
        inspect_cmd = ['bowtie2-inspect', '--summary', str(genome_index_prefix)]
        if subprocess.run(inspect_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
            return 0
        
        if genome_source_dir is None:
            logger.error('Genome index not found and genome_path is not provided.')
            return -1
            
        source_path = Path(genome_source_dir).joinpath(host_file_name)
        if not source_path.is_file():
            logger.error('Host genome source file %s not found in %s.' % (host_file_name, genome_source_dir))
            return -1
        return 0
    except Exception as e:
        logger.error('Genome check error: %s.' % str(e))
        return -1


def setup_genomes(host_file_name, genome_source_dir):
    try:
        if sys_deps_check(['bowtie2-inspect', 'bowtie2-build']) == -1:
            return -1

        # 索引路徑與前綴（直接使用檔名作為前綴）
        genome_index_prefix = Path("/app/genomes").joinpath(host_file_name)
        
        # 1. 檢查索引是否已存在 (使用 bowtie2-inspect 驗證完整性)
        inspect_cmd = ['bowtie2-inspect', '--summary', str(genome_index_prefix)]
        if subprocess.run(inspect_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
            logger.info('Host genome index for %s already exists and is valid.' % host_file_name)
            return 0

        # 2. 若索引不存在，則需要從來源建立
        if genome_source_dir is None:
            logger.error('Genome index not found and genome_path is not provided.')
            return -1
            
        source_path = Path(genome_source_dir).joinpath(host_file_name)
        if not source_path.is_file():
            logger.error('Host genome source file %s not found in %s.' % (host_file_name, genome_source_dir))
            return -1

        logger.info('Preparing genome index for %s from %s' % (host_file_name, source_path))
        Path.mkdir(Path("/app/genomes"), parents=True, exist_ok=True)
        
        # 處理是否需要解壓縮 (建立暫存檔)
        temp_fasta_path = Path("/app/genomes").joinpath(host_file_name + ".temp.fna")
        is_gz = host_file_name.endswith('.gz')
        
        if is_gz:
            logger.info('Decompressing genome file to temporary FASTA')
            with open(temp_fasta_path, "w") as f:
                subprocess.run(
                    ['gunzip', '-c', str(source_path)],
                    check=True,
                    stdout=f
                )
            build_in_file = str(temp_fasta_path)
        else:
            build_in_file = str(source_path)
        
        # 建立 Bowtie2 索引
        logger.info('Indexing genome file (this may take a while)...')
        subprocess.run(
            [
                'bowtie2-build',
                '--threads', '6',
                build_in_file,
                str(genome_index_prefix)
            ],
            check=True,
            cwd='/app/genomes'
        )
        
        # 清除暫存檔
        if is_gz and temp_fasta_path.exists():
            os.remove(temp_fasta_path)
            
        return 0

    except subprocess.CalledProcessError as e:
        logger.error('Genome setup error (subprocess): %s.' % str(e))
        return -1
    except Exception as e:
        logger.error('Genome setup unexpected error: %s.' % str(e))
        return -1


def primary_mapped_from_flagstat(file_path):
    with open(file_path, 'r', encoding='UTF-8') as f:
        for line in f.readlines():
            if 'primary mapped' in line:
                primary_mapped_reads = line.split(' ')[0]
                return primary_mapped_reads
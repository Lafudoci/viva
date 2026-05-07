import hashlib
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
        all_pkg_list = subprocess.run(['conda', 'list'], capture_output=True).stdout.decode(
            encoding='utf-8').split('\n')
        for pkg_string in all_pkg_list:
            # print(pkg_string)
            if not pkg_string.startswith('#'):
                if len(pkg_string.split()) == 4:
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
    if blastdb_name == None:
        return -1
    
    # 只要是 RVDB 格式都支援自動設定
    is_rvdb = blastdb_name.startswith('U-RVDB') or blastdb_name.startswith('C-RVDB')
    
    auto_rvdb_fastagz_md5_dict = {
        'C-RVDBv29.0.fasta.gz': 'deb369751ea32c723f640ee192688e48',
        'U-RVDBv29.0.fasta.gz': 'e6352c74dc691a600e830bceca650c3a'
    }
    
    try:
        if sys_deps_check(['wget', 'gunzip', 'makeblastdb']) == -1:
            return -1
            
        m_env = os.environ.copy()
        if blastdb_path:
            m_env['BLASTDB'] = blastdb_path
        else:
            m_env['BLASTDB'] = "/app/blastdb:" + m_env.get('BLASTDB', '')

        # 1. 檢查是否已經建立索引
        if subprocess.run(['blastdbcmd', '-db', blastdb_name, '-info'], env=m_env, capture_output=True).returncode == 0:
            logger.info(f'blastdb {blastdb_name} is ready and indexed.')
            return
            
        # 2. 如果沒索引，嘗試尋找原始檔進行建立
        search_dir = Path(blastdb_path) if blastdb_path else Path("/app/blastdb")
        fasta_file = search_dir / blastdb_name
        gz_file = search_dir / (blastdb_name + '.gz')
        
        # 如果是預設路徑且檔案不存在，嘗試處理 RVDB 下載
        if not blastdb_path and not fasta_file.exists() and not gz_file.exists():
            if is_rvdb:
                arch_gz = Path("/app/blastdb_arch") / (blastdb_name + ".gz")
                if not arch_gz.is_file():
                    download_rvdb(blastdb_name)
                
                expected_md5 = auto_rvdb_fastagz_md5_dict.get(blastdb_name + '.gz')
                if expected_md5 and md5_check(arch_gz, expected_md5) == -1:
                    return -1
                
                Path.mkdir(search_dir, parents=True, exist_ok=True)
                with open(fasta_file, 'w') as f:
                    subprocess.run(['gunzip', '-c', str(arch_gz)], stdout=f)
            else:
                logger.error(f"BLASTDB {blastdb_name} not found in {search_dir}")
                return -1
        
        # 如果只有 .gz，先解壓
        if not fasta_file.exists() and gz_file.exists():
            logger.info(f"Decompressing {gz_file}...")
            # 解壓到當前目錄
            subprocess.run(['gunzip', '-k', str(gz_file)], check=True)
            
        # 建立索引
        if fasta_file.exists():
            logger.info(f"Building index for {blastdb_name}...")
            subprocess.run([
                'makeblastdb', '-in', str(fasta_file),
                '-dbtype', 'nucl',
                '-out', str(search_dir / blastdb_name),
                '-blastdb_version', '5'
            ], check=True)
            return
        else:
            logger.error(f"Cannot find fasta file to build index for {blastdb_name}")
            return -1

    except Exception as e:
        logger.error(f'Error setting up blastdb: {e}')
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


def setup_genomes(host_name):
    genome_id = 'GRCh38.p14'
    genome_source_table = {
        'human': {
            'bt2_gname': genome_id,
            'arch_name': 'GCF_000001405.40_GRCh38.p14_genomic.fna.gz',
            'file_name': 'GCF_000001405.40_GRCh38.p14_genomic.fna',
            'source_url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz',
            'md5': 'c30471567037b2b2389d43c908c653e1'
        }
    }
    # TEMP return 0 for custom genome input
    
    if sys_deps_check(['wget', 'gzip', 'bowtie2-inspect', 'bowtie2-build']) == -1:
        return -1

    # Custom genome input
    if genome_source_table.get(host_name) == None:
        try:
            # Check if index exists
            if subprocess.run(['bowtie2-inspect', '--summary', host_name], cwd='/app/genomes', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0:
                return 0
            
            # Check if archive exists
            # We assume user provides filename like 'genome.fna.gz' and puts it in /app/genomes_arch/
            arch_path = Path("/app/genomes_arch").joinpath(host_name)
            if not arch_path.is_file():
                # Provide a more specific error message if it looks like a file but is missing
                logger.error('Custom genome index not found and archive %s not found in /app/genomes_arch/.' % host_name)
                return -1

            logger.info('Preparing custom genome index for %s' % host_name)
            Path.mkdir(Path("/app/genomes"), parents=True, exist_ok=True)
            
            # Decompress
            temp_fasta_name = host_name
            if host_name.endswith('.gz'):
                temp_fasta_name = host_name[:-3]
            
            temp_fasta_path = Path("/app/genomes").joinpath(temp_fasta_name)
            
            logger.info('Decompressing genome file')
            with open(temp_fasta_path, "w") as f:
                subprocess.run(
                    ['gunzip', '-c', host_name],
                    check=True,
                    cwd='/app/genomes_arch',
                    stdout=f
                )
            
            # Build Index
            logger.info('Indexing genome file')
            subprocess.run(
                [
                    'bowtie2-build',
                    '--threads',
                    '6',
                    str(temp_fasta_path),
                    host_name
                ],
                check=True,
                cwd='/app/genomes'
            )
            
            # Clean up temp fasta
            if temp_fasta_path.exists():
                os.remove(temp_fasta_path)
                
            return 0

        except subprocess.CalledProcessError as e:
            logger.error('Custom genome setup error: %s.' % str(e))
            return -1
        except Exception as e:
            logger.error('Custom genome setup unexpected error: %s.' % str(e))
            return -1

    try:
        if subprocess.run(['bowtie2-inspect', '--summary', genome_source_table[host_name]['bt2_gname']], cwd='/app/genomes').returncode == 0:
            return
    except Exception:
        pass
    try:
        logger.info('Preparing host genome file')
        Path.mkdir(Path("/app/genomes"), parents=True, exist_ok=True)
        if Path("/app/genomes_arch/"+genome_source_table[host_name]['arch_name']).exists() != True:
            logger.info('Downloading genome file')
            subprocess.run(
                [
                    'wget',
                    genome_source_table[host_name]['source_url'],
                    '-P', '/app/genomes_arch'
                ],
                check=True)
        else:
            logger.info('Genome archive exists.')
        # check md5 hash
        logger.info('Checking md5 hash of genome file')
        if md5_check(Path("/app/genomes_arch/%s" % genome_source_table[host_name]['arch_name']),
                        genome_source_table[host_name]['md5']) == -1:
            return -1
        logger.info('Decompressing genome file')
        with open('/app/genomes/'+genome_source_table[host_name]['file_name'], "w") as f:
            subprocess.run(
                [
                    'gunzip',
                    '-c',
                    genome_source_table[host_name]['arch_name'],
                ],
                check=True,
                cwd='/app/genomes_arch',
                stdout=f)
        logger.info('Indexing genome file')
        subprocess.run(
            [
                'bowtie2-build',
                '--threads',
                '6',
                genome_source_table[host_name]['file_name'],
                genome_id
            ],
            check=True,
            cwd='/app/genomes')
    except subprocess.CalledProcessError as e:
        logger.error('RVDB setup error: %s.' % str(e))
        return -1


def primary_mapped_from_flagstat(file_path):
    with open(file_path, 'r', encoding='UTF-8') as f:
        for line in f.readlines():
            if 'primary mapped' in line:
                primary_mapped_reads = line.split(' ')[0]
                return primary_mapped_reads
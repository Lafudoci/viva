import json
import logging
import os
import textwrap
import time
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
    vcf_dict = {'comments':[], 'column_names':[], 'vc':[]}
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
        f.write('%d\t%s\n'%(int(time.time()), text))


def load_log_file(log_path):
    log_file_path = log_path.joinpath('log.txt')
    with open(log_file_path, 'r') as f:
        log_list = f.readlines()
    return log_list


if __name__ == "__main__":
    load_vcf_file('newtask_202101260045_varscan.vcf')
import mysql.connector
import json
import os
import sys
from datetime import datetime
from pathlib import Path


def read_json(json_path):
    with open(json_path, "r") as file:
        data = json.load(file)
    return data


def read_n_upload(summary_json_path):
    print("Connected to MariaDB.")
    # Establish a connection to the MariaDB server
    connection = mysql.connector.connect(
        host=os.environ['MYSQL_HOST'],
        user=os.environ['MYSQL_USER'],
        password=os.environ['MYSQL_PASSWORD'],
        database=os.environ['VIVA_DB'],
    )
    with connection.cursor() as cursor:
        try:
            data = read_json(summary_json_path)

            if connection.is_connected():
                task_meta_cols = (
                    'create_time',
                    'task_id',
                    'task_name',
                    'task_note',
                    'start_date',
                    'finish_date',
                    'ref_source_contig_name',
                    'ref_source_contig_pident'
                )
                # Prepare insert statement
                insert = f"INSERT INTO task_meta ({','.join(task_meta_cols)}) VALUES ({','.join(['%s'] * len(task_meta_cols))})"
                # Prepare value statement
                values = []
                for col in task_meta_cols:
                    if col.endswith('create_time'):
                        values.append(datetime.now())
                    elif col.endswith('_date'):
                        dt = datetime.strptime(
                            data[col], "%Y-%m-%d %H:%M:%S UTC+8")
                        values.append(dt)
                    else:
                        values.append(data[col])
                # Execute the SQL statement
                cursor.execute(insert, values)

                reads_meta_cols = (
                    'create_time',
                    'task_id',
                    'sample_product_name',
                    'sample_product_lot',
                    'sample_sequencing_date',
                    'sample_note',
                    'reads_file_name_r1',
                    'reads_file_name_r2',
                    'reads_file_md5_r1',
                    'reads_file_md5_r2'
                )
                # Prepare insert statement
                insert = f"INSERT INTO reads_meta ({','.join(reads_meta_cols)}) VALUES ({','.join(['%s'] * len(reads_meta_cols))})"
                # Prepare value statement
                values = []
                for col in reads_meta_cols:
                    if col.endswith('create_time'):
                        values.append(datetime.now())
                    if col == 'task_id':
                        values.append(data[col])
                    if col.startswith('sample_'):
                        values.append(data['reads_meta']['sample_meta'][col])
                    if col == 'reads_file_name_r1':
                        values.append(data['reads_meta']
                                      ['reads_file_meta']['file_name']['r1'])
                    if col == 'reads_file_name_r2':
                        values.append(data['reads_meta']
                                      ['reads_file_meta']['file_name']['r2'])
                    if col == 'reads_file_md5_r1':
                        values.append(data['reads_meta']
                                      ['reads_file_meta']['md5']['r1'])
                    if col == 'reads_file_md5_r2':
                        values.append(data['reads_meta']
                                      ['reads_file_meta']['md5']['r2'])
                # Execute the SQL statement
                cursor.execute(insert, values)

                fastp_meta_cols = (
                    'create_time',
                    'task_id',
                    'before_total_reads',
                    'before_total_bases',
                    'before_total_q30',
                    'before_r1_length',
                    'before_r2_length',
                    'after_total_reads',
                    'after_total_bases',
                    'after_total_q30',
                    'after_r1_length',
                    'after_r2_length',
                    'duplication_rate'
                )
                # Prepare insert statement
                insert = f"INSERT INTO fastp_abs ({','.join(fastp_meta_cols)}) VALUES ({','.join(['%s'] * len(fastp_meta_cols))})"
                # Prepare value statement
                values = []
                for col in fastp_meta_cols:
                    if col.endswith('create_time'):
                        values.append(datetime.now())
                    elif col == 'task_id':
                        values.append(data[col])
                    else:
                        values.append(data['fastp_abs'][col])
                # Execute the SQL statement
                cursor.execute(insert, values)

                remove_genome_cols = (
                    'create_time',
                    'task_id',
                    'genome',
                    'mapped_reads',
                    'remove_percentage'
                )
                # Prepare insert statement
                insert = f"INSERT INTO remove_genome ({','.join(remove_genome_cols)}) VALUES ({','.join(['%s'] * len(remove_genome_cols))})"
                # Prepare value statement
                values = []
                for col in remove_genome_cols:
                    if col.endswith('create_time'):
                        values.append(datetime.now())
                    elif col == 'task_id':
                        values.append(data[col])
                    else:
                        values.append(data['remove_genome'][col])
                # Execute the SQL statement
                cursor.execute(insert, values)

                impurit_filter_cols = (
                    'create_time',
                    'task_id',
                    'order_no',
                    'fasta_header',
                    'fasta_header_escape',
                    'seq_length',
                    'mapped_reads',
                    'remove_percentage',
                    'origin_file_path',
                    'impurities_num'
                )
                # Prepare insert statement
                insert = f"INSERT INTO impurit_filter ({','.join(impurit_filter_cols)}) VALUES ({','.join(['%s'] * len(impurit_filter_cols))})"
                # Prepare value statement
                order_total_no = int(
                    data['impurit_filter_meta'].get('impurities_num', 0))
                order = 1
                for order in range(1, order_total_no+1):
                    values = []
                    for col in impurit_filter_cols:
                        if col.endswith('create_time'):
                            values.append(datetime.now())
                        elif col == 'task_id':
                            values.append(data[col])
                        elif col == 'order_no':
                            values.append(str(order))
                        elif col.startswith('fasta_') or col.startswith('seq_'):
                            values.append(data['impurit_filter_meta']
                                          ['seq_meta'][str(order)][col])
                        elif col == 'origin_file_path':
                            values.append(data['impurit_filter_meta'][col])
                        elif col == 'impurities_num':
                            values.append(data['impurit_filter_meta'][col])
                        elif col == 'mapped_reads':
                            values.append(
                                data['impurit_filter_results'][str(order)][col])
                        elif col == 'remove_percentage':
                            values.append(
                                data['impurit_filter_results'][str(order)][col])
                    # Execute the SQL statement
                    cursor.execute(insert, values)

                ref_aln_cols = (
                    'create_time',
                    'task_id',
                    'order_no',
                    'ref_from_user',
                    'spades_mode',
                    'origin_file_path',
                    'fasta_header',
                    'fasta_header_escape',
                    'seq_length',
                    'ref_source_contig_name',
                    'ref_source_contig_pident',
                    'mapped_rate_bt2',
                    'mapped_reads_bt2',
                    'startpos_bt2',
                    'endpos_bt2',
                    'numreads_bt2',
                    'covbases_bt2',
                    'coverage_bt2',
                    'meandepth_bt2',
                    'meanbaseq_bt2',
                    'meanmapq_bt2',
                    'mapped_rate_bwa',
                    'mapped_reads_bwa',
                    'startpos_bwa',
                    'endpos_bwa',
                    'numreads_bwa',
                    'covbases_bwa',
                    'coverage_bwa',
                    'meandepth_bwa',
                    'meanbaseq_bwa',
                    'meanmapq_bwa'
                )
                # Prepare insert statement
                insert = f"INSERT INTO ref_aln ({','.join(ref_aln_cols)}) VALUES ({','.join(['%s'] * len(ref_aln_cols))})"
                # Prepare value statement
                order_total_no = int(data['ref_meta_dict']['ref_num'])
                order = 1
                for order in range(1, order_total_no+1):
                    values = []
                    for col in ref_aln_cols:
                        if col.endswith('create_time'):
                            values.append(datetime.now())
                        elif col == 'task_id':
                            values.append(data[col])
                        elif col == 'order_no':
                            values.append(str(order))
                        elif col == 'ref_from_user':
                            values.append(data['ref_meta_dict'][col])
                        elif col == 'spades_mode':
                            values.append(data['ref_meta_dict'][col])
                        elif col == 'origin_file_path':
                            values.append(data['ref_meta_dict'][col])
                        elif col == 'fasta_header':
                            values.append(data['ref_meta_dict']
                                          ['seq_meta'][str(order)][col])
                        elif col == 'fasta_header_escape':
                            values.append(data['ref_meta_dict']
                                          ['seq_meta'][str(order)][col])
                        elif col == 'seq_length':
                            values.append(data['ref_meta_dict']
                                          ['seq_meta'][str(order)][col])
                        elif col == 'ref_source_contig_name':
                            values.append(data[col])
                        elif col == 'ref_source_contig_pident':
                            values.append(data[col])
                        elif col.endswith('_bt2'):
                            if col.startswith('mapped_rate'):
                                values.append(
                                    data['aln']['mapped_rate']['bowtie2'][str(order)])
                            elif col.startswith('mapped_reads'):
                                values.append(
                                    data['aln']['mapped_reads']['bowtie2'][str(order)])
                            else:
                                values.append(
                                    data['cov']['bowtie2'][str(order)][col.split('_')[0]])
                        elif col.endswith('_bwa'):
                            if col.startswith('mapped_rate'):
                                values.append(
                                    data['aln']['mapped_rate']['bowtie2'][str(order)])
                            elif col.startswith('mapped_reads'):
                                values.append(
                                    data['aln']['mapped_reads']['bowtie2'][str(order)])
                            else:
                                values.append(data['cov']['bwa']
                                              [str(order)][col.split('_')[0]])
                    # Execute the SQL statement
                    cursor.execute(insert, values)

                draft_meta_cols = (
                    'create_time',
                    'task_id',
                    'conflicts',
                    'snv_list',
                    'error',
                    'file_path'
                )
                # Prepare insert statement
                insert = f"INSERT INTO draft_meta ({','.join(draft_meta_cols)}) VALUES ({','.join(['%s'] * len(draft_meta_cols))})"
                # Prepare value statement
                order_total_no = int(data['ref_meta_dict']['ref_num'])
                order = 1
                for order in range(1, order_total_no+1):
                    values = []
                    for col in draft_meta_cols:
                        if col.endswith('create_time'):
                            values.append(datetime.now())
                        elif col == 'task_id':
                            values.append(data[col])
                        else:
                            values.append(
                                str(data['draft_meta'][str(order)][col]))
                    # Execute the SQL statement
                    cursor.execute(insert, values)

                unmapped_analysis_cols = (
                    'create_time',
                    'task_id',
                    'spades_mode',
                    'BLASTdb_name',
                    'qseqid',
                    'sacc',
                    'pident',
                    'qlen',
                    'length',
                    'evalue',
                    'stitle',
                    'bitscore',
                    'clean_sacc',
                    'clean_stitle',
                    'clean_stitle_org'
                )
                # Prepare insert statement
                insert = f"INSERT INTO unmapped_analysis ({','.join(unmapped_analysis_cols)}) VALUES ({','.join(['%s'] * len(unmapped_analysis_cols))})"
                # Prepare value statement
                order_total_no = len(
                    data['unmapped_analysis']['highly_matched_result'])
                if order_total_no == 0:
                    values = []
                    for col in unmapped_analysis_cols:
                        if col.endswith('create_time'):
                            values.append(datetime.now())
                        elif col == 'task_id':
                            values.append(data[col])
                        elif col == 'spades_mode':
                            values.append(data['unmapped_analysis'][col])
                        elif col == 'BLASTdb_name':
                            values.append(data['unmapped_analysis'][col])
                        else:
                            values.append('N/A')
                else:
                    order = 1
                    for order in range(1, order_total_no+1):
                        values = []
                        for col in unmapped_analysis_cols:
                            if col.endswith('create_time'):
                                values.append(datetime.now())
                            elif col == 'task_id':
                                values.append(data[col])
                            elif col == 'spades_mode':
                                values.append(data['unmapped_analysis'][col])
                            elif col == 'BLASTdb_name':
                                values.append(data['unmapped_analysis'][col])
                            else:
                                values.append(
                                    str(data['unmapped_analysis']['highly_matched_result'][order-1][col]))
                # Execute the SQL statement
                cursor.execute(insert, values)

                versions_cols = (
                    'create_time',
                    'task_id',
                    'name',
                    'version'
                )
                # Prepare insert statement
                insert = f"INSERT INTO versions ({','.join(versions_cols)}) VALUES ({','.join(['%s'] * len(versions_cols))})"
                # Prepare value statement
                version_dict = data['version']
                for name, version in version_dict.items():
                    values = []
                    for col in versions_cols:
                        if col.endswith('create_time'):
                            values.append(datetime.now())
                        elif col == 'task_id':
                            values.append(data[col])
                        elif col == 'name':
                            values.append(name)
                        elif col == 'version':
                            values.append(version)
                    # Execute the SQL statement
                    cursor.execute(insert, values)

                # Commit the changes to the database
                print("Writing to database...")
                connection.commit()
                print("Upload finished.")
            else:
                print(f"Error: Database connection error.")
                sys.exit(-1)

        except mysql.connector.Error as err:
            print(f"Error: {err}")
            connection.rollback()
            sys.exit(-1)


def run(task):
    summary_path = task.path.joinpath(task.id, task.id + '_summary.json')
    read_n_upload(summary_path)


if __name__ == "__main__":
    summary_json_path = Path(
        '/home/leftc/Documents/viva-work/tasks/test_run_202402020614/test_run_202402020614_summary.json')
    read_n_upload(summary_json_path)

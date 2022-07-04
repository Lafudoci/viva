import csv
import json
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def generate_summary_csv(batch_task_id, taks_id_list):
    summaries_dict = {}
    for task_id in taks_id_list:
        with open('/app/tasks/%s/%s_summary.json'%(task_id, task_id), 'r') as f:
            j = json.load(f)
            summaries_dict[task_id] = {
                'task date': j['start_date'],
                'product': j['reads_meta']['sample_meta']['sample_product_name'],
                'lot': j['reads_meta']['sample_meta']['sample_product_lot'],
                'seq date': j['reads_meta']['sample_meta']['sample_sequencing_date'],
                'fastp pre Q30': j['fastp_abs']['before_total_q30'],
                'fastp post Q30': j['fastp_abs']['after_total_q30'],
                'fastp duplication': j['fastp_abs']['duplication_rate'],
                'remove genome name': j['remove_genome']['genome'],
                'remove genome precent': j['remove_genome']['remove_percentage'],
                'aln rate bt2 of 1st ref': j['aln']['mapped_rate']['bowtie2']['1'],
                'aln rate bwa of 1st ref': j['aln']['mapped_rate']['bwa']['1'],
                'cov bt2 of 1st ref': j['cov']['bowtie2']['1']['coverage'],
                'cov bwa of 1st ref': j['cov']['bwa']['1']['coverage'],
                'depth bt2 of 1st ref': j['cov']['bwa']['1']['meandepth'],
                'depth bwa of 1st ref': j['cov']['bwa']['1']['meandepth']
            }

    
    with open('%s_summary.csv'%batch_task_id, 'w', newline='') as csvfile:
        fieldnames = [
        'task date',
        'product',
        'lot',
        'seq date',
        'fastp pre Q30',
        'fastp pre Q30',
        'fastp duplication',
        'remove genome name',
        'remove genome precent',
        'aln rate bt2 of 1st ref',
        'aln rate bwa of 1st ref',
        'cov bt2 of 1st ref',
        'cov bwa of 1st ref',
        'depth bt2 of 1st ref',
        'depth bwa of 1st ref'
        ]
        csv_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        csv_writer.writeheader()
        for s in summaries_dict:
            csv_writer.writerow(s)
        

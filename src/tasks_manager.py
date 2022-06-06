import argparse
import configparser
import logging
import subprocess
import sys
import time
from ast import keyword
from pathlib import Path

import new_task
import utils

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def get_reads_meta_from_ids(read_id_list, file_path):
    db_reads_meta = configparser.ConfigParser()
    db_reads_meta.read(file_path)
    selected_meta = {}
    for read_id in read_id_list:
        if read_id in db_reads_meta:
            selected_meta[read_id] = dict(db_reads_meta.items(read_id)).copy()
        else:
            logger.critical('Read ID: %s was not found in the db.')
            sys.exit(-1)
    return selected_meta


def load_sample_list_file(file_path):
    sample_list = []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            sample_id = line.strip()
            sample_list.append(sample_id)
    return sample_list


def input_paths_check(file_path):
    if Path(file_path) != None and Path(file_path).is_file():
        return True
    else:
        logger.info('Reference sequence file not found.')
        return False


def batch_task_viva(task_sheet_dict):
    queue_dict = {}
    number = 1
    for v in task_sheet_dict.values():
        read_id_list = load_sample_list_file(v['samples_list_filepath'])
        reads_meta_dict = get_reads_meta_from_ids(
            read_id_list, v['read_meta_path'])
        for read_id, read_meta in reads_meta_dict.items():
            if v['batch_tasks_note'] != '':
                batch_task_note = v['batch_tasks_note']
            else:
                batch_task_note = None
            queue_dict[number] = {
                'read_id': read_id,
                'read_meta_dict': read_meta,
                'batch_task_note': batch_task_note,
                'preset_path': v['preset_path'],
                'task_prefix': '%s-%s-%s' % (
                    read_meta['product'],
                    read_meta['lot'],
                    read_meta['seq_date'])
            }
            number += 1
    batch_task_id = "batch_task_%s" % time.strftime(
        "%Y%m%d%H%M", time.localtime())
    batch_task_queue_path = Path.cwd().joinpath('tasks').joinpath('%s_queue.json' % batch_task_id)
    utils.build_json_file(batch_task_queue_path, queue_dict)
    queue_length = len(queue_dict)
    for current_queue_number in range(1, queue_length+1):
        logger.info('Excuting queue No. %d' % current_queue_number)
        task_inputs = queue_dict[current_queue_number].copy()
        prefix = task_inputs['task_prefix']
        ex_r1 = task_inputs['read_meta_dict']['ex_r1']
        ex_r2 = task_inputs['read_meta_dict']['ex_r2']
        preset = task_inputs['preset_path']
        task_note = task_inputs['batch_task_note']
        new_task.main([
            '--prefix', prefix,
            '--ex_r1', ex_r1,
            '--ex_r2', ex_r2,
            '--preset', preset,
            '--task_note', task_note,
            '--sample_product_name', task_inputs['read_meta_dict']['product'],
            '--sample_product_lot', task_inputs['read_meta_dict']['lot'],
            '--sample_sequencing_date', task_inputs['read_meta_dict']['seq_date'],
            '--sample_note', task_inputs['read_meta_dict']['reads_note']
        ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--single_task', help="Run single task.", action='store_true')
    parser.add_argument(
        '--task_sheet', help="Task sheet filename to run with.", default=None)
    args, unknown = parser.parse_known_args()
    if args.single_task:
        new_task.main(sys.argv[1:])
    elif args.task_sheet != None:
        task_sheet_dict = {}
        task_sheet_config = configparser.ConfigParser()
        task_sheet_config.read(args.task_sheet)
        for section in task_sheet_config.sections():
            task_sheet_dict[section] = {}
            for key, val in task_sheet_config.items(section):
                task_sheet_dict[section][key] = val
        batch_task_viva(task_sheet_dict)
    else:
        logger.critical('Must provide a task sheet or use --single_task arg.')
        sys.exit(-1)


if __name__ == "__main__":
    main()

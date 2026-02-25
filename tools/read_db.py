import sqlite3
import os
import sys

import argparse

def get_db_path(args):
    if args.db_path:
        return args.db_path
    
    cwd = os.getcwd()
    # 嘗試多個可能存放 DB 的相對位置
    paths_to_try = [
        os.path.join(cwd, 'tasks', 'viva_results.db'),
        os.path.join(cwd, 'src', 'tasks', 'viva_results.db'),
        os.path.join(cwd, '..', 'tasks', 'viva_results.db')
    ]
    for p in paths_to_try:
        if os.path.exists(p):
            return p
            
    # 若都找不到，預設回傳當前目錄下的 tasks 目錄
    return os.path.join(cwd, 'tasks', 'viva_results.db')

def fetch_task_info(cursor, task_id):
    print(f"\n{'='*60}")
    print(f"Details for Task: {task_id}")
    print(f"{'='*60}")
    
    # 1. 任務主表 (Tasks)
    cursor.execute("SELECT * FROM Tasks WHERE task_id=?", (task_id,))
    task_row = cursor.fetchone()
    if task_row:
        col_names = [description[0] for description in cursor.description]
        print("\n--- 1. Task General Summary (Tasks) ---")
        for col_name, value in zip(col_names, task_row):
            print(f"  {col_name}: {value}")
            
    # 2. 品管與除宿主數據表 (QC_Metrics)
    cursor.execute("SELECT * FROM QC_Metrics WHERE task_id=?", (task_id,))
    qc_row = cursor.fetchone()
    if qc_row:
        col_names = [description[0] for description in cursor.description]
        print("\n--- 2. QC Metrics (QC_Metrics) ---")
        for col_name, value in zip(col_names, qc_row):
            if col_name != 'task_id':
                print(f"  {col_name}: {value}")
    else:
        print("\n--- 2. QC Metrics ---")
        print("  No QC metrics recorded (possibly failed before QC).")

    # 3. 比對結果數據表 (Alignment_Stats)
    cursor.execute("SELECT * FROM Alignment_Stats WHERE task_id=?", (task_id,))
    aln_rows = cursor.fetchall()
    if aln_rows:
        col_names = [description[0] for description in cursor.description]
        print(f"\n--- 3. Alignment Stats ({len(aln_rows)} records) ---")
        for i, row in enumerate(aln_rows, 1):
             print(f"  [Record {i}]:")
             for col_name, value in zip(col_names, row):
                 if col_name not in ('id', 'task_id'):
                     print(f"    {col_name}: {value}")
    else:
        print("\n--- 3. Alignment Stats ---")
        print("  No alignment data recorded.")
                     
    # 4. 變異點表 (Variants)
    cursor.execute("SELECT * FROM Variants WHERE task_id=?", (task_id,))
    var_rows = cursor.fetchall()
    if var_rows:
        col_names = [description[0] for description in cursor.description]
        print(f"\n--- 4. Variants ({len(var_rows)} records total, showing up to 10) ---")
        for i, row in enumerate(var_rows[:10], 1):
            print(f"  [Record {i}]:")
            for col_name, value in zip(col_names, row):
                 if col_name not in ('id', 'task_id'):
                     print(f"    {col_name}: {value}")
        if len(var_rows) > 10:
            print(f"  ... and {len(var_rows) - 10} more variants")
    else:
        print("\n--- 4. Variants ---")
        print("  No variants recorded.")


def main():
    parser = argparse.ArgumentParser(description="Quickly read and display VIVA database summary.")
    parser.add_argument('--db_path', '-d', type=str, help="Path to the viva_results.db file. If not provided, the script will try to find it automatically.", default=None)
    args = parser.parse_args()

    db_path = get_db_path(args)
    if not os.path.exists(db_path):
        print(f"Error: Database file not found at {db_path}")
        print("請確認您已經執行過任務並且產生了 viva_results.db。")
        sys.exit(1)
        
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # 確認 Tasks 表格是否存在
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Tasks';")
        if not cursor.fetchone():
            print("Error: The 'Tasks' table does not exist in the database.")
            sys.exit(1)
            
        # 列出所有任務
        cursor.execute("SELECT task_id, status, start_date FROM Tasks ORDER BY start_date DESC;")
        tasks = cursor.fetchall()
        
        if not tasks:
            print("No tasks found in the database. (空的表)")
            sys.exit(0)
            
        print("\nAvailable Tasks (從新到舊):")
        for i, (task_id, status, start_date) in enumerate(tasks, 1):
            if status == "Completed":
                status_color = "\033[92m" + status + "\033[0m" # Green
            elif status == "Failed":
                status_color = "\033[91m" + status + "\033[0m" # Red
            else:
                status_color = "\033[93m" + status + "\033[0m" # Yellow
                
            print(f" [{i}] {task_id} | Status: {status_color} | Started: {start_date}")
            
        while True:
            try:
                choice = input("\n請輸入想檢視的任務編號 (或輸入 'q' 退出): ")
                if choice.lower() == 'q':
                    break
                idx = int(choice) - 1
                if 0 <= idx < len(tasks):
                    selected_task_id = tasks[idx][0]
                    fetch_task_info(cursor, selected_task_id)
                else:
                    print("不正確的編號，請重新選擇。")
            except ValueError:
                print("請輸入有效的數字格式！")

    except sqlite3.Error as e:
        print(f"SQLite 錯誤: {e}")
    finally:
        if 'conn' in locals() and conn:
            conn.close()

if __name__ == '__main__':
    main()

import argparse
import os
import json
import logging
from pathlib import Path
import sys

# 將 src 目錄加入 sys.path 以便匯入 db_manager
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(project_root, 'src'))
import db_manager

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def import_task(task_dir, db):
    """將單個先前已完成的任務目錄匯入資料庫"""
    task_dir_path = Path(task_dir)
    if not task_dir_path.is_dir():
        logger.warning(f"路徑不是一個目錄, 略過: {task_dir}")
        return False

    task_id = task_dir_path.name
    summary_path = task_dir_path / f"{task_id}_summary.json"
    
    # 檢查是否有 summary.json 來確認這是一個已完成的 VIVA 任務
    if not summary_path.exists():
        logger.warning(f"找不到 summary.json ({summary_path}), 將以 'Unknown' 狀態匯入資料庫。")
        
        # 嘗試從 task_id 解析 task_name 和 start_date
        task_name = task_id
        start_date = 'Unknown'
        parts = task_id.rsplit('_', 1)
        if len(parts) == 2:
            task_name = parts[0]
            date_part = parts[1]
            if len(date_part) == 12 and date_part.isdigit():
                start_date = f"{date_part[:4]}-{date_part[4:6]}-{date_part[6:8]} {date_part[8:10]}:{date_part[10:12]}:00"
        
        db.create_task(
            task_id=task_id,
            task_name=task_name,
            start_date=start_date,
            task_note="Imported without summary.json"
        )
        db.update_task_status(task_id, 'Unknown', 'Missing summary.json')
        return True

    try:
        with open(summary_path, 'r', encoding='utf-8') as f:
            s_dict = json.load(f)
            
        task_name = s_dict.get('task_name', 'Unknown')
        start_date = s_dict.get('start_date', 'Unknown')
        preset_id = s_dict.get('preset', {}).get('preset_id') if s_dict.get('preset') else None
        task_note = s_dict.get('task_note')
        
        # 讀取樣本 metadata
        sample_meta = s_dict.get('reads_meta', {}).get('sample_meta', {})
        product = sample_meta.get('sample_product_name')
        lot = sample_meta.get('sample_product_lot')
        seq_date = sample_meta.get('sample_sequencing_date')

        logger.info(f"開始處理任務: {task_id}")
        
        # 1. 建立任務 (這會在資料庫中建立該筆任務記錄，狀態為 Running)
        logger.info(f"建立任務記錄...")
        db.create_task(
            task_id=task_id,
            task_name=task_name,
            start_date=start_date,
            preset_id=preset_id,
            task_note=task_note,
            product=product,
            lot=lot,
            seq_date=seq_date
        )
        
        # 2. 將狀態更新為 Completed
        logger.info(f"更新任務狀態為 Completed...")
        db.update_task_status(task_id, 'Completed')
        
        # 3. 匯入最終整合數據 (這會覆寫正確的 finish_date 並寫入 QC_Metrics, Alignment_Stats, Variants 等資料表)
        logger.info(f"儲存最終進階數據到各個資料表...")
        db.save_final_summary(task_id, s_dict)
        
        logger.info(f"成功匯入任務: {task_id}")
        return True

    except Exception as e:
        logger.error(f"匯入任務 {task_id} 時發生錯誤: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="批次將先前已完成的 VIVA 分析結果目錄匯入 local SQLite 資料庫(viva_results.db)")
    parser.add_argument('task_dirs', nargs='+', help="一個或多個已完成分析的任務目錄路徑。支援使用萬用字元 (如 /path/to/tasks/*)")
    parser.add_argument('-d', '--db_path', type=str, default=None, help="手動指定 viva_results.db 資料庫路徑 (預設為當前執行目錄下的 tasks/viva_results.db)")
    parser.add_argument('--rebuild', action='store_true', help="如果設定此選項，匯入前會先刪除舊的資料庫檔案重建 Schema")
    
    args = parser.parse_args()
    
    # 決定資料庫路徑
    db_path = args.db_path
    if db_path is None:
        db_path = os.path.join(os.getcwd(), 'tasks', 'viva_results.db')
        
    if args.rebuild and os.path.exists(db_path):
        logger.info(f"選項 --rebuild 已啟用，刪除舊資料庫: {db_path}")
        os.remove(db_path)
    
    # 初始化資料庫連線
    db = db_manager.VIVADatabase(db_path=db_path)
    logger.info(f"使用資料庫路徑: {db.db_path}")
    
    success_count = 0
    total = len(args.task_dirs)
    
    # 逐一匯入指定的任務目錄
    for task_dir in args.task_dirs:
        if import_task(task_dir, db):
            success_count += 1
            
    logger.info(f"批次匯入完成: 成功匯入 {success_count}/{total} 個任務。")

if __name__ == '__main__':
    main()

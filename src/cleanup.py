import os
import shutil
import argparse
import json
from pathlib import Path

# 定義要清理的項目路徑（相對於任務根目錄）
CLEANUP_RULES = {
    "dirs": [
        "reads/original",
        "impurities_prefilter/bt2_alignment",
        "impurities_prefilter/bwa_alignment",
        "assembly",  # reference_prepare 的組裝暫存
    ],
    "spades_internal_dirs": [
        "K21", "K33", "K55", "tmp", "corrected"
    ],
    "file_patterns": [
        "alignment/bowtie2/*.bt2",
        "alignment/bwa/*.bwt",
        "alignment/bwa/*.sa",
        "alignment/bwa/*.ann",
        "alignment/bwa/*.amb",
        "alignment/bwa/*.pac",
        "alignment/bwa/*_unmapped_R*.fastq.gz",
    ]
}

def get_dir_size(path):
    total = 0
    try:
        for entry in os.scandir(path):
            if entry.is_file():
                total += entry.stat().st_size
            elif entry.is_dir():
                total += get_dir_size(entry.path)
    except (PermissionError, FileNotFoundError):
        pass
    return total

def cleanup_task(task_path, force=False):
    task_path = Path(task_path)
    if not (task_path / "log.txt").exists():
        return 0, 0  # 不是 VIVA 任務目錄

    print(f"\n[處理任務] {task_path.name}")
    total_freed = 0
    items_removed = 0

    # 1. 刪除指定目錄
    for d_rel in CLEANUP_RULES["dirs"]:
        d_path = task_path / d_rel
        if d_path.is_dir():
            size = get_dir_size(d_path)
            print(f"  - 移除目錄: {d_rel} ({size / 1024 / 1024:.2f} MB)")
            if force:
                shutil.rmtree(d_path)
            total_freed += size
            items_removed += 1

    # 2. 深入清理 SPAdes 目錄內部 (保留 contigs.fasta 和 log)
    unmapped_dir = task_path / "unmapped_analysis"
    if unmapped_dir.is_dir():
        for spades_dir in unmapped_dir.glob("*_spades_*"):
            if spades_dir.is_dir():
                for sub_d in CLEANUP_RULES["spades_internal_dirs"]:
                    target = spades_dir / sub_d
                    if target.is_dir():
                        size = get_dir_size(target)
                        print(f"  - 移除組裝暫存: {spades_dir.name}/{sub_d} ({size / 1024 / 1024:.2f} MB)")
                        if force:
                            shutil.rmtree(target)
                        total_freed += size
                        items_removed += 1

    # 3. 刪除特定模式檔案
    for pattern in CLEANUP_RULES["file_patterns"]:
        for f_path in task_path.glob(pattern):
            if f_path.is_file():
                size = f_path.stat().st_size
                print(f"  - 移除中間檔: {f_path.relative_to(task_path)} ({size / 1024 / 1024:.2f} MB)")
                if force:
                    f_path.unlink()
                total_freed += size
                items_removed += 1

    return total_freed, items_removed

def main():
    parser = argparse.ArgumentParser(description="VIVA 任務目錄清理工具")
    parser.add_argument("--dir", required=True, help="目標根目錄 (例如 tasks/ 或 NAS 掛載路徑)")
    parser.add_argument("--force", action="store_true", help="確認執行刪除 (不加此參數僅預覽)")
    args = parser.parse_args()

    target_root = Path(args.dir)
    if not target_root.is_dir():
        print(f"錯誤: 找不到目錄 {args.dir}")
        return

    if not args.force:
        print("!!! 注意: 目前為預覽模式 (Dry-run)，不會真正刪除檔案。加上 --force 執行刪除。 !!!")

    # 判斷是單一任務還是多任務目錄
    if (target_root / "log.txt").exists():
        tasks = [target_root]
    else:
        tasks = [d for d in target_root.iterdir() if d.is_dir()]

    grand_total_freed = 0
    total_items = 0
    
    for task in tasks:
        freed, items = cleanup_task(task, force=args.force)
        grand_total_freed += freed
        total_items += items

    print("\n" + "="*40)
    status = "實際釋放" if args.force else "預計可釋放"
    print(f"總計: {status} {grand_total_freed / 1024 / 1024:.2f} MB")
    print(f"清理項目數: {total_items}")
    if not args.force and total_items > 0:
        print("若確認無誤，請執行: python3 src/cleanup.py --dir <路徑> --force")

if __name__ == "__main__":
    main()

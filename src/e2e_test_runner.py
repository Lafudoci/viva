import argparse
import configparser
import logging
import os
import shutil
import sys
import time
from pathlib import Path
from mock_generator import MockDataGenerator
import new_task
from e2e_verifier import E2EVerifier

logger = logging.getLogger("e2e_test_runner")

def parse_preset(preset_path):
    config = configparser.ConfigParser()
    config.read(preset_path)
    if 'PRESET' not in config:
        raise ValueError(f"Invalid preset file: {preset_path}")
    
    return {
        "ref": config.get('PRESET', 'ref', fallback=None),
        "host": config.get('PRESET', 'remove_host', fallback=None),
        "impurities": config.get('PRESET', 'remove_impurities', fallback=None),
        "blastdb": config.get('PRESET', 'unmapped_blastdb', fallback=None)
    }

def get_host_path(host_name):
    # 這裡假設宿主序列在 /app/genomes/ 下，或者我們需要下載它
    # 為了測試，我們先找尋 genomes/ 目錄下的 .fna 或 .fasta
    if not host_name: return None
    
    # 常見的內建名稱轉換
    if host_name == "human":
        # 如果是 human，可能需要下載，這裡我們先拋出警告或找尋 dummy
        logger.warning("Scenario requires 'human' host. Ensure it exists in genomes/.")
        
    potential_path = Path("/app/genomes") / f"{host_name}.fasta"
    if not potential_path.exists():
        potential_path = Path("/app/genomes") / f"{host_name}.fna"
    
    return str(potential_path) if potential_path.exists() else None

def run_test(scenario_name, config, total_reads, output_dir, preset_path=None):
    logger.info(f"========== Starting Scenario: {scenario_name} ==========")
    
    # 1. 產生模擬讀序
    gen = MockDataGenerator(output_dir)
    r1, r2 = gen.generate_scenario(config, total_reads=total_reads)
    
    # 2. 執行 VIVA
    task_prefix = f"e2e_{scenario_name}_{int(time.time())}"
    viva_args = [
        "--single_task",
        "--prefix", task_prefix,
        "--ex_r1", str(r1),
        "--ex_r2", str(r2),
        "--threads", "8",
        "--auto_cleanup", "False" # 保留中間檔以便驗證
    ]
    
    if preset_path:
        viva_args.extend(["--preset", str(preset_path)])
    else:
        # 如果沒提供 preset，至少要給 --ref
        if config.get("target", {}).get("path"):
            viva_args.extend(["--ref", config["target"]["path"]])
    
    logger.info(f"Running VIVA for {task_prefix}...")
    task_id = new_task.main(viva_args)
    
    # 3. 驗證結果
    task_path = Path.cwd() / "tasks"
    expected_gt_path = output_dir / "ground_truth.json"
    
    verifier = E2EVerifier(task_id, task_path, expected_results_path=str(expected_gt_path))
    # 注意：這裡需要修改 e2e_verifier.py 來支援 ground_truth.json 的邏輯
    # 目前 e2e_verifier.py 是針對硬編碼的 expected_results.json
    results = verifier.verify()
    
    logger.info(f"Scenario {scenario_name} Result: {'PASS' if results['passed'] else 'FAIL'}")
    print(verifier.get_markdown_summary())
    
    return results['passed']

def run_lod_sweep(target_config, contaminant_config, host_config, total_reads, output_dir, preset_path=None):
    """執行 LOD Sweep 測試，從 10^-3 掃描到 10^-7"""
    logger.info("========== Starting LOD Sweep Test ==========")
    ratios = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
    results = []
    
    for ratio in ratios:
        s_name = f"lod_{ratio:.1e}"
        s_out = output_dir / s_name
        
        # 調整污染物比例，其餘給 host 或 target
        current_config = {
            "target": target_config.copy(),
            "contaminants": [
                {**contaminant_config, "ratio": ratio}
            ]
        }
        
        # 重新分配剩餘比例
        remaining = 1.0 - ratio
        if host_config:
            current_config["host"] = host_config.copy()
            current_config["host"]["ratio"] = remaining * 0.9
            current_config["target"]["ratio"] = remaining * 0.1
        else:
            current_config["target"]["ratio"] = remaining
            
        passed = run_test(s_name, current_config, total_reads, s_out, preset_path)
        results.append((ratio, passed))
        
    return results

def main():
    parser = argparse.ArgumentParser(description="VIVA E2E In Silico Test Runner")
    parser.add_argument("--preset", help="Path to .ini preset file")
    parser.add_argument("--ref", help="Manual path to reference FASTA (if no preset)")
    parser.add_argument("--host", help="Manual host name or path (if no preset)")
    parser.add_argument("--impurities", help="Manual path to impurities FASTA (if no preset)")
    parser.add_argument("--total_reads", type=int, default=100000, help="Total reads to simulate")
    parser.add_argument("--scenario", choices=["pure_target", "host_contamination", "targeted_analysis", "non_targeted_analysis", "lod", "all"], default="pure_target")
    parser.add_argument("--output_dir", default="./e2e_tests", help="Directory for mock data")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    
    preset_info = parse_preset(args.preset) if args.preset else {
        "ref": args.ref or str(Path(__file__).parent / "test_data" / "AC_000008.1.fasta"),
        "host": args.host,
        "impurities": args.impurities
    }
    
    # 如果預設的參考序列不存在且使用者沒給，報錯
    if not os.path.exists(preset_info["ref"]):
        logger.error(f"Reference file not found: {preset_info['ref']}. Please provide --ref.")
        sys.exit(1)

    # 定義場景設定
    scenarios = {}
    
    # 基礎配置
    target_cfg = {"path": preset_info["ref"], "ratio": 1.0} if preset_info.get("ref") else None
    host_path = get_host_path(preset_info.get("host"))
    host_cfg = {"path": host_path, "ratio": 0.9} if host_path else None
    
    # 場景 1: Pure Target
    if target_cfg:
        scenarios["pure_target"] = {"target": target_cfg}
        
    # 場景 2: Host Contamination (Targeted)
    if target_cfg and host_cfg:
        scenarios["host_contamination"] = {
            "target": {**target_cfg, "ratio": 0.1},
            "host": host_cfg
        }

    # 場景 3: Targeted Analysis (with Impurities)
    if target_cfg and preset_info.get("impurities"):
        scenarios["targeted_analysis"] = {
            "target": {**target_cfg, "ratio": 0.8},
            "host": {**host_cfg, "ratio": 0.1} if host_cfg else None,
            "contaminants": [
                {"path": preset_info["impurities"], "ratio": 0.1, "name": "Impurity"}
            ]
        }
    
    # 場景 4: Non-targeted Analysis
    if target_cfg:
        contaminant_path = preset_info.get("impurities") or preset_info.get("ref")
        scenarios["non_targeted_analysis"] = {
            "target": {**target_cfg, "ratio": 0.5},
            "contaminants": [
                {"path": contaminant_path, "ratio": 0.1, "name": "DiscoveryVirus"}
            ]
        }

    # 執行測試
    output_base = Path(args.output_dir)
    success = True
    
    if args.scenario == "lod":
        if not target_cfg:
            logger.error("LOD test requires a reference in preset.")
            sys.exit(1)
        contaminant_path = preset_info.get("impurities") or preset_info.get("ref")
        results = run_lod_sweep(target_cfg, {"path": contaminant_path, "name": "LOD_Virus"}, host_cfg, args.total_reads, output_base, args.preset)
        print("\nLOD Sweep Results:")
        for r, p in results:
            print(f"Ratio {r:.1e}: {'✅ PASS' if p else '❌ FAIL'}")
        return

    if args.scenario == "all":
        for s_name, s_config in scenarios.items():
            s_out = output_base / s_name
            if not run_test(s_name, s_config, args.total_reads, s_out, args.preset):
                success = False
    elif args.scenario in scenarios:
        s_config = scenarios[args.scenario]
        s_out = output_base / args.scenario
        success = run_test(args.scenario, s_config, args.total_reads, s_out, args.preset)
    else:
        logger.error(f"Scenario {args.scenario} not available with current preset.")
        sys.exit(1)
        
    if not success:
        logger.error("Some E2E tests FAILED.")
        sys.exit(1)
    else:
        logger.info("All E2E tests PASSED.")

if __name__ == "__main__":
    main()

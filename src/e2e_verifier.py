import json
import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)

class E2EVerifier:
    def __init__(self, task_id, task_path, expected_results_path=None):
        self.task_id = task_id
        self.task_path = Path(task_path)
        self.summary_path = self.task_path.joinpath(task_id, f"{task_id}_summary.json")
        self.report_md_path = self.task_path.joinpath(task_id, f"{task_id}_report.md")
        self.expected_results_path = expected_results_path
        self.results = {
            'completeness': [],
            'consistency': [],
            'passed': True
        }

    def verify(self):
        logger.info(f"Starting verification for task: {self.task_id}")
        
        # 1. 完整性檢查 (Completeness)
        self._check_completeness()
        
        # 2. 一致性檢查 (Consistency)
        if self.expected_results_path and os.path.exists(self.expected_results_path):
            self._check_consistency()
        else:
            logger.warning("Expected results file not found, skipping consistency check.")
            self.results['consistency'].append({"item": "Baseline Comparison", "status": "SKIPPED", "detail": "No expected_results.json"})

        self.results['passed'] = all(r['status'] == 'PASS' for r in self.results['completeness']) and \
                                 all(r['status'] in ('PASS', 'SKIPPED') for r in self.results['consistency'])
        
        return self.results

    def _check_completeness(self):
        # 檢查關鍵檔案是否存在
        critical_files = [
            (self.summary_path, "Summary JSON"),
            (self.report_md_path, "Markdown Report"),
            (self.task_path.joinpath(self.task_id, 'reads', 'fastp.json'), "fastp.json"),
            (self.task_path.joinpath(self.task_id, 'alignment', 'flagstat.json'), "flagstat.json")
        ]
        
        for file_path, label in critical_files:
            exists = file_path.exists()
            self.results['completeness'].append({
                "item": label,
                "status": "PASS" if exists else "FAIL",
                "detail": f"File exists: {file_path}" if exists else "File missing"
            })

        # 檢查 Summary 內容是否完整
        if self.summary_path.exists():
            try:
                with open(self.summary_path, 'r') as f:
                    summary = json.load(f)
                
                required_keys = ['fastp_abs', 'aln', 'cov', 'vc', 'version']
                for key in required_keys:
                    val = summary.get(key)
                    is_complete = val is not None and val != {}
                    self.results['completeness'].append({
                        "item": f"Summary Field: {key}",
                        "status": "PASS" if is_complete else "FAIL",
                        "detail": f"Field '{key}' has data" if is_complete else f"Field '{key}' is missing or empty"
                    })
            except Exception as e:
                self.results['completeness'].append({
                    "item": "Summary JSON Parse",
                    "status": "FAIL",
                    "detail": str(e)
                })

    def _check_consistency(self):
        try:
            with open(self.expected_results_path, 'r') as f:
                expected_data = json.load(f)
            
            with open(self.summary_path, 'r') as f:
                current_data = json.load(f)
            
            # 根據任務名稱前綴尋找對應的期望值 (例如 test_run)
            # 因為 task_id 包含時間戳，我們取前綴
            task_prefix = '_'.join(self.task_id.split('_')[:-1])
            
            # 判斷是否為 ground_truth.json (in silico 產出的)
            if "sources" in expected_data and "total_reads" in expected_data:
                self._check_ground_truth_consistency(expected_data, current_data)
                return

            expected = expected_data.get(task_prefix)
            
            if not expected:
                self.results['consistency'].append({
                    "item": f"Baseline for {task_prefix}",
                    "status": "SKIPPED",
                    "detail": "No baseline entry for this test type"
                })
                return

            # 比對 Mapping Rate (Bowtie2)
            if 'mapping_rate_bt2' in expected:
                # summary['aln']['mapped_rate']['bowtie2']['1'] -> "99.85 %"
                curr_rate_str = current_data.get('aln', {}).get('mapped_rate', {}).get('bowtie2', {}).get('1', '0 %')
                curr_rate = float(curr_rate_str.replace('%', '').strip())
                exp_rate = expected['mapping_rate_bt2']
                # 容許誤差 1%
                passed = abs(curr_rate - exp_rate) < 1.0
                self.results['consistency'].append({
                    "item": "Mapping Rate (Bowtie2)",
                    "status": "PASS" if passed else "FAIL",
                    "detail": f"Expected: {exp_rate}%, Found: {curr_rate}%"
                })

            # 比對 Variant Count
            if 'variant_count' in expected:
                curr_vc = 0
                vc_data = current_data.get('vc', {}).get('lofreq', {}).get('1', {})
                for pos in vc_data:
                    curr_vc += len(vc_data[pos].get('SNV', {}))
                
                exp_vc = expected['variant_count']
                passed = curr_vc == exp_vc
                self.results['consistency'].append({
                    "item": "Variant Count (LoFreq)",
                    "status": "PASS" if passed else "FAIL",
                    "detail": f"Expected: {exp_vc}, Found: {curr_vc}"
                })

        except Exception as e:
            self.results['consistency'].append({
                "item": "Consistency Check Execution",
                "status": "FAIL",
                "detail": str(e)
            })

    def _check_ground_truth_consistency(self, gt, current):
        """比對模擬產生的 Ground Truth 與實際執行結果"""
        # 1. Target Mapping Rate
        if "target" in gt["sources"]:
            exp_ratio = gt["sources"]["target"]["ratio"] * 100
            # summary['aln']['mapped_rate']['bowtie2']['1'] -> "99.85 %"
            curr_rate_str = current.get('aln', {}).get('mapped_rate', {}).get('bowtie2', {}).get('1', '0 %')
            curr_rate = float(curr_rate_str.replace('%', '').strip())
            
            # 容許誤差 5% (模擬過程中的隨機性以及去宿主/過濾的影響)
            passed = abs(curr_rate - exp_ratio) < 5.0
            self.results['consistency'].append({
                "item": "Target Mapping Rate (BT2)",
                "status": "PASS" if passed else "FAIL",
                "detail": f"Expected: ~{exp_ratio:.2f}%, Found: {curr_rate}%"
            })

        # 2. Host Removal Efficiency
        if "host" in gt["sources"]:
            host_reads = gt["sources"]["host"]["expected_reads"]
            removed_reads = current.get('remove_genome', {}).get('mapped_reads', 0)
            # 檢查是否至少移除了 90% 的模擬宿主讀序
            efficiency = (removed_reads / host_reads) * 100 if host_reads > 0 else 100
            passed = efficiency > 90.0
            self.results['consistency'].append({
                "item": "Host Removal Efficiency",
                "status": "PASS" if passed else "FAIL",
                "detail": f"Removed: {removed_reads}/{host_reads} ({efficiency:.2f}%)"
            })

        # 3. Non-targeted Discovery (BLAST)
        contaminants = [s for s in gt["sources"] if s.startswith("contaminant_")]
        if contaminants:
            unmapped_hits = current.get('unmapped_analysis', {})
            # 只要有找到任何一個模擬的污染物關鍵字在 BLAST 結果中
            found_count = 0
            for c_name in contaminants:
                clean_name = c_name.replace("contaminant_", "")
                for contig in unmapped_hits:
                    if clean_name.lower() in str(unmapped_hits[contig]).lower():
                        found_count += 1
                        break
            
            passed = found_count > 0
            self.results['consistency'].append({
                "item": "Non-targeted Discovery",
                "status": "PASS" if passed else "FAIL",
                "detail": f"Found {found_count}/{len(contaminants)} simulated contaminants"
            })

    def get_markdown_summary(self):
        lines = [f"### Test Results for {self.task_id}"]
        lines.append(f"**Overall Status: {'✅ PASS' if self.results['passed'] else '❌ FAIL'}**")
        
        lines.append("\n#### Completeness Check")
        lines.append("| Item | Status | Detail |")
        lines.append("| --- | --- | --- |")
        for r in self.results['completeness']:
            lines.append(f"| {r['item']} | {'✅ PASS' if r['status'] == 'PASS' else '❌ FAIL'} | {r['detail']} |")
            
        lines.append("\n#### Consistency Check")
        lines.append("| Item | Status | Detail |")
        lines.append("| --- | --- | --- |")
        for r in self.results['consistency']:
            status_icon = '✅ PASS' if r['status'] == 'PASS' else ('❌ FAIL' if r['status'] == 'FAIL' else '⚠️ SKIPPED')
            lines.append(f"| {r['item']} | {status_icon} | {r['detail']} |")
            
        return "\n".join(lines)

def run_verification(task_id, task_path, expected_results_path=None):
    verifier = E2EVerifier(task_id, task_path, expected_results_path)
    report = verifier.verify()
    md = verifier.get_markdown_summary()
    return report, md

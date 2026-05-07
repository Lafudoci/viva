import os
import subprocess
import json
import logging
import random
from pathlib import Path

logger = logging.getLogger(__name__)

class MockDataGenerator:
    def __init__(self, output_dir, seed=42):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.seed = seed
        self.ground_truth = {
            "total_reads": 0,
            "sources": {},
            "seed": seed
        }

    def _run_random_reads(self, ref_fasta, num_reads, output_prefix, read_len=150, insert_size=300):
        """呼叫 bbmap 的 randomreads.sh 產生模擬讀序"""
        r1_out = self.output_dir / f"{output_prefix}_R1.fastq"
        r2_out = self.output_dir / f"{output_prefix}_R2.fastq"
        
        cmd = [
            "randomreads.sh",
            f"ref={ref_fasta}",
            f"out={r1_out}",
            f"out2={r2_out}",
            f"reads={num_reads}",
            f"length={read_len}",
            "paired=t",
            f"seed={self.seed}",
            f"mininsert={insert_size-50}",
            f"maxinsert={insert_size+50}",
            "adderrors=t",
        ]
        
        logger.info(f"Generating {num_reads} reads from {ref_fasta}...")
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return r1_out, r2_out
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running randomreads.sh: {e.stderr}")
            raise

    def _generate_noise_ref(self, length=10000):
        """產生隨機序列作為 noise 來源"""
        noise_path = self.output_dir / "noise_ref.fasta"
        bases = ['A', 'C', 'G', 'T']
        seq = ''.join(random.choice(bases) for _ in range(length))
        with open(noise_path, "w") as f:
            f.write(">unmapped_noise\n")
            f.write(seq + "\n")
        return noise_path

    def generate_scenario(self, scenario_config, total_reads=100000):
        """
        根據場景設定產生混合讀序
        scenario_config: {
            "target": {"path": "...", "ratio": 0.4},
            "host": {"path": "...", "ratio": 0.5},
            "contaminants": [
                {"path": "...", "ratio": 0.05, "name": "Virus_A"},
                {"path": "...", "ratio": 0.05, "name": "Virus_B"}
            ]
        }
        """
        temp_files = []
        self.ground_truth["total_reads"] = total_reads
        self.ground_truth["sources"] = {}

        # 1. 建立來源列表
        sources_list = []
        if "target" in scenario_config:
            sources_list.append(("target", scenario_config["target"]))
        if "host" in scenario_config:
            sources_list.append(("host", scenario_config["host"]))
        if "contaminants" in scenario_config:
            for i, c in enumerate(scenario_config["contaminants"]):
                sources_list.append((f"contaminant_{c.get('name', i)}", c))

        # 檢查比例，不足 1.0 則補雜訊 (模擬 unmapped reads)
        total_ratio = sum(c["ratio"] for _, c in sources_list)
        noise_ref = None
        if total_ratio < 1.0:
            noise_ratio = 1.0 - total_ratio
            noise_ref = self._generate_noise_ref()
            sources_list.append(("unmapped_noise", {"path": str(noise_ref), "ratio": noise_ratio}))

        # 2. 計算各來源所需的讀序數並產生
        for name, config in sources_list:
            # BBMap randomreads.sh paired=t 時，reads=N 會產出 N 個 pairs (2N sequences)
            # 這裡的 total_reads 代表總序列數，故需除以 2
            reads_count = int(total_reads * config["ratio"] / 2)
            if reads_count == 0: continue
            
            ref_path = config["path"]
            r1, r2 = self._run_random_reads(ref_path, reads_count, name)
            temp_files.append((r1, r2))
            
            self.ground_truth["sources"][name] = {
                "ref": str(ref_path),
                "expected_reads": reads_count * 2,
                "ratio": config["ratio"]
            }

        # 3. 合併與洗牌 (Shuffle and Merge)
        final_r1 = self.output_dir / "mock_R1.fastq.gz"
        final_r2 = self.output_dir / "mock_R2.fastq.gz"
        
        self._merge_files([f[0] for f in temp_files], final_r1)
        self._merge_files([f[1] for f in temp_files], final_r2)

        # 4. 儲存 Ground Truth
        with open(self.output_dir / "ground_truth.json", "w") as f:
            json.dump(self.ground_truth, f, indent=4)

        # 5. 清理暫存檔
        for r1, r2 in temp_files:
            os.remove(r1)
            os.remove(r2)
        if noise_ref and noise_ref.exists():
            os.remove(noise_ref)

        return final_r1, final_r2

    def _merge_files(self, file_list, output_path):
        """合併並壓縮檔案"""
        logger.info(f"Merging files into {output_path}...")
        with open(output_path, "wb") as outfile:
            # 這裡可以使用 gzip 壓縮，但為了速度我們先產出 fastq
            # 或者直接呼叫 cat | gzip
            cat_cmd = ["cat"] + [str(f) for f in file_list]
            gzip_cmd = ["gzip", "-c"]
            
            p1 = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(gzip_cmd, stdin=p1.stdout, stdout=outfile)
            p1.stdout.close()
            p2.communicate()

if __name__ == "__main__":
    # 簡單的測試
    logging.basicConfig(level=logging.INFO)
    gen = MockDataGenerator("./test_mock_output")
    # 這裡需要實際存在的 FASTA 路徑才能執行
    # gen.generate_scenario(...)

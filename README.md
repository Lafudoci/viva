# VIVA 工具說明文件

VIVA (Virus Integrity and Variant Analyzer) 是一個自動化流程分析工具，專門用於病毒變異點鑑定與未知序列分析。本文件詳細介紹 VIVA 工具的運作流程原理、Docker 映像檔的建立與使用說明、快速使用範例、各項參數的詳細用法，以及批次分析功能的使用教學。

## 1. VIVA 運作流程原理 (Operation Flow / Principles)

VIVA 主要的分析流程如下（依執行先後順序）：

1. **Reads Preprocessing (讀序前處理與去宿主)** (`reads_preprocess.py`)
   - 匯入原始 R1/R2 Fastq 檔案，並計算 MD5 雜湊值。
   - 透過 [fastp](https://github.com/OpenGene/fastp) 進行序列品質過濾 (Quality Control) 與切除頭尾品質較差的序列 (`--global_trimming`)。
   - （選擇性）若提供去宿主選項 (`--remove_host`)，則利用 [Bowtie2](https://github.com/BenLangmead/bowtie2) 對指定宿主基因體進行讀序定位。隨後利用 `samtools` 萃取出**未定位上宿主的讀序** (unmapped reads)，有效移除宿主序列 (Host Sequences) 的干擾。
2. **Reference Preparation (參考序列準備)** (`reference_prepare.py`)
   - 將目標病原體參考序列 (Reference FASTA, `--ref`) 建立 `bwa` 與 `bowtie2` 索引 (index)，供後續讀序定位使用。
3. **Impurities Prefiltering (不純物序列過濾)** (`impurities_prefilter.py`)
   - （選擇性）若有提供不純物序列檔 (`--remove_impurities`)，則將經過前處理的讀序對此不純物資料庫進行讀序定位。過濾能定位上不純物序列的讀序，僅保留未命中序列以進入主分析。
4. **Reads Alignment (讀序定位)** (`reads_alignment.py`)
   - 將乾淨的目標讀序使用設定的讀序定位軟體 (`--alns`，預設為 `bowtie2,bwa`) 定位至參考序列。
   - 透過 `samtools` 進行格式轉換、排序 (sorting)、產生索引檔 (indexing)，並紀錄精確的覆蓋度 (Coverage) 與深度 (Depth) 等統計資料。
5. **Unmapped Analysis (未定位序列分析)** (`unmapped_analysis.py`)
   - （選擇性）針對無法定位上目標病原體參考序列的讀序 (unmapped reads)，使用 SPAdes (metaSPAdes) 進行 de novo 序列組裝 (`--unmapped_assemble`)，產生可能為其他微生物的重疊群 (contigs)。
   - 利用 [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) 將產生的組裝序列或原始讀序與外部資料庫 (`--unmapped_blastdb` 或 RVDB 等) 進行比對。此步驟搭配自選過濾條件 (長度 `--unmapped_len_filter`, 同源性 `--unmapped_ident_filter`) 找尋潛在的未知病原體。
6. **Variant Calling (變異點鑑定)** (`variant_calling.py`)
   - 主比對步驟完成後，進一步偵測不同病毒株或樣本間的單一核苷酸變異 (SNVs) 及小段序列的插入/缺失 (InDels)。
   - 此步驟使用 `lofreq` 與 `varscan`，並結合給定的篩選閾值（發生頻率 `--vc_threshold` 及品質分數 `--min_vc_score`）去除雜訊並排除偽陽性 (False positive) 點位，列出高可信度的變異位點。
7. **Report & Summary Generation (報告生成)** (`summary_generator.py`, `report_generator.py`)
   - 統整分析過程中所有的產出：包含品管統計結果、覆蓋深度折線圖、未定位序列的 BLAST 物種標定與變異點列表。
   - 最終產出便於查閱的綜合 HTML 報告及 CSV 分析總結。

---

## 2. Docker Image Building 說明

如需自行編譯 VIVA 工具專屬的 Docker Image（例如映像檔名稱為 `viva:v1.11.1`），請在**專案原始碼目錄**下執行以下指令：

```bash
# 1. 開啟終端機並切換到 VIVA 專案原始碼根目錄
cd /home/leftc/Google\ Drive/GitHub-GD/viva  # 視您的實際存放位置而定

# 2. 利用原始碼目錄中的 Dockerfile 建置 Docker Image
#    參數 -t 用以標記產生的 Image 名稱及版號
sudo docker build -t viva:v1.11.1 .
```

---

## 3. 快速使用範例 (Quick Start Example)

本範例為一個**單一檔案分析任務** (`--single_task`) 的執行方式。請在您的**工作目錄**下準備好相應資料後執行：

```bash
# 切換到您的分析工作目錄 (包含任務輸出與相關參考檔案的目錄)
cd /home/ra5/viva

# 以 Docker Container 啟動分析任務
sudo docker run -i --rm \
  -v $(pwd)/tasks:/app/tasks \
  -v $(pwd)/genomes:/app/genomes \
  -v $(pwd)/blastdb:/app/blastdb \
  -v /home:/home \
  viva:v1.11.1 \
  --single_task \
  --prefix TFDA-MPXV-20250918 \
  --ref "/home/ra5/ref-fasta/MPXV/OR030941.1.fasta" \
  --remove_host GCA_023783515.1_ASM2378351v1.1 \
  --blastdb_path /home/ra5/bioapp/blastdb \
  --rvdb_anno_path /home/ra5/bioapp/blastdb/RVDBv30_AnnotationList_Jun2025.tab \
  --unmapped_blastdb "U-RVDBv30.0.fasta" \
  --unmapped_blastdb_extra_list "core_nt nt_prok nt_others 16S_ribosomal_RNA" \
  --unmapped_len_filter 100 \
  --min_vc_score 4 \
  --threads 20 \
  --ex_r1 /home/ra5/NGS-reads/20250918-RSV-MPXV/MPV_S3_L001_R1_001.fastq.gz \
  --ex_r2 /home/ra5/NGS-reads/20250918-RSV-MPXV/MPV_S3_L001_R2_001.fastq.gz \
  --task_note "署內培養 MPXV 分析 (去除宿主基因體序列)"
```

---

## 4. 詳細參數使用教學及舉例

主要參數整理自系統底層 (`tasks_manager.py` 及 `new_task.py` 的 argparse 設定)：

### I. 任務執行模式與輸入/輸出參數
| 參數 | 說明 | 舉例 |
| --- | --- | --- |
| `--single_task` | (必要) 以單個樣本分析模式執行程式。 | `--single_task` |
| `--task_sheet` | (批次模式) 用於批次任務，提供包含設定檔格式 `.ini` 之檔案。 | `--task_sheet batch_job.ini` |
| `--prefix` | 此單次任務的總命名前綴，決定最終分析資料夾或檔案名稱。 | `--prefix TestSample001` |
| `--ex_r1` | 外部 R1 讀序的檔案絕對路徑 (支援 `fastq.gz` 格式)。 | `--ex_r1 /data/sample_R1.fastq.gz` |
| `--ex_r2` | 外部 R2 讀序的檔案絕對路徑。 | `--ex_r2 /data/sample_R2.fastq.gz` |

### II. 目標讀序定位參數 (Reference & Filtration)
| 參數 | 說明 | 舉例 |
| --- | --- | --- |
| `--ref` | 分析欲比對之參考序列 (Reference FASTA file) 路徑。若無提供則自動切換為 De novo 分析模式。 | `--ref /home/ref/RSV.fasta` |
| `--remove_host` | 指定並移除特定的宿主序列以提升分析效能。可使用內建字詞 (`human`, `dog`, `vero`, `chicken`, `rhesus_monkey`)，或提供位於 `/app/genomes/` 下的自訂基因體名稱。 | `--remove_host human`<br>`--remove_host GCA_023783515.1` |
| `--remove_impurities` | 進階雜訊去除功能；提供不純物序列之參考 FASTA 檔，用來過濾對應的讀序。 | `--remove_impurities /data/noise.fasta` |
| `--alns` | 若有提供參考序列，可選擇讀序定位軟體，多重選擇時以逗號隔開 (預設為 `bowtie2,bwa`)。 | `--alns bwa` |

### III. 統轄運算與系統資源分配
| 參數 | 說明 | 舉例 |
| --- | --- | --- |
| `--threads` | 限定此任務分配能使用的最大 CPU 執行緒數目 (預設為 `6`)。 | `--threads 20` |
| `--spades_mem` | 限制 `spades` de novo Assemble 與 unmapped Assemble 所佔用的最大記憶體容量 (以 GB 為單位，預設為 `22`)，以避免主機資源耗竭。 | `--spades_mem 32` |

### IV. 組裝與 BLAST 未定位序列參數
| 參數 | 說明 | 舉例 |
| --- | --- | --- |
| `--unmapped_assemble` | 設定是否針對未定位序列利用 metaSPAdes 進行組裝 (預設為 `True`)。 | `--unmapped_assemble True` |
| `--blastdb_path` | 主機上自建之 BLAST 資料庫存放目錄。 | `--blastdb_path /home/bioapp/blastdb` |
| `--unmapped_blastdb` | 定義針對未定位的組裝 sequences 或 reads 優先進行比對的專用 BLAST 資料庫，檔名須存在於 `--blastdb_path` 目錄下。 | `--unmapped_blastdb "U-RVDBv30.0.fasta"` |
| `--unmapped_blastdb_extra_list`| 若需比對多個資料庫，可提供額外之資料庫名稱字串(各名稱以空格隔開)，系統將循序進行比對查詢。 | `--unmapped_blastdb_extra_list "core_nt nt_prok"` |
| `--unmapped_len_filter` | BLAST 比對後，過濾掉長度低於此閾值的序列結果 (預設 `500` bp)。 | `--unmapped_len_filter 100` |
| `--unmapped_ident_filter`| BLAST 比對後，過濾掉同源性 (Identity) 低於此百分比的序列結果 (預設 `95`%)。 | `--unmapped_ident_filter 90` |
| `--rvdb_anno_path` | 使用 RVDB 資料庫時，給定實體註解的 `.tab` 檔案來輔助擷取完整的生物分類註解並呈現在報告中。 | `--rvdb_anno_path /home/.../RVDBv30.tab` |

### V. 變異點分析擷取條件 (Variant Calling)
| 參數 | 說明 | 舉例 |
| --- | --- | --- |
| `--min_vc_score` | 最低變異點分數品質閥門檻，用來剔除低品質的偽陽性 (False positive) 點位 (預設值 `1`)。 | `--min_vc_score 4` |
| `--vc_threshold` | 用來過濾變異位點上發生頻率 (Allele Frequency, AF) 過低的等位基因 (預設 `0.7` 即為 70%)。 | `--vc_threshold 0.5` |

### VI. 品管控制及其他補充參數
| 參數 | 說明 | 舉例 |
| --- | --- | --- |
| `--global_trimming` | 設定 fastp，將每條 R1 與 R2 都在 5' 與 3' 端一併切除此數字設定的鹼基長度 (預設為 `0`)。 | `--global_trimming 5` |
| `--preset_path` | 讀取已記錄分析流程常態配置的 `.ini` 設定檔，簡化指令輸入 (預設為 `None`)。 | `--preset_path /path/profile.ini` |
| `--task_note` | 任務之註解文字說明，將會呈現在最終報告首頁。(若文字中包含空格請以雙引號包夾)。 | `--task_note "署內培養檢體"` |
| `--sample_product_name` | 樣本名稱詮釋資料 (Metadata)，僅供報告註解用途。 |- |
| `--sample_product_lot` | 樣本批號註記，僅供報告註解用途。 |- |
| `--sample_sequencing_date` | 樣本定序日期註記，僅供報告註解用途。 |- |
| `--sample_note` | 樣本專屬注意事項備註，僅供報告註解用途。 |- |

---

## 5. 批次分析功能 (Batch Analysis) 使用教學

當您需要一次處理多組樣本時，可以透過指定一份 `.ini` 格式的 Batch Task Sheet 來自動化依序執行多個 VIVA 任務。

### 初始化批次任務設定檔 (Task Sheet)

首先，建立一個副檔名為 `.ini` 的批次任務設定檔（例如 `batch_tasks.ini`）。該檔案支援定義單一或多個群組任務：

```ini
[Task_Group_1]
# 第一組批次任務設定
# 指定包含欲分析之樣本 ID 清單的文字檔路徑 (每行一個 ID)
samples_list_filepath = /home/ra5/viva/batches/group1_samples.txt

# 指定記載樣本路徑與 metadata 的總表 .ini 路徑
read_meta_path = /home/ra5/viva/batches/all_reads_meta.ini

# 指定此批次任務統一掛載的分析設定檔 (Preset Profiler)
preset_path = /home/ra5/viva/presets/mpxv_default.ini

# 批次任務註記，會一併反應在報告上的註解區
batch_tasks_note = "Group 1 MPXV Routine Analysis"


[Task_Group_2]
# 第二組批次任務設定
samples_list_filepath = /home/ra5/viva/batches/group2_samples.txt
read_meta_path = /home/ra5/viva/batches/all_reads_meta.ini
preset_path = /home/ra5/viva/presets/rsv_default.ini
batch_tasks_note = "Group 2 RSV Routine Analysis"
```

### 準備 Metadata 總表 (`all_reads_meta.ini`)

在 `read_meta_path` 指定的參考檔案中，需要有以每支樣本 ID 為段落的設定資訊：

```ini
[Sample_001]
product = MPXV
lot = 20250918-01
seq_date = 2025-09-18
reads_note = NA
ex_r1 = /home/ra5/NGS-reads/20250918-RSV-MPXV/Sample_001_R1_001.fastq.gz
ex_r2 = /home/ra5/NGS-reads/20250918-RSV-MPXV/Sample_001_R2_001.fastq.gz
```

### 執行批次分析指令

當設定檔皆部署完畢後，在**工作目錄**下同樣地透過 Docker 執行，但使用 `--task_sheet` 取代 `--single_task` 參數：

```bash
cd /home/ra5/viva

sudo docker run -i --rm \
  -v $(pwd)/tasks:/app/tasks \
  -v $(pwd)/genomes:/app/genomes \
  -v $(pwd)/blastdb:/app/blastdb \
  -v /home:/home \
  viva:v1.11.1 \
  --task_sheet /home/ra5/viva/batches/batch_tasks.ini
```

批次運行時，系統將會建立一個具有 `batch_task_YYYYMMDDHHMM` 命名的序列清單狀態檔案記錄任務排程，並循序依每個樣本自動載入 `preset_path` 的參數，完成所有樣本的比對工作以及最終的匯總 CSV 報告產生 (`batch_task_report.py`)。

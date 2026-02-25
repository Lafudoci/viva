import sqlite3
import os
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

class VIVADatabase:
    def __init__(self, db_path=None):
        """
        初始化 VIVADatabase，預設路徑建立在 tasks 掛載目錄下，
        以便跨任務持久化儲存。
        """
        if db_path is None:
            self.db_path = os.path.join(os.getcwd(), 'tasks', 'viva_results.db')
        else:
            self.db_path = db_path
            
        # Ensure the tasks directory exists
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        self._init_tables()

    def _get_connection(self):
        return sqlite3.connect(self.db_path)

    def _init_tables(self):
        """建立 SQLite 關聯資料表 (如果尚不存在)"""
        conn = self._get_connection()
        cursor = conn.cursor()

        # 1. 任務主表 (Tasks)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS Tasks (
                task_id TEXT PRIMARY KEY,
                task_name TEXT,
                start_date TEXT,
                finish_date TEXT,
                status TEXT, -- 'Running', 'Failed', 'Completed'
                error_log TEXT,
                preset_id TEXT,
                task_note TEXT,
                product_name TEXT,
                product_lot TEXT,
                seq_date TEXT,
                viva_version TEXT,
                preset_version TEXT
            )
        ''')

        # 2. 品管與除宿主數據表 (QC_Metrics)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS QC_Metrics (
                task_id TEXT PRIMARY KEY,
                reads_r1_name TEXT,
                reads_r2_name TEXT,
                fastp_before_reads INTEGER,
                fastp_after_reads INTEGER,
                fastp_after_q30 REAL,
                duplication_rate REAL,
                mean_length_r1 INTEGER,
                mean_length_r2 INTEGER,
                remove_host_name TEXT,
                remove_host_reads INTEGER,
                remove_host_pct REAL,
                FOREIGN KEY(task_id) REFERENCES Tasks(task_id)
            )
        ''')

        # 3. 比對結果數據表 (Alignment_Stats)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS Alignment_Stats (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                task_id TEXT,
                ref_order INTEGER,
                ref_header TEXT,
                aligner TEXT,
                mapped_reads INTEGER,
                mapped_rate REAL,
                coverage REAL,
                mean_depth REAL,
                FOREIGN KEY(task_id) REFERENCES Tasks(task_id)
            )
        ''')

        # 4. 變異點表 (Variants)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS Variants (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                task_id TEXT,
                ref_order INTEGER,
                position INTEGER,
                ref_base TEXT,
                alt_base TEXT,
                caller TEXT, -- 'lofreq', 'varscan'
                bt2_af TEXT,
                bwa_af TEXT,
                FOREIGN KEY(task_id) REFERENCES Tasks(task_id)
            )
        ''')

        # 5. 雜質過濾數據表 (Impurities_Stats)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS Impurities_Stats (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                task_id TEXT,
                impurity_name TEXT,
                aligner TEXT,
                mapped_reads INTEGER,
                coverage REAL,
                FOREIGN KEY(task_id) REFERENCES Tasks(task_id)
            )
        ''')

        # 6. BLAST 結果數據表 (BLAST_Results)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS BLAST_Results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                task_id TEXT,
                db_name TEXT,
                hit_count INTEGER,
                FOREIGN KEY(task_id) REFERENCES Tasks(task_id)
            )
        ''')

        conn.commit()
        conn.close()

    def create_task(self, task_id, task_name, start_date, preset_id=None, task_note=None, product=None, lot=None, seq_date=None):
        """任務初始化時寫入 Running 狀態"""
        try:
            conn = self._get_connection()
            cursor = conn.cursor()
            cursor.execute('''
                INSERT OR REPLACE INTO Tasks 
                (task_id, task_name, start_date, status, preset_id, task_note, product_name, product_lot, seq_date)
                VALUES (?, ?, ?, 'Running', ?, ?, ?, ?, ?)
            ''', (task_id, task_name, start_date, preset_id, task_note, product, lot, seq_date))
            conn.commit()
        except sqlite3.Error as e:
            logger.error(f"DB Error (create_task): {e}")
        finally:
            if conn: conn.close()

    def update_task_status(self, task_id, status, error_log=None):
        """於流程中更新任務狀態，例如 Failed 或記錄 Error"""
        try:
            conn = self._get_connection()
            cursor = conn.cursor()
            
            # 若為 Completd，記錄結束時間
            finish_date = None
            if status in ['Completed', 'Failed']:
               finish_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
               cursor.execute('''
                    UPDATE Tasks 
                    SET status=?, error_log=?, finish_date=?
                    WHERE task_id=?
                ''', (status, str(error_log) if error_log else None, finish_date, task_id))
            else:
                 cursor.execute('''
                    UPDATE Tasks 
                    SET status=?, error_log=?
                    WHERE task_id=?
                ''', (status, str(error_log) if error_log else None, task_id))
            conn.commit()
        except sqlite3.Error as e:
            logger.error(f"DB Error (update_task_status): {e}")
        finally:
            if conn: conn.close()
            
    def save_final_summary(self, task_id, s_dict):
        """讀取最終的 summary.json dict 並寫入資料庫表"""
        try:
            conn = self._get_connection()
            cursor = conn.cursor()
            
            # 1. Update task general info that might be missing
            viva_version = s_dict.get('version', {}).get('viva')
            preset_version = s_dict.get('preset', {}).get('version')
            
            cursor.execute('''
                UPDATE Tasks 
                SET finish_date=?, viva_version=?, preset_version=?
                WHERE task_id=?
            ''', (s_dict.get('finish_date'), viva_version, preset_version, task_id))

            # 2. Insert QC Metrics
            r_meta = s_dict.get('reads_meta', {}).get('reads_file_meta', {}).get('file_name', {})
            f_abs = s_dict.get('fastp_abs', {})
            rm_host = s_dict.get('remove_genome', {})
            
            cursor.execute('''
                INSERT OR REPLACE INTO QC_Metrics
                (task_id, reads_r1_name, reads_r2_name, fastp_before_reads, fastp_after_reads, 
                fastp_after_q30, duplication_rate, mean_length_r1, mean_length_r2,
                remove_host_name, remove_host_reads, remove_host_pct)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                task_id,
                r_meta.get('r1'), r_meta.get('r2'),
                f_abs.get('before_total_reads'), f_abs.get('after_total_reads'),
                f_abs.get('after_total_q30'), f_abs.get('duplication_rate'),
                f_abs.get('after_r1_length'), f_abs.get('after_r2_length'),
                rm_host.get('genome'), rm_host.get('mapped_reads'), rm_host.get('remove_percentage')
            ))
            
            # 3. Insert Alignment Stats
            ref_num = int(s_dict.get('ref_meta_dict', {}).get('ref_num', 0))
            for order in range(1, ref_num + 1):
                o_str = str(order)
                header = s_dict.get('ref_meta_dict', {}).get('seq_meta', {}).get(o_str, {}).get('fasta_header_escape')
                
                for aligner in ['bowtie2', 'bwa']:
                    m_rate = s_dict.get('aln', {}).get('mapped_rate', {}).get(aligner, {}).get(o_str)
                    m_reads = s_dict.get('aln', {}).get('mapped_reads', {}).get(aligner, {}).get(o_str)
                    cov_pct = s_dict.get('cov', {}).get(aligner, {}).get(o_str, {}).get('coverage')
                    mean_dp = s_dict.get('cov', {}).get(aligner, {}).get(o_str, {}).get('meandepth')
                    
                    if m_reads is not None:
                        cursor.execute('''
                            INSERT INTO Alignment_Stats
                            (task_id, ref_order, ref_header, aligner, mapped_reads, mapped_rate, coverage, mean_depth)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (task_id, order, header, aligner, m_reads, m_rate, cov_pct, mean_dp))
                        
            # 4. Insert Variants (LoFreq)
            for order in range(1, ref_num + 1):
                o_str = str(order)
                lofreq_snv = s_dict.get('vc', {}).get('lofreq', {}).get(o_str, {})
                for pos, v_data in lofreq_snv.items():
                    ref = v_data.get('REF')
                    for alt, aln_result in v_data.get('SNV', {}).items():
                        bt2_af = aln_result.get('bowtie2', {}).get('FREQ')
                        bwa_af = aln_result.get('bwa', {}).get('FREQ')
                        cursor.execute('''
                            INSERT INTO Variants
                            (task_id, ref_order, position, ref_base, alt_base, caller, bt2_af, bwa_af)
                            VALUES (?, ?, ?, ?, ?, 'lofreq', ?, ?)
                        ''', (task_id, order, pos, ref, alt, bt2_af, bwa_af))
                        
                # Varscan
                varscan_snv = s_dict.get('vc', {}).get('varscan', {}).get(o_str, {})
                for pos, v_data in varscan_snv.items():
                    ref = v_data.get('REF')
                    for alt, aln_result in v_data.get('SNV', {}).items():
                        bt2_af = aln_result.get('bowtie2', {}).get('FREQ')
                        bwa_af = aln_result.get('bwa', {}).get('FREQ')
                        cursor.execute('''
                            INSERT INTO Variants
                            (task_id, ref_order, position, ref_base, alt_base, caller, bt2_af, bwa_af)
                            VALUES (?, ?, ?, ?, ?, 'varscan', ?, ?)
                        ''', (task_id, order, pos, ref, alt, bt2_af, bwa_af))
                        
            # 5. Insert Impurities Stats
            im_meta = s_dict.get('impurit_filter_meta', {}).get('seq_meta', {})
            im_results = s_dict.get('impurit_filter_results', {})
            im_cov = s_dict.get('impurit_filter_coverage', {})
            
            for o_str, meta in im_meta.items():
                imp_name = meta.get('fasta_header_escape') or meta.get('fasta_header')
                for aligner in ['bt2', 'bwa']:
                    # some times the aligner keys are 'bowtie2' instead of 'bt2' in some places, 
                    # but looking at summary.json it uses 'bt2' and 'bwa' for impurities.
                    m_reads = im_results.get(o_str, {}).get(aligner, {}).get('mapped_reads')
                    cov_pct = im_cov.get(o_str, {}).get(aligner, {}).get('coverage')
                    
                    if m_reads is not None and m_reads != "N/A":
                        # Convert tool 'bt2' to 'bowtie2' to match alignment stats, or leave as 'bt2'
                        aligner_name = 'bowtie2' if aligner == 'bt2' else 'bwa'
                        cursor.execute('''
                            INSERT INTO Impurities_Stats
                            (task_id, impurity_name, aligner, mapped_reads, coverage)
                            VALUES (?, ?, ?, ?, ?)
                        ''', (task_id, imp_name, aligner_name, int(m_reads), float(cov_pct) if cov_pct else None))
                        
            # 6. Insert BLAST Results
            unmapped = s_dict.get('unmapped_analysis', {})
            for contig_name, blast_info in unmapped.items():
                if contig_name == "N/A":
                    continue
                db_name = blast_info.get('BLASTdb_name')
                if not db_name or db_name == "N/A":
                    continue
                
                # Check if this db_name is already in our temp dict to count hits
                # In summary.json, unmapped_analysis is keyed by contig name (e.g. "NODE_1_...").
                # Let's count hits per db_name.
                pass
            
            # Since unmapped_analysis has { "NODE_1...": {"BLASTdb_name": "nt", ...}, ... }
            # Let's aggregate hit count per db_name for this sample
            blast_counts = {}
            unmapped = s_dict.get('unmapped_analysis', {})
            for contig_name, blast_info in unmapped.items():
                if contig_name == "N/A":
                    continue
                db_name = blast_info.get('BLASTdb_name')
                if not db_name or db_name == "N/A":
                    continue
                blast_counts[db_name] = blast_counts.get(db_name, 0) + 1
                
            for db_name, hit_count in blast_counts.items():
                cursor.execute('''
                    INSERT INTO BLAST_Results
                    (task_id, db_name, hit_count)
                    VALUES (?, ?, ?)
                ''', (task_id, db_name, hit_count))

            conn.commit()
        except sqlite3.Error as e:
            logger.error(f"DB Error (save_final_summary): {e}")
        finally:
             if conn: conn.close()

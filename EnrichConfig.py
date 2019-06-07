import kobas.config as kobas_config
import os


PEP_URL_TEMP = 'ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD/seq_pep/{}.pep.fasta.gz'
DB_URL_TEMP = 'ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD/sqlite3/{}.db.gz'
BLAST_DEFAULT = '-evalue 1e-5 -outfmt 6 -max_target_seqs 1'

KOBASRC = kobas_config.getrc()
BLAST_DIR = KOBASRC['blastout']
PEP_DIR = KOBASRC['blastdb']
DB_DIR = KOBASRC['kobasdb']
KOBAS_RUN = KOBASRC['kobas_run']
KOBAS_PY = KOBASRC['kobas_py'] 
WHEAT_PEP_DB = os.path.join(BLAST_DIR, 'wheat.pep.fa')
TERM_CAT = os.path.join(BLAST_DIR, 'term.cat.txt')

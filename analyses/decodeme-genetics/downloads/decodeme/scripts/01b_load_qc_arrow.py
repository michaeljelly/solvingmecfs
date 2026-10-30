# Leaner version of 01: arrow-based reading/filtering of remaining GWAS files to QC-passed variants
import pyarrow as pa, pyarrow.csv as pc, pyarrow.compute as C, pyarrow.parquet as pq, sys
W='/tmp/work_decodeme/'
qc=pc.read_csv(W+'ss/gwas_qced.var.gz').column('ID').combine_chunks()
cols=['CHROM','GENPOS','ID','ALLELE0','ALLELE1','A1FREQ','A1FREQ_CASES','N','N_CASES','N_CONTROLS','BETA','SE','LOG10P']
for g in sys.argv[1:]:
    t=pc.read_csv(W+f'ss/{g}.regenie.gz',parse_options=pc.ParseOptions(delimiter=' '),
        convert_options=pc.ConvertOptions(include_columns=cols,column_types={'CHROM':pa.int64()}))
    t=t.filter(C.is_in(t['ID'],value_set=qc))
    n=t.num_rows; sig=C.sum(C.greater(t['LOG10P'],7.30103)).as_py()
    print(g,'rows QC',n,'Ncases',C.max(t['N_CASES']).as_py(),'Nctrl',C.max(t['N_CONTROLS']).as_py(),'GW sig',sig,flush=True)
    pq.write_table(t,W+f'pq/{g}.parquet'); del t

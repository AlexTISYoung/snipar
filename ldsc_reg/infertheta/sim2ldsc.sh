ldsc_repo="/disk/genetics2/pub/repo/ssgac/ldsc_mod"

python3 ldsc_reg/sim2ldsc_ldscores.py
python3 ${ldsc_repo}/munge_sumstats.py --sumstats ldsc_reg/ldscores/simdata2ldsc_dir.csv --out ldsc_reg/ldscores/simdata.sumstats
python3 ${ldsc_repo}/ldsc.py --ref-ld-chr ldsc_reg/ldscores/ --out ldsc_reg/ldscores/h2_est --h2 ldsc_reg/ldscores/simdata.sumstats.sumstats.gz --w-ld-chr ldsc_reg/ldscores/


ldsc_repo="/disk/genetics2/pub/repo/ssgac/ldsc_mod" 
eur_w_ld_chr="/var/genetics/pub/data/ld_ref_panel/eur_w_ld_chr/"
merge_alleles="/disk/genetics2/pub/data/PH3_Reference/w_hm3.snplist"

# Generating CSV's
python3 ldsc_reg/infertheta/gen_output.py

# Munging generated CSV's
python3 ${ldsc_repo}/munge_sumstats.py --sumstats ldsc_reg/infertheta/ldsc_reg/simulated_data_dir.csv --out ldsc_reg/infertheta/ldsc_reg/simulated_data_dir

# Run LDSC Reg
python3 ${ldsc_repo}/ldsc.py --ref-ld-chr ${eur_w_ld_chr} --out ldsc_reg/infertheta/ldsc_reg/rg_delta_beta --h2 ldsc_reg/infertheta/ldsc_reg/simulated_data_dir.sumstats.gz --w-ld-chr ${eur_w_ld_chr}
 

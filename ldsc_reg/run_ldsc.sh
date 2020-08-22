repo="/disk/genetics2/pub/repo/ssgac/ldsc_mod" 
eur_w_ld_chr="/var/genetics/pub/data/ld_ref_panel/eur_w_ld_chr/"

# Generating CSV's
python3 ldsc_reg/gen_output.py

# Munging generated CSV's
python3 ${ldsc_repo}/munge_sumstats.py --sumstats ldsc_reg/simulated_data_dir.csv --out ldsc_reg/simulated_data_dir

python3 ${ldsc_repo}/munge_sumstats.py --sumstats ldsc_reg/simulated_data_par.csv --out ldsc_reg/simulated_data_par

# Run LDSC Reg
python3 ${ldsc_repo}/ldsc.py --ref-ld-chr ${eur_w_ld_chr} --out ldsc_reg/rg_delta_beta --rg ldsc_reg/simulated_data_dir.sumstats.gz,ldsc_reg/simulated_data_par.sumstats.gz --w-ld-chr ${eur_w_ld_chr}
 


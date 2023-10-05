rm outputs/tmp/*
rm snipar/tests/test_data/*
fams=3000
nsnps=1000
for i in {1..2}
do
python snipar/simulate/simulate_pop.py ${nsnps} 0.5 ${fams} $((fams/3)) $((fams/3)) 0.5 "outputs/tmp/t__t${i}" --chrom ${i}
python -c "import pandas as pd
from snipar.pedigree import create_pedigree
print('********************')
print('outputs/tmp/t__t${i}_fams.ped')
print('********************')
ped = pd.read_csv('outputs/tmp/t__t${i}_fams.ped', delim_whitespace=True);
ped.to_csv('snipar/tests/test_data/sample_parent_sib${i}.ped', index = False, sep = ' ');
offsprings = ped[~ (ped['IID'].str.endswith('_M') | ped['IID'].str.endswith('_P'))]
offsprings.to_csv('snipar/tests/test_data/sample_sib${i}.ped', index = False, sep = ' ');
remove = pd.read_csv('outputs/tmp/t__t1_remove.txt', delim_whitespace=True, names = ['counter', 'id']);
ped = ped[~ped['IID'].isin(remove['id'])];
ped.to_csv('snipar/tests/test_data/sample${i}.ped', index = False, sep = ' ');
sibships = ped.groupby(['FID','FATHER_ID', 'MOTHER_ID']).agg({'IID':list}).reset_index();
ids = set(ped['IID'].tolist());
kinship = [];
for index, row in sibships.iterrows():
    for i in range(1,len(row['IID'])):
        for j in range(i):
            kinship.append([row['FID'], row['IID'][i], row['FID'], row['IID'][j], 'FS']);
for index, row in ped.iterrows():
    if row['FATHER_ID'] in ids:
        kinship.append([row['FID'], row['IID'], row['FID'], row['FATHER_ID'],'PO']);
    if row['MOTHER_ID'] in ids:
        kinship.append([row['FID'], row['IID'], row['FID'], row['MOTHER_ID'],'PO']);
kinship_df = pd.DataFrame(kinship, columns = ['FID1', 'ID1', 'FID2', 'ID2', 'InfType']);
kinship_df.to_csv('snipar/tests/test_data/sample${i}.king', sep = '\t', index = False);
"
cp outputs/tmp/t__t${i}.agesex snipar/tests/test_data/sample${i}.agesex
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g outputs/tmp/t__t${i}.haps -s outputs/tmp/t__t${i}.sample -og snipar/tests/test_data/sample${i}.bgen -os snipar/tests/test_data/sample${i}.sample -filetype shapeit_haplotypes -ofiletype bgen
plink2 --bgen snipar/tests/test_data/sample${i}.bgen ref-last --sample snipar/tests/test_data/sample${i}.sample --make-bed --out snipar/tests/test_data/sample${i} --oxford-single-chr ${i}
python snipar/simulate/simulate_trait_quad.py snipar/tests/test_data/sample${i}.bed outputs/tmp/t__t${i}_fams.ped 0.8 snipar/tests/test_data/h2_quad_0.8${i} --no_sib --dncor 0.5
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g snipar/tests/test_data/sample${i}.bgen -s snipar/tests/test_data/sample${i}.sample -og snipar/tests/test_data/sample_reduced${i}.bgen -os snipar/tests/test_data/sample_reduced${i}.sample -filetype bgen -ofiletype bgen -excl-samples outputs/tmp/t__t${i}_remove.txt
plink2 --bgen snipar/tests/test_data/sample_reduced${i}.bgen ref-last --sample snipar/tests/test_data/sample_reduced${i}.sample --make-bed --out snipar/tests/test_data/sample_reduced${i} --oxford-single-chr ${i}
cp outputs/tmp/t__t${i}.our.segments.gz snipar/tests/test_data/sample${i}.our.segments.gz
cp outputs/tmp/t__t${i}.king.segments.gz snipar/tests/test_data/sample${i}.king.segments.gz
cp outputs/tmp/t__t${i}.kingallsegs.txt snipar/tests/test_data/sample${i}.kingallsegs.txt
cp outputs/tmp/t__t${i}.direct_effects snipar/tests/test_data/sample${i}.direct_effects
cp outputs/tmp/t__t${i}.indirect_effects snipar/tests/test_data/sample${i}.indirect_effects
cp outputs/tmp/t__t${i}.father_phen snipar/tests/test_data/sample${i}.father_phen
cp outputs/tmp/t__t${i}.mother_phen snipar/tests/test_data/sample${i}.mother_phen
cp outputs/tmp/t__t${i}.sib_phen snipar/tests/test_data/sample${i}.sib_phen

awk 'NR % 4 == 1 || NR % 4 == 2 {print $1 "\t" $2}' snipar/tests/test_data/sample${i}.fam > snipar/tests/test_data/sample_sib_ids${i}
plink2 --keep-fam snipar/tests/test_data/sample_sib_ids${i}  --bfile snipar/tests/test_data/sample${i} --make-bed --out snipar/tests/test_data/sample_sib${i}
paste snipar/tests/test_data/sample_sib${i}.fam snipar/tests/test_data/sample${i}.sib_phen  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7}' > outputs/tmp/tmp
cp outputs/tmp/tmp snipar/tests/test_data/sample_sib${i}.fam

awk 'NR % 4 == 3 {print $1 "\t" $2}'  snipar/tests/test_data/sample${i}.fam > snipar/tests/test_data/sample_father_ids${i}
plink2 --keep-fam snipar/tests/test_data/sample_father_ids${i}  --bfile snipar/tests/test_data/sample${i} --make-bed --out snipar/tests/test_data/sample_father${i}
paste snipar/tests/test_data/sample_father${i}.fam snipar/tests/test_data/sample${i}.father_phen  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7}' > outputs/tmp/tmp
cp outputs/tmp/tmp snipar/tests/test_data/sample_father${i}.fam

awk 'NR % 4 == 0 {print $1 "\t" $2}'  snipar/tests/test_data/sample${i}.fam > snipar/tests/test_data/sample_mother_ids${i}
plink2 --keep-fam snipar/tests/test_data/sample_mother_ids${i}  --bfile snipar/tests/test_data/sample${i} --make-bed --out snipar/tests/test_data/sample_mother${i}
paste snipar/tests/test_data/sample_mother${i}.fam snipar/tests/test_data/sample${i}.mother_phen  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7}' > outputs/tmp/tmp
cp outputs/tmp/tmp snipar/tests/test_data/sample_mother${i}.fam

plink1 --bfile snipar/tests/test_data/sample_mother${i} --bmerge snipar/tests/test_data/sample_father${i} --out snipar/tests/test_data/sample_parent${i}
done

cp snipar/tests/test_data/sample1.agesex snipar/tests/test_data/sample.agesex
cp snipar/tests/test_data/sample1.ped snipar/tests/test_data/sample.ped
cp snipar/tests/test_data/sample1.king snipar/tests/test_data/sample.king
python -c "
import numpy as np
import pandas as pd
np.random.seed(100)
mask = np.random.choice(range(3000), 100, replace=False)
ped = pd.read_csv('snipar/tests/test_data/sample.ped', delim_whitespace=True)
ped.loc[ped['FID'].isin(mask)].to_csv('snipar/tests/test_data/pedigree_creation_sample.ped', sep = '\t', index = False)

agesex = pd.read_csv('snipar/tests/test_data/sample.agesex', delim_whitespace=True)
agesex.loc[agesex['FID'].isin(mask)].to_csv('snipar/tests/test_data/pedigree_creation_sample.agesex', sep = '\t', index = False)

king = pd.read_csv('snipar/tests/test_data/sample.king', delim_whitespace=True)
king.loc[king['FID1'].isin(mask)].to_csv('snipar/tests/test_data/pedigree_creation_sample.king', sep = '\t', index = False)
"
cp outputs/tmp/t__t1phen.npy snipar/tests/test_data/phen.npy
for fname in outputs/tmp/*.png; do cp ${fname} snipar/tests/test_data/${fname:17} ; done
for fname in outputs/tmp/*.png; do cp ${fname} graphs/${fname:17} ; done
cp outputs/tmp/t__t1rel_change*.png snipar/tests/test_data/rel_change.png
zcat snipar/tests/test_data/sample1.king.segments.gz > outputs/tmp/sample.king.segments
zcat snipar/tests/test_data/sample2.king.segments.gz | tail -n +2 >> outputs/tmp/sample.king.segments
gzip outputs/tmp/sample.king.segments
cp outputs/tmp/sample.king.segments.gz snipar/tests/test_data/

zcat snipar/tests/test_data/sample1.our.segments.gz > outputs/tmp/sample.our.segments
zcat snipar/tests/test_data/sample2.our.segments.gz | tail -n +2 >> outputs/tmp/sample.our.segments
gzip outputs/tmp/sample.our.segments
cp outputs/tmp/sample.our.segments.gz snipar/tests/test_data/

cat snipar/tests/test_data/sample1.kingallsegs.txt > snipar/tests/test_data/sample.kingallsegs.txt
cat snipar/tests/test_data/sample2.kingallsegs.txt | tail -n +2 >> snipar/tests/test_data/sample.kingallsegs.txt

echo "ID1     ID2     IBDType Chr     start_coordinate        stop_coordinate" > snipar/tests/test_data/sample.empty.segments
gzip snipar/tests/test_data/sample.empty.segments


plink1 --bfile snipar/tests/test_data/sample1 --bmerge snipar/tests/test_data/sample2 --out snipar/tests/test_data/sample1_2
plink2 --bfile snipar/tests/test_data/sample1_2 --pca 10  approx --out snipar/tests/test_data/sample1_2_pca
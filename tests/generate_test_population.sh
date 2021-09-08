rm outputs/tmp/*
rm test_data/*
for i in {1..2}
do
python example/simulate_pop.py 1000 0.5 3000 1000 1000 0.5 "outputs/tmp/t__t${i}" --chrom ${i}
python -c "import pandas as pd;
from sibreg.bin.preprocess_data import create_pedigree;
ped = pd.read_csv('outputs/tmp/t__t${i}_fams.ped', delim_whitespace=True);
ped.to_csv('test_data/sample_parent_sib${i}.ped', index = False, sep = ' ');
offsprings = ped[~ (ped['IID'].str.endswith('_M') | ped['IID'].str.endswith('_P'))]
offsprings.to_csv('test_data/sample_sib${i}.ped', index = False, sep = ' ');
remove = pd.read_csv('outputs/tmp/t__t1_remove.txt', delim_whitespace=True, names = ['counter', 'id']);
ped = ped[~ped['IID'].isin(remove['id'])];
ped.to_csv('test_data/sample${i}.ped', index = False, sep = ' ');
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
kinship_df.to_csv('test_data/sample${i}.king', sep = '\t', index = False);
"
cp outputs/tmp/t__t${i}.agesex test_data/sample${i}.agesex
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g outputs/tmp/t__t${i}.haps -s outputs/tmp/t__t${i}.sample -og test_data/sample${i}.bgen -os test_data/sample${i}.sample -filetype shapeit_haplotypes -ofiletype bgen
plink/plink2 --bgen test_data/sample${i}.bgen ref-last --sample test_data/sample${i}.sample --make-bed --out test_data/sample${i} --oxford-single-chr ${i}
python example/simulate_trait_quad.py test_data/sample${i}.bed outputs/tmp/t__t${i}_fams.ped 0.8 test_data/h2_quad_0.8${i} --no_sib --dncor 0.5
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g test_data/sample${i}.bgen -s test_data/sample${i}.sample -og test_data/sample_reduced${i}.bgen -os test_data/sample_reduced${i}.sample -filetype bgen -ofiletype bgen -excl-samples outputs/tmp/t__t${i}_remove.txt
plink/plink2 --bgen test_data/sample_reduced${i}.bgen ref-last --sample test_data/sample_reduced${i}.sample --make-bed --out test_data/sample_reduced${i} --oxford-single-chr ${i}
cp outputs/tmp/t__t${i}.segments.gz test_data/sample${i}.segments.gz
cp outputs/tmp/t__t${i}allsegs.txt test_data/sample${i}allsegs.txt
cp outputs/tmp/t__t${i}.direct_effects test_data/sample${i}.direct_effects
cp outputs/tmp/t__t${i}.indirect_effects test_data/sample${i}.indirect_effects
cp outputs/tmp/t__t${i}.father_phen test_data/sample${i}.father_phen
cp outputs/tmp/t__t${i}.mother_phen test_data/sample${i}.mother_phen
cp outputs/tmp/t__t${i}.sib_phen test_data/sample${i}.sib_phen

awk 'NR % 4 == 1 || NR % 4 == 2 {print $1 "\t" $2}' test_data/sample${i}.fam > test_data/sample_sib_ids${i}
plink/plink2 --keep-fam test_data/sample_sib_ids${i}  --bfile test_data/sample${i} --make-bed --out test_data/sample_sib${i}
paste test_data/sample_sib${i}.fam test_data/sample${i}.sib_phen  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7}' > outputs/tmp/tmp
cp outputs/tmp/tmp test_data/sample_sib${i}.fam

awk 'NR % 4 == 3 {print $1 "\t" $2}'  test_data/sample${i}.fam > test_data/sample_father_ids${i}
plink/plink2 --keep-fam test_data/sample_father_ids${i}  --bfile test_data/sample${i} --make-bed --out test_data/sample_father${i}
paste test_data/sample_father${i}.fam test_data/sample${i}.father_phen  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7}' > outputs/tmp/tmp
cp outputs/tmp/tmp test_data/sample_father${i}.fam

awk 'NR % 4 == 0 {print $1 "\t" $2}'  test_data/sample${i}.fam > test_data/sample_mother_ids${i}
plink/plink2 --keep-fam test_data/sample_mother_ids${i}  --bfile test_data/sample${i} --make-bed --out test_data/sample_mother${i}
paste test_data/sample_mother${i}.fam test_data/sample${i}.mother_phen  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7}' > outputs/tmp/tmp
cp outputs/tmp/tmp test_data/sample_mother${i}.fam

plink1 --bfile test_data/sample_mother${i} --bmerge test_data/sample_father${i} --out test_data/sample_parent${i}
done

cp test_data/sample1.agesex test_data/sample.agesex
cp test_data/sample1.ped test_data/sample.ped
cp test_data/sample1.king test_data/sample.king
cp outputs/tmp/t__t1phen.npy test_data/phen.npy
for fname in outputs/tmp/*.png; do cp ${fname} test_data/${fname:17} ; done
for fname in outputs/tmp/*.png; do cp ${fname} graphs/${fname:17} ; done
cp outputs/tmp/t__t1rel_change*.png test_data/rel_change.png
zcat test_data/sample1.segments.gz > outputs/tmp/sample.segments
zcat test_data/sample2.segments.gz | tail -n +2 >> outputs/tmp/sample.segments
gzip outputs/tmp/sample.segments
cp outputs/tmp/sample.segments.gz test_data/

cat test_data/sample1allsegs.txt > test_data/sampleallsegs.txt
cat test_data/sample2allsegs.txt | tail -n +2 >> test_data/sampleallsegs.txt


plink1 --bfile test_data/sample1 --bmerge test_data/sample2 --out test_data/sample1_2
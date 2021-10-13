# Simulate population genotypes
python example/simulate_pop_simple.py 1000 0.5 3000 1000 1000 0.5 example/sample --chrom 1
# Convert to bgen
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g example/sample.haps -s example/sample.sample -og example/sample.bgen -os example/sample.sample -filetype shapeit_haplotypes -ofiletype bgen
# Convert to bed
plink/plink2 --bgen example/sample.bgen ref-last --sample example/sample.sample --make-bed --out example/sample --oxford-single-chr 1
# Simulate phenotype
python example/simulate_trait_quad.py example/sample.bed example/sample.ped 0.8 example/h2_quad_0.8 --no_sib --dncor 0.5
# Remove some parents
/disk/genetics/ukb/alextisyoung/qctool/build/release/qctool_v2.0.7 -g example/sample.bgen -s example/sample.sample -og example/sample_reduced.bgen -os example/sample_reduced.sample -filetype bgen -ofiletype bgen -excl-samples example/sample_remove.txt
# Convert reduced to bed
plink/plink2 --bgen example/sample_reduced.bgen ref-last --sample example/sample_reduced.sample --make-bed --out example/sample_reduced --oxford-single-chr 1
# Clean up
rm example/sample.haps
mv example/sample_reduced.bgen example/sample.bgen
mv example/sample_reduced.sample example/sample.sample
mv example/sample_reduced.bed example/sample.bed
mv example/sample_reduced.fam example/sample.fam
mv example/sample_reduced.bim example/sample.bim
rm example/sample_remove.txt
rm example/sample.log
rm example/sample_reduced.log





python -c 'import pandas as pd;
from sibreg.bin.preprocess_data import create_pedigree;
ped = pd.read_csv("outputs/tmp/t__t_fams.ped", sep = " ");
remove = pd.read_csv("outputs/tmp/t__t_remove.txt", sep = " ", names = ["counter", "id"]);
ped = ped[~ped["IID"].isin(remove["id"])];
ped.to_csv("example/sample.ped", index = False, sep = " ");
sibships = ped.groupby(["FID","FATHER_ID", "MOTHER_ID"]).agg({"IID":list}).reset_index();
ids = set(ped["IID"].tolist());
kinship = [];
for index, row in sibships.iterrows():
    for i in range(1,len(row["IID"])):
        for j in range(i):
            kinship.append([row["FID"], row["IID"][i], row["FID"], row["IID"][j], "FS"]);
for index, row in ped.iterrows():
    if row["FATHER_ID"] in ids:
        kinship.append([row["FID"], row["IID"], row["FID"], row["FATHER_ID"],"PO"]);
    if row["MOTHER_ID"] in ids:
        kinship.append([row["FID"], row["IID"], row["FID"], row["MOTHER_ID"],"PO"]);
kinship_df = pd.DataFrame(kinship, columns = ["FID1", "ID1", "FID2", "ID2", "InfType"]);
kinship_df.to_csv("example/sample.king", sep = "\t", index = False);
agesex = ped.copy();
agesex["is_father"] = agesex["IID"].isin(agesex["FATHER_ID"]);
agesex["is_mother"] = agesex["IID"].isin(agesex["MOTHER_ID"]);
agesex["sex"] = "F";
agesex.loc[agesex["is_father"],"sex"] = "M";
agesex["age"] = 10;
agesex.loc[agesex["is_father"]|agesex["is_mother"],"age"] = 100;
agesex.to_csv("example/sample.agesex", sep = " ", index = False);
king_chr1 = pd.read_csv("outputs/tmp/t__t.segments.gz", sep = "\t");
king_chr2 = king_chr1.copy();
king_chr2["Chr"] = 2;
king_chr12 = pd.concat([king_chr1, king_chr2]);
king_chr12.to_csv("example/sample.segments.gz", sep = "\t", index = False);
print(king_chr12);
result = create_pedigree("example/sample.king",
                           "example/sample.agesex",
).sort_values(by=["FID", "IID"])
selected_people = [str(i*60)+"_0" for i in range(100)] + [str(i*60)+"_1" for i in range(100)] + [str(i*60)+"_P" for i in range(100)] + [str(i*60)+"_M" for i in range(100)]#set(np.array([[row["IID"], row["FATHER_ID"], row["MOTHER_ID"]] for index, row in result.iterrows() if int(row["FID"])%20==0]).reshape((1,-1))[0])
new_result = result[result["IID"].isin(selected_people)]
new_result.to_csv("example/pedigree_creation_sample.ped", sep = " ", index = False)
king = pd.read_csv("example/sample.king", sep = "\t")
new_king = king[king["ID1"].isin(new_result["IID"]) | king["ID2"].isin(new_result["IID"])]
new_king.to_csv("example/pedigree_creation_sample.king", sep="\t", index=False)
'

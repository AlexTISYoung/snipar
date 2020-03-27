python example/simulate_pop.py 1000 0.05 450 350 50 0 "test_data/t__t"
python -c 'import pandas as pd;
ped = pd.read_csv("test_data/t__t_fams.ped", names = ["FID", "IID", "FATHER_ID", "MOTHER_ID"], sep = " ");
remove = pd.read_csv("test_data/t__t_remove.txt", sep = " ", names = ["counter", "id"]);
ped = ped[~ped["IID"].isin(remove["id"])];
ped.to_csv("test_data/sample.ped", index = False, sep = " ")
sibships = ped.groupby(["FID","FATHER_ID", "MOTHER_ID"]).agg({"IID":list}).reset_index()
ids = set(ped["IID"].tolist())
kinship = []
for index, row in sibships.iterrows():
    for i in range(1,len(row["IID"])):
        for j in range(i):
            kinship.append([row["FID"], row["IID"][i], row["FID"], row["IID"][j], "FS"])
for index, row in ped.iterrows():
    if row["FATHER_ID"] in ids:
        kinship.append([row["FID"], row["IID"], row["FID"], row["FATHER_ID"],"PO"])
    if row["MOTHER_ID"] in ids:
        kinship.append([row["FID"], row["IID"], row["FID"], row["MOTHER_ID"],"PO"])
kinship_df = pd.DataFrame(kinship, columns = ["FID1", "ID1", "FID2", "ID2", "InfType"])
kinship_df.to_csv("test_data/sample.king", sep = "\t", index = False)
agesex = ped.copy()
agesex["is_father"] = agesex["IID"].isin(agesex["FATHER_ID"])
agesex["is_mother"] = agesex["IID"].isin(agesex["MOTHER_ID"])
agesex["sex"] = "F"
agesex.loc[agesex["is_father"],"sex"] = "M"
agesex["age"] = 10
agesex.loc[agesex["is_father"]|agesex["is_mother"],"age"] = 100
agesex.to_csv("test_data/sample.agesex", sep = " ", index = False)
'
mv test_data/t__t.segments.gz test_data/sample.segments.gz
plink/plink --file test_data/t__t --make-bed --out test_data/t__t
plink/plink --bfile test_data/t__t --remove test_data/t__t_remove.txt --make-bed --out test_data/sample1
cp test_data/sample1.bed test_data/sample2.bed 
cp test_data/sample1.bim test_data/sample2.bim 
cp test_data/sample1.fam test_data/sample2.fam
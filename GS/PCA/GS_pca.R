setwd('/disk/genetics/sibling_consortium/GS20k/alextisyoung/HM3/pca')
library(bigsnpr)
library(ggplot2)
library(RColorBrewer)

## Ref: https://privefl.github.io/bigsnpr/articles/bedpca.html

##### Read bed files #####
obj.bed = bed('data/1000G_hm3.bed') # 1KG samples on HapMap3 SNPs

obj.gsbed = bed('../../../aokbay/imputed/HRC/plink/HM3/GS_HRC_HM3.bed') # GS samples on HM3 SNPs

#### Remove relatives from 1kG ####

ind.rel <- match(c(rel$IID1, rel$IID2), obj.bed$fam$sample.ID)  # /!\ use $ID1 instead with old PLINK
ind.norel <- rows_along(obj.bed)[-ind.rel]

##### Projection onto 1kg PCs #####
gs_proj_1kg = bed_projectPCA(obj.bed.ref=obj.bed,obj.bed.new=obj.gsbed,k=20,
                            ind.row.ref=ind.norel, ncores=20)

##### Get superpopulation labels ######

pop_labels = read.table('data/igsr_samples.tsv',sep='\t',header=T)

### Plot GS on 1kg PCs
ref_pcs <- predict(gs_proj_1kg$obj.svd.ref)
dimnames(ref_pcs)[[1]] = obj.bed$.fam[ind.norel,'sample.ID']
ref_superpops = pop_labels[match(dimnames(ref_pcs)[[1]],pop_labels$Sample.name),'Superpopulation.code']

superpops = levels(as.factor(ref_superpops))
legend_labs = c(superpops,'GS')
colours = brewer.pal(length(legend_labs), "Set4")
superpops = cbind(superpops,colours[1:length(superpops)])

pdf('GS_1kg_PCA_proj.pdf',width=10,height=10)
plot(ref_pcs[,1],ref_pcs[,2],col=superpops[match(ref_superpops,superpops[,1]),2],xlab='PC1',ylab='PC2')
legend(y=-50,x=-200,legend=legend_labs,col=colours,pch=c(rep(1,length(superpops)),2))
points(gs_proj_1kg$OADP_proj[,1],gs_proj_1kg$OADP_proj[,2],col=colours[length(legend_labs)],pch=2)
dev.off()


### KNN ###
library(class)
has_label = !is.na(ref_superpops)
knn_1kg = knn(train=ref_pcs[has_label,],test=gs_proj_1kg$OADP_proj,cl=ref_superpops[has_label],prob=TRUE,k=20,l=20)

gs_pop = data.frame(IID=obj.gsbed$.fam$sample.ID,
                    pop=knn_1kg, prob=attributes(knn_1kg)$prob)

write.table(gs_pop,'GS_population_assignment.txt',quote=F,row.names=F)

### Compute PCs ###
setwd('/disk/genetics/sibling_consortium/GS20k/alextisyoung/grandpar/')
library(bigsnpr)
obj.gsbed = bed('haplotypes/bedfiles/autosome.bed') # GS samples on HM3 SNPs

king_unrel = read.table('king_unrelunrelated.txt')
king_unrel = sapply(king_unrel[,2],function(x) strsplit(x,'->')[[1]][2])
ind.norel <- match(c(king_unrel), obj.gsbed$fam$sample.ID)  # /!\ use $ID1 instead with old PLINK
ind.norel=ind.norel[!is.na(ind.norel)]
obj.svd <- bed_autoSVD(obj.gsbed, ind.row = ind.norel, k = 20,
                       ncores = nb_cores())

PCs <- matrix(NA, nrow(obj.gsbed), ncol(obj.svd$u))
PCs[ind.norel, ] <- predict(obj.svd)
proj <- bed_projectSelfPCA(obj.svd, obj.gsbed, ind.row=rows_along(obj.gsbed)[-ind.norel])
PCs[-ind.norel, ] <- proj$OADP_proj
dimnames(PCs)[[1]] = obj.gsbed$fam$sample.ID
dimnames(PCs)[[2]] = paste('PC',1:20,sep='')
write.table(PCs,'GS_PCs.txt',quote=F)
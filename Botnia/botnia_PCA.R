setwd('/ludc/Active_Projects/BOTNIA_AYoung_analysis/Private/haplotypes/bedfiles/pca/')
install.packages('bigsnpr')
library(bigsnpr)

## Ref: https://privefl.github.io/bigsnpr/articles/bedpca.html

bedfile <- download_1000G("data")
##### Read bed files #####
obj.bed = bed("data/1000G_phase3_common_norel.bed") # 1KG samples on HapMap3 SNPs
# Botnia samples genotypes
obj.gsbed = bed('/ludc/Active_Projects/BOTNIA_AYoung_analysis/Private/haplotypes/bedfiles/autosome.bed') 


##### Projection onto 1kg PCs #####
gs_proj_1kg = bed_projectPCA(obj.bed.ref=obj.bed,obj.bed.new=obj.gsbed,k=20, ncores=20)

##### Get superpopulation labels ######
pop_labels = read.table('data/igsr_samples.tsv.1',sep='\t',header=T)

### put reference onto inferred PCs and label with superpop
ref_pcs <- predict(gs_proj_1kg$obj.svd.ref)
dimnames(ref_pcs)[[1]] = obj.bed$.fam[,'sample.ID']
ref_superpops = pop_labels[match(dimnames(ref_pcs)[[1]],pop_labels$Sample.name),'Superpopulation.code']

#superpops = levels(as.factor(ref_superpops))
#legend_labs = c(superpops,'GS')
#colours = brewer.pal(length(legend_labs), "Set4")
#superpops = cbind(superpops,colours[1:length(superpops)])

#pdf('GS_1kg_PCA_proj.pdf',width=10,height=10)
#plot(ref_pcs[,1],ref_pcs[,2],col=superpops[match(ref_superpops,superpops[,1]),2],xlab='PC1',ylab='PC2')
#legend(y=-50,x=-200,legend=legend_labs,col=colours,pch=c(rep(1,length(superpops)),2))
#points(gs_proj_1kg$OADP_proj[,1],gs_proj_1kg$OADP_proj[,2],col=colours[length(legend_labs)],pch=2)
#dev.off()


### KNN ###
library(class)
has_label = !is.na(ref_superpops)
knn_1kg = knn(train=ref_pcs[has_label,],test=gs_proj_1kg$OADP_proj,cl=ref_superpops[has_label],prob=TRUE,k=20,l=20)

inferred_pop = data.frame(IID=obj.gsbed$.fam$sample.ID,
                    pop=knn_1kg, prob=attributes(knn_1kg)$prob)

write.table(inferred_pop,'population_assignment.txt',quote=F,row.names=F)
non_eur = inferred_pop[inferred_pop$pop!='EUR' | inferred_pop$prob<1,1]
# match to first iteration .bgen ids
bgen_ids = read.table('../../chr_22.sample')
bgen_indices = match(non_eur[,1],bgen_ids[,2])-1
bgen_non_eur = paste('sample',bgen_indices,sep='_')
write.table(non_eur,'non_european_samples_bgen.txt',quote=F,row.names=F,col.names=F)

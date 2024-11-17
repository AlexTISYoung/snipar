library(ggplot2)
library(dplyr)
library(latex2exp)


# Non-sampling variance

readfile <- read.csv("paper_materials/sp_non_sampling_var_ref_to_standardGWAS_at_Fst0.001.csv", sep='\t')
readfile$method <- mapvalues(readfile$method, 
          from=c("standard", 'sib diff', 'robust'), 
          to=c("standard GWAS", 'robust / sib-difference', 'non-transmitted'))
ord <- c('robust / sib-difference', 'non-transmitted', 'Young et al. 2022', 'unified', 'standard GWAS')
my_colors <- c("deepskyblue", "deeppink", "darkslategray3", "#569756", "darksalmon")
readfile <- readfile %>%  mutate(method =  factor(method, levels = ord)) %>%
  arrange(method)
options(repr.plot.width=12, repr.plot.height=8)
readfile[,c(1,3,2,4)]
ggplot(readfile, aes(x=factor(Fst), y=non_sampling_var, fill=method)) + 
    geom_point(aes(color=method),position=position_dodge(width=0.75), size=1.5)+
    scale_color_manual(values = my_colors) +
    geom_errorbar(aes(ymin=non_sampling_var-non_sampling_var_std*1.96, ymax=non_sampling_var+non_sampling_var_std*1.96,
  color=method), width=0.5,
                 position=position_dodge(0.75),size=1) +
            labs(x=TeX("$F_{st}$"), y =TeX('bias (non-sampling variance)')) +
            geom_hline(yintercept=0, linetype="dashed", 
                color = "black",size=0.2) +
            theme(text=element_text(size=15), legend.position=c(0.2, 0.8),legend.box.background = element_rect(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.8), legend.key = element_rect(fill = "white"),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave('paper_materials/sp_non_sampling_vs_Fst.pdf', width = 8,height = 6)



readfile <- readfile[which(readfile$method!='standard GWAS'),]
options(repr.plot.width=12, repr.plot.height=8)
readfile[,c(1,3,2,4)]
ggplot(readfile, aes(x=factor(Fst), y=non_sampling_var, fill=method)) + 
    geom_point(aes(color=method),position=position_dodge(width=0.75), size=1.5)+
    scale_color_manual(values = my_colors) +
    geom_errorbar(aes(ymin=non_sampling_var-non_sampling_var_std*1.96, ymax=non_sampling_var+non_sampling_var_std*1.96,
  color=method), width=0.5,
                 position=position_dodge(0.75),size=1) +
            labs(x=TeX("$F_{st}$"), y = TeX('bias (non-sampling variance)'))
            geom_hline(yintercept=0, linetype="dashed", 
                color = "black",size=0.2) +
            theme(text=element_text(size=15), 
           panel.border = element_rect(colour = "black", fill=NA, size=0.8),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
legend.position=c(0.2, 0.8),legend.box.background = element_rect(colour = "black"),legend.key = element_rect(fill = "white"))
ggsave('paper_materials/sp_non_sampling_vs_Fst_zoomed.pdf', width = 8,height = 6)



# Z^2

readfile <- read.csv("paper_materials/sp_Z_var.csv", sep='\t')
readfile$method <- mapvalues(readfile$method, 
          from=c("standard", 'sib diff', 'robust'), 
          to=c("standard GWAS", 'robust / sib-difference', 'non-transmitted'))
ord <- c('robust / sib-difference', 'non-transmitted', 'Young et al. 2022', 'unified', 'standard GWAS')
readfile <- readfile %>%  mutate(method =  factor(method, levels = ord)) %>%
  arrange(method)
options(repr.plot.width=12, repr.plot.height=8)
readfile[,c(1,3,2,4)]
ggplot(readfile, aes(x=factor(Fst), y=non_sampling_var, fill=method)) + 
    geom_point(aes(color=method),position=position_dodge(width=0.75), size=1.5)+
    scale_color_manual(values = my_colors) +
    geom_errorbar(aes(ymin=non_sampling_var-non_sampling_var_std*1.96, ymax=non_sampling_var+non_sampling_var_std*1.96,
  color=method), width=0.5,
                 position=position_dodge(0.75),size=1) +
            labs(x=TeX("$F_{st}$"), y = TeX("mean $Z^2$")) +
            geom_hline(yintercept=1, linetype="dashed", 
                color = "black",size=0.2) +
            theme(text=element_text(size=15), legend.position=c(0.2, 0.8),legend.box.background = element_rect(colour = "black"),
           panel.border = element_rect(colour = "black", fill=NA, size=0.8),legend.key = element_rect(fill = "white"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave('paper_materials/sp_Zvar_vs_Fst.pdf', width = 8,height = 6)

readfile <- readfile[which(readfile$method!='standard GWAS'),]
options(repr.plot.width=12, repr.plot.height=8)
readfile[,c(1,3,2,4)]
ggplot(readfile, aes(x=factor(Fst), y=non_sampling_var, fill=method)) + 
    geom_point(aes(color=method),position=position_dodge(width=0.75), size=1.5)+
    scale_color_manual(values = my_colors) +
    geom_errorbar(aes(ymin=non_sampling_var-non_sampling_var_std*1.96, ymax=non_sampling_var+non_sampling_var_std*1.96,
  color=method), width=0.5,
                 position=position_dodge(0.75),size=1) +
            labs(x=TeX("$F_{st}$"), y = TeX("mean $Z^2$")) +
            geom_hline(yintercept=1, linetype="dashed", 
                color = "black",size=0.2) +
            theme(text=element_text(size=15), legend.position=c(0.2, 0.8),legend.box.background = element_rect(colour = "black"),
           panel.border = element_rect(colour = "black", fill=NA, size=0.8),legend.key = element_rect(fill = "white"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave('paper_materials/sp_Zvar_vs_Fst_zoomed.pdf', width = 8,height = 6)

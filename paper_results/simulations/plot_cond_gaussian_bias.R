library(ggplot2)
library(dplyr)
library(latex2exp)
require(plyr)
readfile <- read.csv("impute_cond_gau_results/sp_non_sampling_var.csv", sep='\t')

my_colors <- c("deepskyblue", "deeppink", "darkslategray3", "darkseagreen4", "darksalmon")
ggplot(readfile, aes(x=factor(Fst), y=non_sampling_var, fill=method)) + 
    geom_point(aes(color=method),position=position_dodge(width=0.75), size=1.5)+
    scale_color_manual(values = my_colors) +
    geom_errorbar(aes(ymin=non_sampling_var-non_sampling_var_std*1.96, ymax=non_sampling_var+non_sampling_var_std*1.96,
  color=method), width=0.25,
                 position=position_dodge(0.75),size=1) +
            labs(x=TeX("$F_{st}$"), y = TeX("Mean $Z^2$")) +
            geom_hline(yintercept=1, linetype="dashed", 
                color = "black",size=0.2) +
               scale_y_continuous(breaks=c(1,50,100,150)) + 
            theme(text=element_text(size=15), legend.position=c(0.2, 0.8),legend.box.background = element_rect(colour = "black"),
           panel.border = element_blank(), legend.key = element_rect(fill = "white"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
axis.line = element_line(colour = "black"))
ggsave('impute_cond_gau_results/mean_Z2_vs_Fst.pdf',width = 8,height = 6)
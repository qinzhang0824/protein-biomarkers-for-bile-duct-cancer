library(rio)
library(ggplot2)
library(ggsci)
library("ggpubr")
library(forestplot)

mr.or <-fread("final results/Discovery cohort/mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls",sep='\t',header=T)
mr.or1 <-fread('final results/Validation cohort 1/mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls',header = T)
mr.or2 <-fread('final results/Validation cohort 2/mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls',header = T)

out_multi.finn <-mr.or[mr.or$pval<(0.05/1178),]
out_multi.4915 <-mr.or1[mr.or1$id.exposure %in% c("NCAN","SERPINA1"),]
out_multi.g3859 <-mr.or2[mr.or2$id.exposure %in% c("NCAN","SERPINA1"),]

out_multi <-rbind(out_multi.finn,out_multi.4915,out_multi.g3859)

hz <- paste(round(out_multi$or,3),
            "(",round(out_multi$or_lci95,3),
            "-",round(out_multi$or_uci95,3),")",sep = "")


tabletext <- cbind(c(NA,"Exposure",out_multi$id.exposure),
                   c(NA,"Outcome",out_multi$outcome),
                   c(NA,"Methods",out_multi$method),
                   c(NA,"P value",ifelse(out_multi$pval<0.001,"P < 0.001",round(out_multi$pval,3))),
                   #c(NA,"P value",out_multi$pval),
                   c(NA,"OR(95% CI)",hz))


pforest <- forestplot(labeltext=tabletext, 
                      graph.pos=4
                      col=fpColors(box="#00468B", lines="black", zero = "red"),
                      mean=c(NA,NA,out_multi$or),
                      lower=c(NA,NA,out_multi$or_lci95), 
                      upper=c(NA,NA,out_multi$or_uci95), 
                      xlab="Odds Ratio(95%CI)",
                      boxsize=0.2,lwd.ci=0.8,  
                      ci.vertices.height = 0.08,ci.vertices=TRUE,
                      zero=1,lwd.zero=1,     
                      colgap=unit(14,"mm"),   
                      xticks = c(0.3,0.5,0.6,0.7,0.8,0.9,1,1.1,1.3,1.5,5,7), 
                      lwd.xaxis=2,         
                      lineheight = unit(1,"cm"), 
                      graphwidth = unit(0.2,"npc"), 
                      #cex=0.9, fn.ci_norm = fpDrawCircleCI, 
                      hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                                      "3" = gpar(lwd=2, col="black")),           
                      clip=c(0,3),
                      txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                     ticks=gpar(cex=1.1),
                                     xlab=gpar(cex = 1.4),
                                     title=gpar(cex = 1.2))
)
pdf("final results/Figure3_Plasma_Protein_vs_Discovery_and_validation_Forest.pdf", width=12, height=6)
pforest
dev.off()

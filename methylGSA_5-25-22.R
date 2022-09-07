if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sesame")

BiocManager::install("methylGSA")
BiocManager::install("tidyverse")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(methylGSA)
library(tidyverse)
library(sesame)


getwd()
setwd("../")
getwd()

res = read.csv("./5-20-22/DML_res5-17-22.csv")

head(res)
cpg.pval <- res$Pval_groupMDD
head(cpg.pval)
names(cpg.pval)= res$Est_Probe_ID ### convert the data frame into a vector 
head(cpg.pval)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

res1 = methylglm(cpg.pval = cpg.pval, minsize = 200, 
                 maxsize = 500, GS.type = "KEGG") ## repeat what I did on the 5th floor computer, worked.

#### methylgometh: methylgometh calls gometh or gsameth function in missMethyl package to adjust number 
####### of CpGs in gene set testing. gometh modifiesgoseq method by fitting a probability weighting 
####### function and resampling from Wallenius non-central hypergeometric distribution. 

res_gometh_min100_sig0.001 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001, topDE = NULL, array.type = "EPIC",
                    GS.type = "GO", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_min100_sig0.001) ### nothing sig
write.csv(res_gometh_min100_sig0.001, file = "./res_gometh_min100_sig0.001.csv")


res_gometh_GO_min100_sig0.05 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.05, topDE = NULL, array.type = "EPIC",
                                          GS.type = "GO", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_GO_min100_sig0.05)
nrow(res_gometh_GO_min100_sig0.05[res_gometh_GO_min100_sig0.05$padj<0.05,]) ## 554 significant
#write.csv(res_gometh_GO_min100_sig0.05, file = "./methyGSA/res_gometh_GO_min100_sig0.05.csv")


res_gometh_GO_min100_sig0.01 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.01, topDE = NULL, array.type = "EPIC",
                                            GS.type = "GO", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)

head(res_gometh_GO_min100_sig0.01)
nrow(res_gometh_GO_min100_sig0.01[res_gometh_GO_min100_sig0.01$padj<0.05,]) ## 1 significant
write.csv(res_gometh_GO_min100_sig0.01, file = "./methyGSA/res_gometh_GO_min100_sig0.01.csv")


res_gometh_min100_top20 = methylgometh(cpg.pval = cpg.pval,  topDE =20, array.type = "EPIC",
                                          GS.type = "GO", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_min100_top20) ## one significant
write.csv(res_gometh_min100_top20, file = "./methyGSA/res_gometh_min100_top20.csv")

res_gometh_min100_top10 = methylgometh(cpg.pval = cpg.pval,  topDE =10, array.type = "EPIC",
                                       GS.type = "GO", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)

head(res_gometh_min100_top10) ### nothing sig


res_gometh_min100_top5 = methylgometh(cpg.pval = cpg.pval,  topDE =5, array.type = "EPIC",
                                       GS.type = "GO", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_min100_top5) ###n nothing sig


res_gometh_min100_top50 = methylgometh(cpg.pval = cpg.pval,  topDE =50, array.type = "EPIC",
                                      GS.type = "GO", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)

head(res_gometh_min100_top50) ### nothing sig

#####
res_gometh_kegg_min100_sig0.001 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001, topDE = NULL, array.type = "EPIC",
                                            GS.type = "KEGG", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_kegg_min100_sig0.001) ### nothing significant

res_gometh_kegg_min100_sig0.01 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.01, topDE = NULL, array.type = "EPIC",
                                               GS.type = "KEGG", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)

head(res_gometh_kegg_min100_sig0.01) ### 1 significant 05022 nuerodegenerative

write.csv(res_gometh_kegg_min100_sig0.01, file ="./methyGSA/res_gometh_kegg_min100_sig0.01.csv")


res_gometh_kegg_min100_sig0.05 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.05, topDE = NULL, array.type = "EPIC",
                                              GS.type = "KEGG", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)

head(res_gometh_kegg_min100_sig0.05) 
nrow(res_gometh_kegg_min100_sig0.05[res_gometh_kegg_min100_sig0.05$padj<0.05,]) ### 52 significant

write.csv(res_gometh_kegg_min100_sig0.05, file ="./methyGSA/res_gometh_kegg_min100_sig0.05.csv")


res_gometh_kegg_min100_top20= methylgometh(cpg.pval = cpg.pval, topDE = 20, array.type = "EPIC",
                                              GS.type = "KEGG", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_kegg_min100_top20) ### nothing significant


res_gometh_kegg_min100_top50= methylgometh(cpg.pval = cpg.pval, topDE = 50, array.type = "EPIC",
                                           GS.type = "KEGG", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_kegg_min100_top50) ### nothing sig

res_gometh_kegg_min100_top5= methylgometh(cpg.pval = cpg.pval, topDE = 5, array.type = "EPIC",
                                           GS.type = "KEGG", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_kegg_min100_top5) ### nothing sig

###

res_gometh_reactome_min100_sig0.001 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001, topDE = NULL, array.type = "EPIC",
                                              GS.type = "Reactome", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_reactome_min100_sig0.001) ### nothing significant


res_gometh_reactome_min100_sig0.01 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.01, topDE = NULL, array.type = "EPIC",
                                                   GS.type = "Reactome", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_reactome_min100_sig0.01) ### 2 significant
write.csv(res_gometh_reactome_min100_sig0.01, file = "./methyGSA/res_gometh_reactome_min100_sig0.01.csv")



res_gometh_reactome_min100_sig0.05 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.05, topDE = NULL, array.type = "EPIC",
                                                  GS.type = "Reactome", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_reactome_min100_sig0.05)
nrow(res_gometh_reactome_min100_sig0.05[res_gometh_reactome_min100_sig0.05$padj<0.05,]) ## 107 significant

##
res_gometh_reactome_min100_top20 = methylgometh(cpg.pval = cpg.pval, topDE = 20, array.type = "EPIC",
                                                  GS.type = "Reactome", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_reactome_min100_top20) ### nothing sig


res_gometh_reactome_min100_top50 = methylgometh(cpg.pval = cpg.pval, topDE = 50, array.type = "EPIC",
                                                GS.type = "Reactome", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)
head(res_gometh_reactome_min100_top50) ### nothing sig

### 

res_gometh_reactome_min100_top100 = methylgometh(cpg.pval = cpg.pval, topDE = 100, array.type = "EPIC",
                                                GS.type = "Reactome", GS.idtype = "SYMBOL", minsize = 100, maxsize = 500)


head(res_gometh_reactome_min100_top100) ### nothing sig





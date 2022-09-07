###### 4-7-22 Sesame

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sesame")

browseVignettes("sesame")

library(sesame)
getwd()
load("sesame_se.Rdata")

smry = DML(sesame_se$betas, ~sesame_se$group) ## DIDN'T WORK!)

### Install new version

BiocManager::install("zwdzwd/sesame")

package.version("sesame")
library(sesame)
library(SummarizedExperiment)
### test data
sesameDataCacheAll()

data <- sesameDataGet('HM450.76.TCGA.matched')
smry <- DML(data$betas[1:1000,], ~type, meta=data$sampleInfo)
test_res = summaryExtractTest(smry)
head(test_res)

### my data
#smry = DML(sesame_se@colData@listData$betas, ~sesame_se@colData@listData$group) 
se2 = SummarizedExperiment(assays=list(betas=as.matrix(assays(sesame_se)$betas)),
                           colData=colData(sesame_se))
smry <- DML(se2[1:1000], ~group)
res1000 = summaryExtractTest(smry)
head(res1000)

smry <- DML(se2, ~group, mc.cores = 8)
res= summaryExtractTest(smry)
head(res)

write.table(res, file = "./dml.res.group.csv", row.names = TRUE, sep = ",", quote = F)

smry2 = DML(se2, ~group+age+sex, mc.cores = 8)
res2= summaryExtractTest(smry2)
head(res2)

write.table(res2, file = "./dml_res_group_w_sexage.csv", row.names = TRUE, sep = ",", quote = F)

#colnames(attr(smry2), "model.matrix")

##merged_segs <- DMR(se2, smry)  didn't work!


###########  clean up the above code a little bit

library(cgwtools)
library(tidyverse)
getwd()
pdata = read.csv("pdata61.csv")
load("beta61_sesame.RData")
lsdata("beta61_sesame.RData")
head(pdata)

pdata = pdata %>% column_to_rownames(., var = "samples")
head(pdata)
class(pdata$group)
pdata$group = as.factor(pdata$group)

### remove NAs, DML doesn't work when have NAs
keep = rowSums(is.na(beta_61)) ==0

beta_61_fil = beta_61[keep,]
se = SummarizedExperiment(assays=list(betas=as.matrix(beta_61_fil)),
                          colData=pdata)
head(se)
save(se, file = "./se.RData")
smry = DML(se, ~group, mc.cores = 8)

res= summaryExtractTest(smry)
head(res)

### Add covariates
smry2 = DML(se, ~group+age+sex, mc.cores = 8)

res2= summaryExtractTest(smry2)
head(res2)

df = res2 %>% as.data.frame()
head(df)
df = df[order(df$Pval_groupMDD),]
head(df)
nrow(df)
sum(df$Pval_groupMDD <0.05)

### plots
library(ggplot2)
df = data.frame(Age = colData(se)$age,
                BetaValue = assay(se)[rownames(res2)[nrow(res2)],])
ggplot(df, aes(Age, BetaValue)) + geom_smooth(method="lm") + geom_point()

#### check if DNA methylation different between MDD and CT
ggplot(res2) + geom_point(aes(Est_groupMDD, -log10(Pval_groupMDD)))

### DMR

dmContrasts(smry2)   ### show contrast

merged = DMR(se, smry2, "groupMDD")  

head(merged)
sum(merged$Seg_Pval_adj <0.05)
save(merged, file = "./DMR_res.RData")
write.table(merged, file = "./DMR_res.csv", row.names = TRUE, sep = ",", quote = F)



#####
load("./sesame_se.Rdata")
rm(sesame_se)

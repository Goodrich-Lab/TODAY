# R code for burden score calculator ----
# Load necessary libraries
library(nhanesA)  # For accessing NHANES data
library(ltm)      # For latent trait models
library(Hmisc)    # For statistical methods and data analysis

# Set a seed for reproducibility
set.seed(12345)

# Decile calculation using isomers ------
# Fetch NHANES dataset for PFAS_J
fr.data = nhanes('PFAS_J')
fr.data.labels = nhanes('PFAS_J', includelabels = TRUE)
# Select relevant columns from the dataset
fr.data = fr.data[, c(1, 2, seq(4, 20, by = 2))]
fr.item = fr.data
# Remove rows with missing data
fr.item.complete = fr.item[which(complete.cases(fr.item)),]
# Reverse coding for the PFAS data (1 - original value)
fr.item.complete[, 3:11] = 1 - fr.item.complete[, 3:11]

# Exclude pregnant women from the analysis
pregnancy.data = nhanes('UCPREG_J')
merge.fr.p = merge(fr.item.complete, pregnancy.data, by = "SEQN", all = FALSE)
# Identify and remove pregnant women's data
sub.fr.p = fr.item.complete[-c(which(merge.fr.p$URXPREG == 1)), ]

# Fetch and prepare PFAS concentration data
pfas.cont = nhanes('PFAS_J')
pfas.cont = pfas.cont[, c(1, seq(3, 19, by = 2))]
sub.fr.p = merge(x = sub.fr.p, y = pfas.cont, by = "SEQN", all.x = TRUE, all.y = FALSE)

# Calculate weighted quantiles for different PFAS compounds
wtd_quantile_PAH <- wtd.quantile(sub.fr.p$LBXMPAH, weights = sub.fr.p$WTSB2YR, 
                                 probs = c(0, 0.4, 0.7, 0.9, 1), normwt = FALSE)

# Repeat the same process for other compounds (FUA, PFDE, etc.)
wtd_quantile_PAH <-wtd.quantile(sub.fr.p$LBXMPAH, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, 0.4, 0.7, 0.9, 1), normwt=FALSE)
wtd_quantile_FUA<-wtd.quantile(sub.fr.p$LBXPFUA, weights=sub.fr.p$WTSB2YR,
                               probs=c(0, 0.4, 0.7, 0.9, 1), normwt=FALSE)
wtd_quantile_PFDE<-wtd.quantile(sub.fr.p$LBXPFDE, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, 0.1, 0.4, 0.8, 0.9, 1), normwt=FALSE)

wtd_quantile_PFDE<-wtd.quantile(sub.fr.p$LBXPFDE, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, .5, 1), normwt=FALSE)

wtd_quantile_PFHS<-wtd.quantile(sub.fr.p$LBXPFHS, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
wtd_quantile_PFNA<-wtd.quantile(sub.fr.p$LBXPFNA, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
wtd_quantile_NFOA<-wtd.quantile(sub.fr.p$LBXNFOA, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
wtd_quantile_NFOS<-wtd.quantile(sub.fr.p$LBXNFOS, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
wtd_quantile_MFOS<-wtd.quantile(sub.fr.p$LBXMFOS, weights=sub.fr.p$WTSB2YR,
                                probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)


# Assign infinite values to the last quantile and negative values to the first quantile for each compound
# This is likely done to handle edge cases in the data
wtd_quantile_PAH[length(wtd_quantile_PAH)]<-Inf
wtd_quantile_FUA[length(wtd_quantile_FUA)]<-Inf
wtd_quantile_PFDE[length(wtd_quantile_PFDE)]<-Inf
wtd_quantile_PFHS[length(wtd_quantile_PFHS)]<-Inf
wtd_quantile_PFNA[length(wtd_quantile_PFNA)]<-Inf
wtd_quantile_NFOA[length(wtd_quantile_NFOA)]<-Inf
wtd_quantile_NFOS[length(wtd_quantile_NFOS)]<-Inf
wtd_quantile_MFOS[length(wtd_quantile_MFOS)]<-Inf

wtd_quantile_PAH[1]<- -0.1
wtd_quantile_FUA[1]<- -0.1
wtd_quantile_PFDE[1]<- -0.1
wtd_quantile_PFHS[1]<- -0.1
wtd_quantile_PFNA[1]<- -0.1
wtd_quantile_NFOA[1]<- -0.1
wtd_quantile_NFOS[1]<- -0.1
wtd_quantile_MFOS[1]<- -0.1

# Categorize data into deciles based on the quantiles calculated
sub.fr.p$PAH = cut(sub.fr.p$LBXMPAH, wtd_quantile_PAH,right=FALSE, 
                   labels=c(1:4))
sub.fr.p$PFUA = cut(sub.fr.p$LBXPFUA, wtd_quantile_FUA,right=FALSE,
                    labels=c(1:4))
sub.fr.p$PFDE = cut(sub.fr.p$LBXPFDE, wtd_quantile_PFDE,right=FALSE, 
                    labels = c(1:5))
sub.fr.p$PFHS = cut(sub.fr.p$LBXPFHS, wtd_quantile_PFHS,right=FALSE,
                    labels = c(1:10))
sub.fr.p$PFNA = cut(sub.fr.p$LBXPFNA, wtd_quantile_PFNA,right=FALSE,
                    labels = c(1:9))
sub.fr.p$NFOA = cut(sub.fr.p$LBXNFOA, wtd_quantile_NFOA,right=FALSE,
                    labels = c(1:10))
sub.fr.p$NFOS = cut(sub.fr.p$LBXNFOS, wtd_quantile_NFOS,right=FALSE,
                    labels = c(1:10))
sub.fr.p$MFOS = cut(sub.fr.p$LBXMFOS, wtd_quantile_MFOS,right=FALSE,
                    labels = c(1:10))

# Analyze Sb-PFOA compound
# Creating a new category based on detection/non-detection
sub.fr.p$BFOA <- sub.fr.p$LBDBFOAL + 1

# Prepare a new matrix for IRT models and perform the analysis
newmat.d = sub.fr.p[, c(21:29, 1)]
names(newmat.d) = c("me-PFOSA-AcOH", "PFUA", "PFDeA", "PFHxS", "PFNA", "nPFOA", "n-PFOS", "Sm-PFOS", "Sb-PFOA", "SEQN")
irt.grm.d0 = grm(newmat.d[, 1:9])
irt.grm.d = grm(newmat.d[,1:8]) # used this model, without Sb-PFOA
eap.d <- ltm::factor.scores(irt.grm.d, method="EAP", resp.patterns =
                              newmat.d[,1:8])$score.dat$z1
summary(eap.d)
score.d = cbind(newmat.d, eap.d)
score.d = score.d[,c("SEQN", "eap.d")]



# # Decile summing isomers -----
# fr.data=nhanes('PFAS_J')
# fr.data=fr.data[, c(1, 2, seq(4, 20, by=2))]
# fr.item = fr.data
# fr.item.complete = fr.item[which(complete.cases(fr.item)),]
# # 204 missing PFAS data
# fr.item.complete[, 3:11] = 1 - fr.item.complete[,3:11]
# # Exclude pregnant women (14)
# pregnancy.data=nhanes('UCPREG_J')
# merge.fr.p = merge(fr.item.complete, pregnancy.data, by = "SEQN", all=FALSE)
# which(merge.fr.p$URXPREG == 1)
# merge.fr.p$SEQN[which(merge.fr.p$URXPREG == 1)]
# sub.fr.p = fr.item.complete[-c(which(merge.fr.p$URXPREG == 1)), ]
# pfas.cont = nhanes('PFAS_J')
# pfas.cont = pfas.cont[,c(1, seq(3, 19, by=2))]
# sub.fr.p = merge(x=sub.fr.p, y=pfas.cont, by="SEQN", all.x=TRUE, all.y = FALSE)
# sub.fr.p$PFOA = sub.fr.p$LBXNFOA + sub.fr.p$LBXBFOA
# sub.fr.p$PFOS = sub.fr.p$LBXNFOS + sub.fr.p$LBXMFOS
# wtd_quantile_PAH <-wtd.quantile(sub.fr.p$LBXMPAH, weights=sub.fr.p$WTSB2YR,
#                                 probs=c(0, 0.4, 0.7, 0.9, 1), normwt=FALSE)
# wtd_quantile_FUA<-wtd.quantile(sub.fr.p$LBXPFUA, weights=sub.fr.p$WTSB2YR,
#                                probs=c(0, 0.4, 0.7, 0.9, 1), normwt=FALSE)
# wtd_quantile_PFDE<-wtd.quantile(sub.fr.p$LBXPFDE, weights=sub.fr.p$WTSB2YR,
#                                 probs=c(0, 0.1, 0.4, 0.8, 0.9, 1), normwt=FALSE)
# wtd_quantile_PFHS<-wtd.quantile(sub.fr.p$LBXPFHS, weights=sub.fr.p$WTSB2YR,
#                                 probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
# wtd_quantile_PFNA<-wtd.quantile(sub.fr.p$LBXPFNA, weights=sub.fr.p$WTSB2YR,
#                                 probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
# wtd_quantile_PFOA<-wtd.quantile(sub.fr.p$PFOA, weights=sub.fr.p$WTSB2YR,
#                                 probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
# wtd_quantile_PFOS<-wtd.quantile(sub.fr.p$PFOS, weights=sub.fr.p$WTSB2YR,
#                                 probs=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), normwt=FALSE)
# wtd_quantile_PAH[length(wtd_quantile_PAH)]<-Inf
# wtd_quantile_FUA[length(wtd_quantile_FUA)]<-Inf
# wtd_quantile_PFDE[length(wtd_quantile_PFDE)]<-Inf
# wtd_quantile_PFHS[length(wtd_quantile_PFHS)]<-Inf
# wtd_quantile_PFNA[length(wtd_quantile_PFNA)]<-Inf
# wtd_quantile_PFOA[length(wtd_quantile_PFOA)]<-Inf
# wtd_quantile_PFOS[length(wtd_quantile_PFOS)]<-Inf
# wtd_quantile_PAH[1]<--0.1
# wtd_quantile_FUA[1]<--0.1
# wtd_quantile_PFDE[1]<--0.1
# wtd_quantile_PFHS[1]<--0.1
# wtd_quantile_PFNA[1]<--0.1
# wtd_quantile_PFOA[1]<--0.1
# wtd_quantile_PFOS[1]<--0.1
# sub.fr.p$PAH = cut(sub.fr.p$LBXMPAH, wtd_quantile_PAH,right=FALSE, labels=c(1:4))
# sub.fr.p$PFUA = cut(sub.fr.p$LBXPFUA, wtd_quantile_FUA,right=FALSE,
#                     labels=c(1:4))
# sub.fr.p$PFDE = cut(sub.fr.p$LBXPFDE, wtd_quantile_PFDE,right=FALSE, labels =
#                       c(1:5))
# sub.fr.p$PFHS = cut(sub.fr.p$LBXPFHS, wtd_quantile_PFHS,right=FALSE,
#                     labels=c(1:10))
# sub.fr.p$PFNA = cut(sub.fr.p$LBXPFNA, wtd_quantile_PFNA,right=FALSE,
#                     labels=c(1:9))
# sub.fr.p$dPFOA = cut(sub.fr.p$PFOA, wtd_quantile_PFOA,right=FALSE,
#                      labels=c(1:10))
# sub.fr.p$dPFOS = cut(sub.fr.p$PFOS, wtd_quantile_PFOS,right=FALSE,
#                                                          labels=c(1:10))
# newmat.ds = sub.fr.p[,c(23:29,1)]
# names(newmat.ds) = c("me-PFOSA-AcOH", "PFUA", "PFDeA", "PFHxS", "PFNA",
#                      "PFOA", "PFOS", "SEQN")
# irt.grm.ds = grm(newmat.ds[,1:7])
# eap.ds <- ltm::factor.scores(irt.grm.ds, method="EAP", resp.patterns =
#                                newmat.ds[,1:7])$score.dat$z1
# summary(eap.ds)
# score.ds = cbind(newmat.ds, eap.ds)
# score.ds = score.ds[,c("SEQN", "eap.ds")]
# 
# 
# # Decile calculation summing isomers ----
# ####################################
# # Similar steps are repeated for the decile calculation when summing isomers
# fr.data=nhanes('PFAS_J')
# fr.data=fr.data[, c(1, 2, seq(4, 20, by=2))]
# 
# # Summed concentrations
# ####################################
# # Calculate the sum of concentrations for specified compounds
# sum.conc = apply(sub.fr.p[,c(12:17, 19:20)], 1, sum)
# sub.fr.p = cbind(sub.fr.p, sum.conc)
# head(sub.fr.p)
# score.sum = sub.fr.p[,c("SEQN", "sum.conc")]
# 
# # Combine all scores into one matrix and calculate correlations
# score.mat2 = merge(score.d, score.sum, by = "SEQN")
# score.mat3 = merge(score.mat2, score.ds, by = "SEQN")
# cor(score.mat3[, 2:4], method = "pearson")

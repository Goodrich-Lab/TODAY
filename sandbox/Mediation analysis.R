## Meet in the middle analysis
(analysis_pfas_names <- c(pfas_names_all, "score_pfsas", "score_pfcas")) 

overlapping_feature <- mim_res2$feature_name

X1 <-  data_scaled$pfas_pfna
Z1 <-  data_scaled %>% dplyr::select(all_of(covars)) %>% as.matrix()
M1 <- data_scaled %>% dplyr::select(all_of(overlapping_feature))
OT1 <- data_scaled$daystomic
status1 <- data_scaled$mic

# result1 <- HIMA::survHIMA(X = X1,
#                          Z = Z1,
#                          M = M1,
#                          OT = OT1,
#                          status = status1,
#                          scale = FALSE, 
#                          FDRcut = 0.2)
# 
# result_hima <- result %>% 
#   tidylog::left_join(mim_res2 %>% 
#                        dplyr::select(feature_name, EntrezGeneSymbol),
#                      by = c("ID"="feature_name")) %>% 
#   dplyr::select(ID, EntrezGeneSymbol, everything())

## Mediation analysis (independently for each mediator)
d1 <- data_scaled %>% 
  dplyr::select(pfas_pfna, overlapping_feature, covars, daystomic,mic) %>% 
  as.data.frame()

nexp<- 1#define number of exposure analytes in your dataset
nmed<- length(prot_names)#define number of mediators in your dataset
ncovars<- length(covars)#define number of covariates
noutcomes<- 1#define number of outcomes in your dataset
med.results<-as.data.frame(matrix(nrow=(nexp*nmed*noutcomes),ncol=41))
colnames(med.results)<-c('nobs', 'ACME.C','ACME.C.lo','ACME.C.hi','ACME.C.Pval','ACME.T','ACME.T.lo',
                         'ACME.T.hi','ACME.T.pval','ADE.C','ADE.C.lo','ADE.C.hi','ADE.C.Pval','ADE.T',
                         'ADE.T.lo','ADE.T.hi','ADE.T.pval','PMed.C','PMed.C.lo','PMed.C.hi','PMed.C.pval',
                         'PMed.T','PMed.T.lo','PMed.T.hi','PMed.T.pval','TE','TE.lo','TE.hi','TE.pval',
                         'ACME.avg','ACME.avg.lo','ACME.avg.hi','ACME.avg.pval','ADE.avg','ADE.avg.lo',
                         'ADE.avg.hi','ADE.avg.pval','PMed.avg','PMed.avg.lo','PMed.avg.hi','PMed.avg.pval')

#Loop to conduct pairwise mediation with multiple exposures and mediators
#Loop repeatedly subsets dataset d1 into d2 for individual pairs of exposures and biomarkers
k = 1
for(name in prot_names){
  
  d2<-data_scaled %>% 
    dplyr::select(pfas_pfna, name, all_of(covars), daystomic, mic) %>% as.data.frame() %>%
    rename(mediator = name)
  set.seed(111)

  med<- mediation::mediate(
    data = d2,
    model.m= lm(mediator~pfas_pfna+sex_male+agebase+eGFR+dxtime,
                  data=d2),
    model.y = survreg(Surv(daystomic,mic)~pfas_pfna+mediator+sex_male+agebase+eGFR+dxtime,
                      data=d2),
    treat="pfas_pfna",
    mediator = "mediator")
  model.y = survreg(Surv(daystomic,mic)~pfas_pfna+mediator+sex_male+agebase+eGFR+dxtime,
                    data=d2)
  # med <- CMAverse::cmest(data = d2, 
  #              model = "gformula",
  #              outcome = "mic",
  #              exposure = "pfas_pfna",
  #              mediator = "mediator", 
  #              event = "daystomic",
  #              basec = c("sex_male","agebase","eGFR","dxtime"),
  #              EMint = TRUE,
  #              mreg = list("linear"), 
  #              yreg = "coxph",
  #              astar = 0,
  #              a = 1,
  #              mval = list(1), 
  #              estimation = "imputation",
  #              inference = "bootstrap",
  #              nboot = 10)
  
  
  
  med.results[k,]<-cbind(nobs(model.y),med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p, med$d1, med$d1.ci[1],
                         med$d1.ci[2], med$d1.p, med$z0, med$z0.ci[1],med$z0.ci[2], med$z0.p, med$z1, 
                         med$z1.ci[1], med$z1.ci[2], med$z1.p, med$n0, med$n0.ci[1], med$n0.ci[2], 
                         med$n0.p, med$n1, med$n1.ci[1], med$n1.ci[2], med$n1.p, med$tau.coef, med$tau.ci[1], 
                         med$tau.ci[2], med$tau.p, med$d.avg, med$d.avg.ci[1], med$d.avg.ci[2], med$d.avg.p,
                         med$z.avg, med$z.avg.ci[1], med$z.avg.ci[2], med$z.avg.p, med$n.avg, med$n.avg.ci[1], 
                         med$n.avg.ci[2], med$n.avg.p)
        # rownames(med.results)[k]<-paste(colnames(d2)[1],name,colnames(d2)[2+ncovars+1],sep='_')
        rownames(med.results)[k] <- name
        print(k)
        print(rownames(med.results)[k])
        k=k+1
      }

## Adding name to the result
med_result <- med.results %>% rownames_to_column("AptName") %>%
  tidylog::left_join(meta_pro %>% dplyr::select(AptName, EntrezGeneSymbol)) %>%
  dplyr::select(AptName, EntrezGeneSymbol, everything())


write_csv(med_result,fs::path(dir_results,'pairwise_mediation_result.csv'))


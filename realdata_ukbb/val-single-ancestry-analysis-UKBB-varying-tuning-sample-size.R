rm(list=ls())
library(readr)
library(data.table)
library(mvtnorm)
library(devtools)
library(lavaan)
library(gdata)
library(xtable)
library(MASS) # for the ginv
library(data.table)
library(corpcor) #for pseudoinverse
library(parallel)
library(dplyr)
library(R.utils) # for gzip
library(stringr) # for str_detect
library(genio) # a package to facilitate reading and writing genetics data. The focus of this vignette is processing plink BED/BIM/FAM files.
library(data.table)
library(pROC)
library(bigsnpr)
#library(rms)
library(DescTools)

traits = c(paste0('rfMRI',1:76), paste0('abdominal',1:41), paste0('eye_oct',1:46), paste0('heart',1:82), 
           paste0('rfMRIe',1:1701), paste0('sMRI', 1:1437), paste0('dMRI', 1:675), paste0('tfMRI', 1:16),
           paste0('olinkmay',1:1463), paste0('olinknov',1:1458))
races = rep('EUR', 6995)

opt = list(LDrefpanel = '1kg', partitions = '0.8,0.2', delta = '0.001,0.01,0.1,1', nlambda = 30, 
           lambda.min.ratio = 0.01, alpha = '0.7,1.0,1.4', p_seq = '1.0e-05,3.2e-05,1.0e-04,3.2e-04,1.0e-03,3.2e-03,1.0e-02,3.2e-02,1.0e-01,3.2e-01,1.0e+00', 
           sparse = FALSE, kb = 500, Pvalthr = '5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01', R2 = '0.1', ensemble = T, type = "logical",
           verbose = 1)


# Example Manual Input (for testing the code): ---------
userID = 'jin'
submissionID = 'singleans'
methods = c('C+T', 'LDpred2')
LDrefpanel = '1kg' # '1kg' or 'ukbb', default: '1kg'
# Parameters for subsampling
k <- 2 # Default = 2 for 2-fold Monte Carlo Cross Validation (MCCV) to select tuning parameters 
MEs = c('CT', 'LDpred2')

# Optional input parameters:
partitions <- opt$partitions
homedir = '/path/to/PennPRS/Files/'
plink_path = '/path/to/PennPRS/software/'
type = 'single-ancestry' # 'multi-ancestry' or 'single-ancestry'
ld_path <- paste0(homedir) # set to the /LD folder
PUMAS_path = paste0(homedir,'code/')
threads = 1
ensemble = opt$ensemble

Indx = list()
methods = c('C+T', 'LDpred2')
for (m in 1:length(methods)){
  Indx[[m]] = data.frame(matrix(NA,length(traits),length(races)))
  rownames(Indx[[m]]) = traits; colnames(Indx[[m]]) = races
}
results = expand.grid(Method = c(methods, 'Ensemble', 'Ensemble2'), Race='EUR', Trait=traits, R2.ind = 0, R2.sum = 0, 
                      ct.p.ind = 0, ct.p.sum = 0, lassosum2.delta.ind = 0, lassosum2.lambda.ind = 0, 
                      lassosum2.delta.sum = 0, lassosum2.lambda.sum = 0, LDpred2.p.ind = 0, LDpred2.h2.ind = 0, 
                      LDpred2.sparse.ind = 0, LDpred2.p.sum = 0, LDpred2.h2.sum = 0, LDpred2.sparse.sum = 0)

updates = NULL

Nvals = c(100, 200, 400, 600, 800, 1000)
for (Nval in Nvals){
  for (i in c(1:245)){
    trait = traits[i]
    race = races[i]
    # Input: ---------
    # trait = 'HDL'
    # race = 'AFR'
    trait_name = paste0(race,'_',trait)
    ld_path0 <- paste0(ld_path, 'LD_1kg/') # set to the /LD_1kg folder under /LD/
    if (LDrefpanel == '1kg'){
      eval_ld_ref_path <- paste0(ld_path, '/1KGref_plinkfile/',race, '/') # set to the /1KGref_plinkfile folder under /LD/
      path_precalLD <- paste0(ld_path, '/LDpred2_lassosum2_corr_1kg/') # set to the /LDpred2_lassosum2_corr_1kg folder under /LD/
    } 
    # Job name/ID: e.g., trait_race_method_userID_submissionID
    jobID = paste(c(trait,race, paste0(methods,collapse = '.'), userID,submissionID), collapse = '_')
    # Create a job-specific (trait, race, methods, userID, jobID) directory to save all the outputs, set the working directory to this directory
    tempdir = '/path/to/PennPRS/'
    workdir = paste0(tempdir,jobID,'/')
    setwd(workdir) 
    
    
    # Parameters
    if ('C+T' %in% methods){
      # Parameters for the Clumping step
      kb = opt$kb # SNPs within 500kb of the index SNP are considered for clumping
      p.ldclump = 0.5 # P-value threshold for a SNP to be included as an index SNP
      Pvalthr = as.numeric(str_split(opt$Pvalthr,",")[[1]]) # default: c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01)
      R2 = as.numeric(str_split(opt$R2,",")[[1]])
      params.ct = expand.grid(r2 = R2, pvalthr = Pvalthr)
    }
    if ('lassosum2' %in% methods){
      delta = as.numeric(str_split(opt$delta,",")[[1]]) # candidate values of the shrinkage parameter in L2 regularization
      nlambda = opt$nlambda # number of different candidate values for lambda (shrinkage parameter in the L1 regularization). Default in lassosum2 pipeline: 30, which may lead to issues when using PUMAS subsampling to tune parameters
      lambda.min.ratio = opt$lambda.min.ratio # Ratio between the lowest and highest candidate values of lambda. Candidate values in (0,Inf), divided by comma
    }
    if ('LDpred2' %in% methods){
      h2.ratio = as.numeric(str_split(opt$alpha,",")[[1]])
      p_seq <- as.numeric(str_split(opt$p_seq,",")[[1]])  # Default
      sparse.option = opt$sparse # If TRUE: generate sparse effect size estimates. Options: subset of {FALSE, TRUE}
    }
    gwas_path <- paste0(workdir, 'sumdata/')
    output_path <- paste0(workdir, 'output/')
    input_path <- paste0(workdir, 'input_for_eval/')
    PennPRS_finalresults_path <- paste0(workdir, 'PennPRS_results/')
    eval_ld_ref = paste0(eval_ld_ref_path,LDrefpanel,'_hm3_',race,'_ref') # or hm3+mega
    # Create a separate directory 'PRS_model_training/' to store input for training PRS models
    prsdir0 = paste0(workdir, 'PRS_model_training/')
    output_path_eval = paste0(workdir, 'output_for_eval/')
    
    # -------- Calculate PRS:
    evaldir = '/path/to/PennPRS_backup/eval/'
    temdir = paste0(evaldir,trait)
    if (!dir.exists(temdir)){dir.create(temdir)}
    temdir = paste0(evaldir,trait,'/',race)
    if (!dir.exists(temdir)){dir.create(temdir)}
    temdir = paste0(evaldir,trait,'/',race,'/effect/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    temdir = paste0(evaldir,trait,'/',race,'/score/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    outdir = paste0(evaldir,trait,'/',race,'/')
    for (method in methods){
      prsdir = paste0(prsdir0, method,'/')
      outscorefile = paste0(outdir, 'score/',method,'.sscore')
      if (!file.exists(outscorefile)){
        prsf = paste0(prsdir, trait_name,'.',method,'.full.txt')
        if (file.exists(prsf)){
          sc = bigreadr::fread2(prsf) # only contains non-zero-effect SNPs
          prs.file = sc[,c(2,3,5:ncol(sc))]
          colnames(prs.file) = c('SNP', 'A1', paste0('BETA',1:(ncol(sc)-4)))
          #write.table(prs.file,file = paste0(outdir,'effect/',method,'.txt'), col.names = T,row.names = F,quote=F)
          write_delim(prs.file, paste0(outdir,'effect/',method,'.txt'), delim='\t')
          if (sum(sapply(1:3, function(x){grepl(c('MRI','abdominal','heart')[x], trait, fixed = TRUE)})) > 0) bfile = paste0('/path/to/PennPRS/bfile/test_subjects/ukb_imp_allchr_v3_40k_fmri_test_subjects')
          if (sum(sapply(1, function(x){grepl(c('eye')[x], trait, fixed = TRUE)})) > 0) bfile = paste0('/path/to/PennPRS/bfile/test_subjects/ukb_imp_allchr_v3_80k_eye_jan2022_test_subjects')
          if (sum(sapply(1, function(x){grepl(c('olink')[x], trait, fixed = TRUE)})) > 0) bfile = paste0('/path/to/PennPRS/bfile/test_subjects/ukb_imp_allchr_v3_50k_protein_may2023_test_subjects')
          prscode = paste(paste0(plink_path, 'plink2'),
                          paste0('--score ',  outdir,'effect/',method,'.txt'),
                          'cols=+scoresums,-scoreavgs',
                          paste0('--score-col-nums 3-',ncol(prs.file)),
                          # paste0(''),
                          paste0(' --bfile ', bfile),
                          " --threads 1",
                          paste0(' --out ', outdir, 'score/',method))
          system(prscode)
          print(paste('Complete calculating PRS for',trait,race,method))
        }
      }
    }
    
    method = 'LDpred2'
    prsdir = paste0(prsdir0, method,'/')
    outscorefile = paste0(outdir, 'score/',method,'.sscore')
    if (file.exists(outscorefile)){
      if ((grepl(c('rfMRI'), trait, fixed = TRUE)) & (!grepl(c('rfMRIe'), trait, fixed = TRUE))){
        validatetable = read.csv(paste0('/path/to/PennPRS/pheno/rfMRI/node/rfMRI_2_node_5mad_scale_resid_white_test_Jun20_2024.csv'))
        trait.indx = as.numeric(substr(trait,6,12))
      } 
      if (grepl(c('abdominal'), trait, fixed = TRUE)){
        validatetable = read.csv(paste0('/path/to/PennPRS/pheno/abdominal/ukb_abdominal_idp_41_042923_all_5mad_scale_resid_white_test_Jun20_2024.csv'))
        trait.indx = as.numeric(substr(trait,10,12))
      } 
      if (grepl(c('eye_oct'), trait, fixed = TRUE)){
        validatetable = read.csv(paste0('/path/to/PennPRS/pheno/eye/oct/eye_oct_and_others_april25_2022_v2_QC_5mad_scale_resid_white_test_Jun20_2024.csv'))
        trait.indx = as.numeric(substr(trait,8,10))
      } 
      validatetable = validatetable[,c(1,2+trait.indx)]
      if (grepl(c('olink'), trait, fixed = TRUE)) validatetable = validatetable[,c(1,1+trait.indx)]
      colnames(validatetable) = c('id','y')
      validatetable = validatetable[complete.cases(validatetable$y),]
      validatetable$id = as.character(validatetable$id)
      n.tuning = numeric() # c(8,120,33)
      for (m in 1:length(methods)){
        method = methods[m]
        prsfile = paste0(outdir, 'score/',method,'.sscore')
        tem = bigreadr::fread2(prsfile)
        tem = tem[,c(1,5:ncol(tem))]
        n.tuning[m] = ncol(tem)-1
        colnames(tem) = c('id',paste0('prs_',MEs[m],'_',1:n.tuning[m]))
        if (m == 1){
          preds = tem;
        } 
        if (m > 1){
          preds = merge(preds, tem); 
        } 
      }
      rownames(preds) = preds$id
      
      validatetable = merge(validatetable, preds, by='id')
      rownames(validatetable) = as.character(validatetable$id)
      #---------------------------------------#---------------------------------------
      set.seed(2024)
      ids1 = sample(rownames(validatetable), floor(nrow(validatetable)/2), replace = F)
      ids2 = setdiff(rownames(validatetable), ids1)
      ids1 = sample(ids1, size = Nval, replace = F)
      traindat = validatetable[validatetable$id %in% ids1,]
      # ----------------------- Train weighted sum -----------------------
      
      # Results from single method tuning based on individual data:
      output = list()
      for (m in 1:length(methods)){
        method = methods[m]
        output[[m]] = matrix(0,n.tuning[m],3)
        colnames(output[[m]]) = c('R2 Adjusted','Regression Coef','P-value')
        
        for(i in 1:n.tuning[m]){
          prstable = traindat[,c('y',paste0('prs_',MEs[m],'_',i))]
          colnames(prstable) = c('y','prs')
          prstable = prstable[complete.cases(prstable),]
          if ((nrow(prstable)>0)&(sum(prstable$prs) != 0)){
            # get residual:
            fit = lm(y~prs, data=prstable)
            output[[m]][i,'Regression Coef'] = coefficients(fit)['prs']
            output[[m]][i,'R2 Adjusted'] = summary(fit)$r.squared
            output[[m]][i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
            #print(i)
          }
        }
        if (method == 'C+T'){
          Indx[[m]][trait,race] = which.max(output[[m]][,'R2 Adjusted'])
        }  
        if (method == 'lassosum2'){
          Indx[[m]][trait,race] = which.max(output[[m]][,'R2 Adjusted'])
        } 
        if (method == 'LDpred2'){
          Indx[[m]][trait,race] = which.max(output[[m]][,'R2 Adjusted'])
        } 
      }
      
      set.seed(1)
      formula.wprs = formula(paste(paste('y', paste(c(sapply(1:length(methods), function(x){paste0('prs_',MEs[x],'_',Indx[[x]][trait,race])})),collapse="+"), sep='~')))
      s0 = lm(formula.wprs, data=traindat)
      
      # Individual level data results:
      # Validation
      valdat = validatetable[validatetable$id %in% ids2,]
      valdat$prs = predict(s0, valdat)
      fit = lm(y~prs, data=valdat)
      results[(results$Method == 'Ensemble') & (results$Race == race) & (results$Trait == trait), 'R2.ind'] = results[(results$Method == 'Ensemble2') & (results$Race == race) & (results$Trait == trait), 'R2.ind'] = summary(fit)$r.squared
      
      # Single best PRS:
      results[(results$Method == 'C+T') & (results$Race == race) & (results$Trait == trait), 'R2.ind'] = cor(valdat$y, valdat[,paste0('prs_',MEs[1],'_',Indx[[1]][trait,race])])^2
      results[(results$Method == 'lassosum2') & (results$Race == race) & (results$Trait == trait), 'R2.ind'] = cor(valdat$y, valdat[,paste0('prs_',MEs[2],'_',Indx[[2]][trait,race])])^2
      results[(results$Method == 'LDpred2') & (results$Race == race) & (results$Trait == trait), 'R2.ind'] = cor(valdat$y, valdat[,paste0('prs_',MEs[2],'_',Indx[[2]][trait,race])])^2
      
      # PUMAS data results:
      optim.indx = numeric()
      for (m in 1:length(methods)){
        method = methods[m]
        prsdir = paste0(prsdir0, method,'/')
        optimf = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData')
        if (file.exists(optimf)){
          load(optimf)
          results[(results$Method == method) & (results$Race == race) & (results$Trait == trait), 'R2.sum'] = cor(valdat$y, valdat[,paste0('prs_',MEs[m],'_',optimal.indx)])^2
          optim.indx[m] = optimal.indx
        }
      }
      pumas.w = bigreadr::fread2(paste0(PennPRS_finalresults_path, trait_name, '.', paste0(methods, collapse = '.'), '.omnibus.weights.txt'))
      pumas.w = colMeans(pumas.w)
      nonzero.indx = which(pumas.w!=0)
      if (length(nonzero.indx) > 0){
        valdat$prspumas = as.matrix(valdat[, sapply(c(1:length(methods))[nonzero.indx], function(x){paste0('prs_',MEs[x],'_',optim.indx[x])})]) %*% matrix(as.numeric(pumas.w[nonzero.indx]),ncol = 1)
        results[(results$Method == 'Ensemble') & (results$Race == race) & (results$Trait == trait), 'R2.sum'] = cor(valdat$y, valdat$prspumas)^2
      }
      
      pumas.w = bigreadr::fread2(paste0(PennPRS_finalresults_path, trait_name, '.', paste0(methods, collapse = '.'), '.omnibus.weights.alternative.txt'))
      pumas.w = colMeans(pumas.w)
      nonzero.indx = which(pumas.w!=0)
      if (length(nonzero.indx) > 0){
        valdat$prspumas = as.matrix(valdat[, sapply(c(1:length(methods))[nonzero.indx], function(x){paste0('prs_',MEs[x],'_',optim.indx[x])})]) %*% matrix(as.numeric(pumas.w[nonzero.indx]),ncol = 1)
        results[(results$Method == 'Ensemble2') & (results$Race == race) & (results$Trait == trait), 'R2.sum'] = cor(valdat$y, valdat$prspumas)^2
      }
      print(results[(results$Race == race) & (results$Trait == trait), c('Method', 'Race', 'Trait','R2.ind', 'R2.sum')])
      
      save(updates, trait, results, Indx, file=paste0("/path/to/PennPRS/Files/results/R2-UKBB-otherMRI-Nval=", Nval, ".RData"))
    }
  }
}




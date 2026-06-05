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

# traits = c(paste0('rfMRI',1:76), paste0('abdominal',1:41), paste0('eye_oct',1:46), paste0('heart',1:82), 
#          paste0('rfMRIe',1:1701), paste0('sMRI', 1:1437), paste0('dMRI', 1:675), paste0('tfMRI', 1:16),
#          paste0('olinkmay',1:1463), paste0('olinknov',1:1458))
# races = rep('EUR', 6995)
#rfMRI: 1:76
#abdominal: 77:117
#eye_oct: 118:163
#heart: 164:245
#rfMRIe: 246:1946
#sMRI: 1947:3383
#dMRI: 3384:4058
#tfMRI: 4059:4074
#olinkmay: 4075:5537
#olinknov: 5538:6995

opt = list(LDrefpanel = '1kg', partitions = '0.8,0.2', delta = '0.001,0.01,0.1,1', nlambda = 30,
           lambda.min.ratio = 0.01, alpha = '0.7,1.0,1.4', p_seq = '1.0e-05,3.2e-05,1.0e-04,3.2e-04,1.0e-03,3.2e-03,1.0e-02,3.2e-02,1.0e-01,3.2e-01,1.0e+00',
           sparse = FALSE, kb = 500, Pvalthr = '5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01', R2 = '0.1', ensemble = TRUE, type = "logical",
           verbose = 1)


# Example Manual Input (for testing the code): ---------
userID = 'user1'
submissionID = 'single_ans'
methods = c('C+T', 'lassosum2', 'LDpred2'); MEs = c('CT', 'lassosum2', 'LDpred2')
LDrefpanel = '1kg' # '1kg' or 'ukbb', default: '1kg'
# Parameters for subsampling
k <- 2 # Default = 2 for 2-fold Monte Carlo Cross Validation (MCCV) to select tuning parameters 

# Optional input parameters:
partitions <- opt$partitions
homedir = '/path/to/PennPRS/Files/'
plink_path = '/path/to/PennPRS/software/'
type = 'single-ancestry' # 'multi-ancestry' or 'single-ancestry'
ld_path <- paste0(homedir) # set to the /LD folder
PUMAS_path = paste0(homedir,'code/')
threads = 1
ensemble = opt$ensemble

trait = '1'
race = 'EUR'
trait_name = paste0(race,'_',trait)
ld_path0 <- paste0(ld_path, 'LD_1kg/') # set to the /LD_1kg folder under /LD/
jobID = paste(c(trait,race, paste0(methods,collapse = '.'), userID,submissionID), collapse = '_')
tempdir = '/path/to/PennPRS/summary/'
workdir = paste0(tempdir,jobID,'/')
suppressWarnings(dir.create(workdir))
setwd(workdir) 

Indx = list()
for (m in 1:length(methods)){
  Indx[[m]] = data.frame(matrix(NA,length(trait),length(race)))
  rownames(Indx[[m]]) = trait; colnames(Indx[[m]]) = race
}
results = expand.grid(Method = c(methods, 'Ensemble', 'Ensemble2'), Race='EUR', Trait=trait, R2.ind = 0, R2.sum = 0)
# results = expand.grid(Method = c(methods), Race='EUR', Trait=trait, R2.ind = 0, R2.sum = 0)



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
prsdir0 = paste0(workdir, 'PRS_model_training/')
output_path_eval = paste0(workdir, 'output_for_eval/')

# -------- Calculate PRS:
evaldir = '/path/to/PennPRS/eval/'
temdir = paste0(evaldir,trait)
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0(evaldir,trait,'/',race)
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0(evaldir,trait,'/',race,'/effect/')
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0(evaldir,trait,'/',race,'/score/')
if (!dir.exists(temdir)){dir.create(temdir)}
outdir = paste0(evaldir,trait,'/',race,'/')
# Calculate PRS for validation individuals using PLINK:
for (method in methods){
  prsdir = paste0(prsdir0, method,'/')
  outscorefile = paste0(outdir, 'score/',method,'.sscore')
  # if (!file.exists(outscorefile)){
    prsf = paste0(prsdir, trait_name,'.',method,'.full.txt')
    # if (file.exists(prsf)){
      sc = bigreadr::fread2(prsf)
      prs.file = sc[,c(2,3,5:ncol(sc))]
      colnames(prs.file) = c('SNP', 'A1', paste0('BETA',1:(ncol(sc)-4)))
      write_delim(prs.file, paste0(outdir,'effect/',method,'.txt'), delim='\t')
      bfile = paste0('/path/to/data/etc/projects/aging/gcta/imaging/non_eye_grm/ukb_imp_allchr_v3_80k_phase6_may2024')
      prscode = paste(paste0(plink_path, 'plink2'),
                      paste0('--score ',  outdir,'effect/',method,'.txt'),
                      'cols=+scoresums,-scoreavgs',
                      paste0('--score-col-nums 3-',ncol(prs.file)),
                      paste0(' --bfile ', bfile),
                      " --threads 1",
                      paste0(' --out ', outdir, 'score/',method))
      system(prscode)
      print(paste('Complete calculating PRS for',trait,race,method))
    # }
  # }
}


validatetable = read.csv(paste0('/path/to/PennPRS/individual/updated_DXA_traits_July2024_fixed.csv'), check.names = FALSE)
# Dynamically select the trait column by name and the ID column ('eid')
validatetable <- validatetable[, c('eid', trait)]
# Rename columns for compatibility with downstream processing
colnames(validatetable) <- c('id', 'y')
# Remove rows with missing values in the 'y' (trait) column
validatetable <- validatetable[complete.cases(validatetable$y), ]
# Ensure 'id' is treated as character data
validatetable$id <- as.character(validatetable$id)

# Initialize tuning parameter tracker (if necessary)
n.tuning <- numeric() # Example: c(8, 120, 33) for methods
for (m in 1:length(methods)){
  method = methods[m]
  prsfile = paste0(outdir, 'score/',method,'.sscore')
  tem = bigreadr::fread2(prsfile)
  tem = tem[,c(1,5:ncol(tem))]
  n.tuning[m] = ncol(tem)-1
  colnames(tem) = c('id',paste0('prs_',MEs[m],'_',1:n.tuning[m]))
  if (m == 1) preds = tem;
  if (m > 1) preds = merge(preds, tem); 
}
rownames(preds) = preds$id

validatetable = merge(validatetable, preds, by='id')
rownames(validatetable) = as.character(validatetable$id)
#---------------------------------------#---------------------------------------
set.seed(2024)
# randomly select half of the validation individuals to tune parameters for individual-level data-based PRS training
ids1 = sample(rownames(validatetable), floor(nrow(validatetable)/2), replace = F)
ids2 = setdiff(rownames(validatetable), ids1)
traindat = validatetable[validatetable$id %in% ids1,]
# ----------------------- Train weighted sum -----------------------

# Single method tuning based on individual data:
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
  Indx[[m]][trait,race] = which.max(output[[m]][,'R2 Adjusted'])
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
for (method in methods){
  method.indx = which(methods == method)
  results[(results$Method == method) & (results$Race == race) & (results$Trait == trait), 'R2.ind'] = cor(valdat$y, valdat[,paste0('prs_',MEs[method.indx],'_',Indx[[method.indx]][trait,race])])^2
}

# PUMAS results:
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
ensemble.weight.file = paste0(PennPRS_finalresults_path, trait_name, '.', paste0(methods, collapse = '.'), '.omnibus.weights.txt')
results[(results$Method == 'Ensemble') & (results$Race == race) & (results$Trait == trait), 'R2.sum'] = NA_real_
if (file.exists(ensemble.weight.file)){
  pumas.w = bigreadr::fread2(ensemble.weight.file)
  pumas.w = colMeans(pumas.w)
  nonzero.indx = which((pumas.w != 0) & !is.na(optim.indx))
  if (length(nonzero.indx) > 0){
    valdat$prspumas = as.matrix(valdat[, sapply(c(1:length(methods))[nonzero.indx], function(x){paste0('prs_',MEs[x],'_',optim.indx[x])})]) %*% matrix(as.numeric(pumas.w[nonzero.indx]),ncol = 1)
    results[(results$Method == 'Ensemble') & (results$Race == race) & (results$Trait == trait), 'R2.sum'] = cor(valdat$y, valdat$prspumas)^2
  }
} else {
  print(paste0('Ensemble weight file not found at ', ensemble.weight.file, '. Skipping Ensemble R2.sum. Run single-ancestry-step2-DXA.R if you need ensemble validation results.'))
}

ensemble.weight.file.alt = paste0(PennPRS_finalresults_path, trait_name, '.', paste0(methods, collapse = '.'), '.omnibus.weights.alternative.txt')
results[(results$Method == 'Ensemble2') & (results$Race == race) & (results$Trait == trait), 'R2.sum'] = NA_real_
if (file.exists(ensemble.weight.file.alt)){
  pumas.w = bigreadr::fread2(ensemble.weight.file.alt)
  pumas.w = colMeans(pumas.w)
  nonzero.indx = which((pumas.w != 0) & !is.na(optim.indx))
  if (length(nonzero.indx) > 0){
    valdat$prspumas = as.matrix(valdat[, sapply(c(1:length(methods))[nonzero.indx], function(x){paste0('prs_',MEs[x],'_',optim.indx[x])})]) %*% matrix(as.numeric(pumas.w[nonzero.indx]),ncol = 1)
    results[(results$Method == 'Ensemble2') & (results$Race == race) & (results$Trait == trait), 'R2.sum'] = cor(valdat$y, valdat$prspumas)^2
  }
} else {
  print(paste0('Alternative ensemble weight file not found at ', ensemble.weight.file.alt, '. Skipping Ensemble2 R2.sum. Run single-ancestry-step2-DXA.R if you need ensemble validation results.'))
}
print(results[(results$Race == race) & (results$Trait == trait),])

save(results, file=paste0("/path/to/PennPRS/Files/results/R2-UKBB-",trait,".RData"))

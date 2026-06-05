library(optparse)
library(parallel)
library(readr)
library(bigreadr)
library(bigsnpr)
library(data.table)
library(dplyr)
library(scales)
library(stringr) # for str_split

standardize_sumstats <- function(sumraw) {
  colnames_upper <- toupper(colnames(sumraw))

  find_col <- function(candidates, required = TRUE) {
    hit <- match(candidates, colnames_upper)
    hit <- hit[!is.na(hit)]
    if (length(hit) > 0) return(colnames(sumraw)[hit[1]])
    if (required) {
      stop(paste0(
        "Input summary statistics are missing required column(s): ",
        paste(candidates, collapse = "/"),
        ". Available columns are: ",
        paste(colnames(sumraw), collapse = ", ")
      ))
    }
    NA_character_
  }

  chr_col <- find_col(c("CHR", "#CHROM", "CHROM"))
  snp_col <- find_col(c("SNP", "RSID", "ID", "MARKERNAME"))
  a1_col <- find_col(c("A1", "EA", "EFFECT_ALLELE", "ALLELE1"))
  a2_col <- find_col(c("A2", "NEA", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "ALLELE2"))
  beta_col <- find_col(c("BETA", "EFFECT", "ESTIMATE"))
  se_col <- find_col(c("SE", "STDERR", "SEBETA"))
  p_col <- find_col(c("P", "PVAL", "PVALUE", "P_VALUE"))
  n_col <- find_col(c("N", "N_EFF", "NEFF", "TOTAL_N", "OBS_CT"))
  af_col <- find_col(
    c("AF1", "A1FREQ", "A1_AF", "EAF", "EAF1", "ALT_FREQ", "FRQ", "MAF"),
    required = FALSE
  )

  if (is.na(af_col)) {
    stop(paste0(
      "Input summary statistics must contain an allele-frequency column such as AF1 or MAF. ",
      "Available columns are: ", paste(colnames(sumraw), collapse = ", ")
    ))
  }

  standardized <- data.frame(
    CHR = sumraw[[chr_col]],
    SNP = sumraw[[snp_col]],
    A1 = sumraw[[a1_col]],
    A2 = sumraw[[a2_col]],
    BETA = sumraw[[beta_col]],
    SE = sumraw[[se_col]],
    P = sumraw[[p_col]],
    N = sumraw[[n_col]],
    MAF = sumraw[[af_col]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  standardized
}

# Add option parsing to input trait as a command-line argument
option_list = list(
  make_option(c("--trait"), type="character", default=NULL, help="Trait ID", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
option = parse_args(opt_parser)

opt = list(LDrefpanel = '1kg', partitions = '0.8,0.2', delta = '0.001,0.01,0.1,1', nlambda = 30,
           lambda.min.ratio = 0.01, alpha = '0.7,1.0,1.4', p_seq = '1.0e-05,3.2e-05,1.0e-04,3.2e-04,1.0e-03,3.2e-03,1.0e-02,3.2e-02,1.0e-01,3.2e-01,1.0e+00',
           sparse = FALSE, kb = 500, Pvalthr = '5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01', R2 = '0.1', ensemble = T, type = "logical",
           verbose = 1)


# Example Manual Input (default for testing): ---------
userID = 'user1'
submissionID = 'single_ans'
methods = c('C+T', 'lassosum2', 'LDpred2')
trait = option$trait # Set trait dynamically from command-line input
race = 'EUR'
LDrefpanel = '1kg'
# ----------------------------------------------------
# Parameters for subsampling
k <- 2 # Default = 2 for 2-fold Monte Carlo Cross Validation (MCCV) to select tuning parameters
ensemble = FALSE # consider training ensemble PRS combining PRS trained by the selected methods
# ----------------


if (length(methods) > 1) ensemble = TRUE

# Optional input parameters:
partitions <- opt$partitions
homedir = '/path/to/PennPRS/Files/'
plink_path = '/path/to/PennPRS/software/'
type = 'single-ancestry' # 'multi-ancestry' or 'single-ancestry'
ld_path <- paste0(homedir) # set to the /LD folder
PUMAS_path = paste0(homedir,'code/')
threads = 1



trait_name = paste0(race,'_',trait)
ld_path0 <- paste0(ld_path, 'LD_1kg/') # set to the /LD_1kg folder under /LD/
if (LDrefpanel == '1kg'){
  eval_ld_ref_path <- paste0(ld_path, '/1KGref_plinkfile/',race,"/") # set to the /1KGref_plinkfile folder under /LD/
  path_precalLD <- paste0(ld_path, '/LDpred2_lassosum2_corr_1kg/') # set to the /LDpred2_lassosum2_corr_1kg folder under /LD/
}
# Job name/ID: e.g., trait_race_method_userID_submissionID
jobID = paste(c(trait,race, paste0(methods,collapse = '.'), userID,submissionID), collapse = '_')
# Create a job-specific (trait, race, methods, userID, jobID) directory to save all the outputs, set the working directory to this directory
tempdir = '/path/to/PennPRS/summary/'
workdir = paste0(tempdir,jobID,'/')
suppressWarnings(dir.create(workdir))
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
  sp.temp = toupper(str_split(opt$sparse,",")[[1]])
  sparse.option = ifelse(sp.temp == 'TRUE', TRUE, FALSE) # opt$sparse # If TRUE: generate sparse effect size estimates. Options: subset of {FALSE, TRUE}
}

#source(paste0(PUMAS_path, 'PennPRS_functions.R')) # please save the PennPRS_functions.R file to the /PUMAS/code/ directory
gwas_path <- paste0(workdir, 'sumdata/')
output_path <- paste0(workdir, 'output/')
input_path <- paste0(workdir, 'input_for_eval/')
PennPRS_finalresults_path <- paste0(workdir, 'PennPRS_results/')
eval_ld_ref = paste0(eval_ld_ref_path,LDrefpanel,'_hm3_',race,'_ref') # or hm3+mega
dir.create(gwas_path, showWarnings = F)
dir.create(output_path, showWarnings = F)
dir.create(input_path, showWarnings = F)
dir.create(PennPRS_finalresults_path, showWarnings = F)
# Create a separate directory 'PRS_model_training/' to store input for training PRS models
prsdir0 = paste0(workdir, 'PRS_model_training/')
if (!dir.exists(prsdir0)) dir.create(prsdir0)
output_path_eval = paste0(workdir, 'output_for_eval/')
dir.create(output_path_eval, showWarnings = F)
# Create a separate directory 'PRS_model_training/' to store input for training PRS models
for (method in methods){
  prsdir = paste0(prsdir0, method,'/')
  if (!dir.exists(prsdir)) dir.create(prsdir)
}
if (ensemble) ensemble.methods = methods

# copy the input GWAS summary data, {Ancestry}_{Trait}.txt, to the /sumdata/ folder
# example code to move GWAS summary data:
system(paste0('cp ', '/path/to/PennPRS/summary/','/',trait_name,'.txt ', workdir, 'sumdata/'))

# # Ancestry: choose from EUR,AFR,AMR,EAS,SAS
# # for eye OCT data, you can use the following code to download and extract the summary data to the /sumdata/ folder:
# setwd('/path/to/PennPRS/inputfiles/')
# hm3snps = bigreadr::fread2('/path/to/PennPRS/Files/hapmap3rsid.txt')
# if (grepl(c('eye_oct'), trait, fixed = TRUE)){
#   i = as.numeric(substr(trait,8,10))
#   system(paste0('zcat /path/to/PennPRS/gwas/eye/oct/eye_oct_80k_march10_2022_pheno',i,'.fastGWA.gz > /path/to/PennPRS/inputfiles/eye_oct',i,'.fastGWA'))
#   a=bigreadr::fread2(paste0('eye_oct',i,'.fastGWA'))
#   a = a[a$SNP %in% hm3snps$SNP,c('CHR', 'SNP', 'A1', 'A2', 'AF1', 'BETA', 'SE', 'P', 'N')]
#   colnames(a) = c('CHR', 'SNP', 'A1', 'A2', 'MAF', 'BETA', 'SE', 'P', 'N')
#   write_delim(a, paste0('EUR_eye_oct',i,'.txt'), delim = '\t')
#   print(i)
#   system(paste0('mv /path/to/PennPRS/inputfiles/', trait_name,'.txt ', workdir, 'sumdata/'))
#   system(paste0('rm -rf /path/to/PennPRS/inputfiles/eye_oct',i,'.fastGWA'))
#   rm(a)
# }
#############################################################################################



######## QC for GWAS Summary Data:
cat(paste0("\n********************************************"))
cat(paste0("\n**** Step 0: QC for the input GWAS data ****"))
cat(paste0("\n********************************************\n"))

sumraw = bigreadr::fread2(paste0(workdir, 'sumdata/', trait_name, '.txt'))
sumraw = standardize_sumstats(sumraw)
sumraw$BETA = as.numeric(sumraw$BETA)
sumraw$SE = as.numeric(sumraw$SE)
sumraw$MAF = as.numeric(sumraw$MAF)
sumraw$P = as.numeric(sumraw$P)
sumraw$N = as.numeric(sumraw$N)
# 0. Are there any SNP that have reasonable z-score?
chi2_thr = 30
remaining.SNPs = which(abs(sumraw$BETA/sumraw$SE) < sqrt(chi2_thr))
if (length(remaining.SNPs) < 5){
  stop(paste0("[Terminated] Job is terminated because less than 5 SNPs have z-score < sqrt(30), suggesting issues with the input GWAS data."))
}

n.na = sum(!complete.cases(sumraw))
if (n.na > 0){
  sumraw = sumraw[complete.cases(sumraw), ]
  if (n.na == 1) print(paste0('* 1 SNP has missing GWAS summary-level information and is removed.'))
  if (n.na > 1) print(paste0('* ', n.na, ' SNPs have missing GWAS summary-level information and are removed.'))
}


# 1. Remove SNPs with problematic BETA
beta.thr = 1e3
rm.indx1 = which(abs(sumraw$BETA) > beta.thr)
if (length(rm.indx1) > 0){
  if (length(rm.indx1) == 1) print(paste0('* 1 SNP has problematic GWAS summary statistic with abs(BETA) > ', beta.thr, ' and is removed.'))
  if (length(rm.indx1) > 1) print(paste0('* ', length(rm.indx1), ' SNPs have problematic GWAS summary statistics with abs(BETA) > ', beta.thr, ' and are removed.'))
} 

# 2. Remove SNPs with problematic p-values
rm.indx2 = which( ((sumraw$P) > 1) | (sumraw$P < 0))
if (length(rm.indx2) > 0){
  if (length(rm.indx2) == 1) print(paste0('* 1 SNP has p-value > 1 or < 0 and is removed.'))
  if (length(rm.indx2) > 1) print(paste0('* ', length(rm.indx2), ' SNPs have p-value > 1 or < 0 and are removed.'))
} 

# 3. Remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size.
rm.indx3 = numeric()
# N.90percentile = quantile(sumraw$N, 0.1)
# rm.indx3 = which(sumraw$N < N.90percentile)
# if (length(rm.indx3) > 0){
#   if (length(rm.indx3) == 1) print(paste0('* 1 SNP has an effective sample size less than 0.67 times the 90th percentile of the total sample size and is removed.'))
#   if (length(rm.indx3) > 1) print(paste0('* ', length(rm.indx3), ' SNPs have an effective sample size less than 0.67 times the 90th percentile of the total sample size and are removed.'))
# } 

# 4. Remove SNPs with extremely large effect sizes (z^2> 100) 
chi2.thr = 1e3
rm.indx4 = which((sumraw$BETA/sumraw$SE)^2 > chi2.thr)
if (length(rm.indx4) > 0){
  if (length(rm.indx4) == 1) print(paste0('* 1 SNP has an extremely large effect size  (z-score^2 > ', chi2.thr, ') and is removed.'))
  if (length(rm.indx4) > 1) print(paste0('* ', length(rm.indx4), ' SNPs have extremely large effect sizes  (z-score^2 > ', chi2.thr, ') and are removed.'))
} 

# 5. Remove SNPs with zero SE 
rm.indx5 = which(sumraw$SE == 0)
if (length(rm.indx5) > 0){
  if (length(rm.indx5) == 1) print(paste0('* 1 SNP has SE = 0 and is removed.'))
  if (length(rm.indx5) > 1) print(paste0('* ', length(rm.indx5), ' SNPs have SE = 0 and are removed.'))
} 
rm.indx = unique(c(rm.indx1, rm.indx2, rm.indx3, rm.indx4, rm.indx5))

if (length(rm.indx) > 0){
  sumraw = sumraw[-rm.indx, ]
  if (nrow(sumraw) == 0){
    stop(paste0("[Terminated] 0 SNPs remaining after QC. Job terminated.\n * Please check the quality of the input GWAS summary data and make sure the columns are in correct format."))
  }
  if (nrow(sumraw) > 0){
    write_delim(sumraw, path = paste0(workdir, 'sumdata/', trait_name, '.txt'), delim='\t')
    if (length(rm.indx) == 1) print(paste0('* 1 problematic SNP removed. QC step completed.'))
    if (length(rm.indx) > 1) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. ', length(rm.indx), ' problematic SNPs removed.'))
  }
}
if (length(rm.indx) == 0) print(paste0('* QC step completed. ', nrow(sumraw), ' SNPs remaining. No SNP was removed.'))


# --------------------------------------------------------------------
# --------------------- Step 1: PUMAS Subsampling ---------------------
# --------------------------------------------------------------------
# Different from PRS-CS (auto) which doesn't have tuning parameters, C+T, lassosum2, and LDpred2 have tuning parameters, 
# and thus we need to conduct k=4 fold MCCV.
# Note!!! For a relatively large tuning sample size (here in our case ~300K), then there is no need to do MCCV and k=2 is sufficient.
# This is usually the case for EUR
# But for non-EUR races, it's very common that the tuning sample size is only several thousands
# If tuning sample size is below 2000 we use k=4 MCCV?, o.w. we just use k=2
if ( opt$verbose >= 1 ) {
  cat(paste0("\n***********************************"))
  cat(paste0("\n**** Step 1: PUMAS Subsampling ****"))
  cat(paste0("\n***********************************\n"))
}


ld_file <- paste0(ld_path0,race,'_LD_hm3.RData')
rs_file <- paste0(ld_path0,race,'_rs_hm3.RData')
pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.subsampling.customized.R '),
                  paste0('--k ',k), 
                  paste0('--partitions ',partitions),
                  paste0('--trait_name ',trait_name),
                  paste0('--gwas_path ', gwas_path),
                  paste0('--ld_file ', ld_file),
                  paste0('--rs_file ', rs_file),
                  paste0('--output_path ', output_path),
                  paste0('--threads ', threads))
system(pumascode)



# --------------------------------------------------------------------------------------------
# --------------------- Step 2: Train PRS models using different methods ---------------------
# --------------------------------------------------------------------------------------------
if ( opt$verbose >= 1 ) {
  cat(paste0("\n**********************************************************"))
  cat(paste0("\n**** Step 2: Train PRS models using different methods ****"))
  cat(paste0("\n**********************************************************\n"))
}
if ('C+T' %in% methods){
  method = 'C+T'
  NCORES = 1
  prsdir = paste0(prsdir0, method,'/')
  
  # Submit k separate jobs (for ite in 1:k) to the server and run them in parallel. 
  for (ite in 1:k){
    # --------------------- Step 2.1: Preparation for input files for C+T ---------------------
    pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
    if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.R.'))
    if (file.exists(pumasout)){
      sumraw = bigreadr::fread2(pumasout)
      sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
      colnames(sumstats) <- c("CHR", "SNP", "REF", "ALT", "BETA", "SE", "P", "N", "FRQ") # MAF or REF_FRQ
      rownames(sumstats) = sumstats$SNP
      # REF: effect allele
      print(paste0('Complete Loading GWAS summary data for iteration ', ite))
    }
    
    temdir = paste0(prsdir,'snplist/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    temdir = paste0(prsdir,'clumped/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    temdir = paste0(prsdir,'filtered/')
    if (!dir.exists(temdir)){dir.create(temdir)}
    
    # Create base data (summary statistic) file containing the P-value information:
    sumstats_input <- sumstats[,c('SNP','P')] #,'REF_FRQ'
    fwrite(sumstats_input, file=paste0(prsdir,'snplist/',trait_name,'_ite',ite,'.txt'),row.names = F, quote = F, sep=' ')
    
    # --------------------- Step 2.2: Run C+T ----------------------------
    set.seed(2023)
    for (r2 in R2){
      # --------------------- LD Clumping (C) ---------------------
      ldclumpcode <- paste0(plink_path, 'plink --bfile ', eval_ld_ref,
                            ' --clump ',prsdir,'snplist/',trait_name,'_ite',ite,'.txt',
                            ' --clump-p1 ',p.ldclump,
                            # ' --clump-p2 ',pc,
                            ' --clump-r2 ',r2,
                            ' --clump-kb ',kb,
                            ' --threads 1',
                            ' --silent',
                            ' --out ',prsdir,'clumped/',trait_name,'_ite',ite,'_r2=',r2)
      system(ldclumpcode)
      # --------------------- Thresholding (T) ---------------------
      LD <- bigreadr::fread2(paste0(prsdir,'clumped/',trait_name,'_ite',ite,'_r2=',r2,'.clumped'))
      clumped.snp <- LD[,3,drop=F][,1]
      sumstats.clumped <- sumstats[clumped.snp,]
      
      for (pvalthr in Pvalthr){
        keep.SNP = sumstats.clumped[sumstats.clumped$P <= pvalthr,c('SNP')]
        sumstats[,paste0('r2_', r2,'_p_',pvalthr)] = sumstats$BETA
        dump.SNP = which(!sumstats$SNP %in% keep.SNP)
        sumstats[dump.SNP, paste0('r2_', r2,'_p_',pvalthr)] = 0
      }
    }
    SCORE = sumstats[,c('CHR','SNP','REF','ALT',  sapply(1:nrow(params.ct), function(x){paste0('r2_', params.ct[x,'r2'],'_p_',params.ct[x,'pvalthr'])}))]
    colnames(SCORE)[3:4] = c('A1','A2')
    write_delim(SCORE, path = paste0(prsdir, trait_name,'.',method,'.ite',ite,".txt"), delim='\t')
  }
  rm(SCORE)
}





if (('lassosum2' %in% methods) | ('LDpred2' %in% methods)){
  NCORES = 17
  map_ldref <- readRDS(paste0(ld_path, 'map_',LDrefpanel,'_ldref_',race,'.rds'))
  
  # Submit k separate jobs (for ite in 1:k) to the server and run them in parallel. 
  # Each job will require < 20G, the memory depends on NCORES (if NCORES = 3 then perhaps 10G is enough).
  for (ite in 1:k){
    # --------------------- Step 2.1: Preparation for input files for each PRS method (here: lassosum2) ---------------------
    # Reformat summary data to use as the input data for lassosum2:
    pumasout = paste0(output_path, trait_name, '.gwas.ite', ite, '.txt')
    if (!file.exists(pumasout)) print(paste0('Subsampling failed to generate ', ite,'-th fold of the MCCV summary data. Rerun pumas.subsampling.R.'))
    if (file.exists(pumasout)){
      sumraw = bigreadr::fread2(pumasout)
      sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
      names(sumstats) <- c("chr", "rsid", "a1", "a0", "beta", "beta_se", "p", "n_eff", "a1_sumdata_af")
      # a0: effect allele
      
      info_snp <- snp_match(sumstats, map_ldref, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
      info_snp <- tidyr:: drop_na(tibble::as_tibble(info_snp))
      sd_ldref <- with(info_snp, sqrt(2 * a1_af * (1 - a1_af)))
      sd_ss <- with(info_snp, sqrt(2 * a1_sumdata_af * (1 - a1_sumdata_af)))
      is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
      df_beta <- info_snp[!is_bad, ]
      print(paste0('Complete pre-processing GWAS summary data for iteration ', ite))
    }
    
    if (ite == 1){
      td = paste0(prsdir0, 'temporary_LDpred2_lassosum2_ite',ite)
      if (!dir.exists(td)) dir.create(td)
      setwd(td)
      tmp <- tempfile(tmpdir = td)
      
      for (chr in 1:22) {
        cat(chr, ".. ", sep = "")
        ## indices in 'df_beta'
        ind.chr <- which(df_beta$chr == chr)
        ## indices in 'map_ldref'
        ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
        ## indices in 'corr0'
        ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
        if (length(ind.chr3) > 0){
          # corr0
          corr0 <- readRDS(paste0(path_precalLD, race, '/LD_ref_chr', chr, '.rds'))[ind.chr3, ind.chr3]
          if (chr == 1) {
            ld <- Matrix::colSums(corr0^2)
            corr <- as_SFBM(corr0, tmp, compact = TRUE)
          } else {
            if (length(corr0) == 1) corr0 =  as(1, "sparseMatrix") # as(corr0, "sparseMatrix") # as.matrix(corr0, 1, 1)
            ld <- c(ld, Matrix::colSums(corr0^2))
            corr$add_columns(corr0, nrow(corr))
          }
          print(paste0('Complete calculating LD for CHR ', chr))
          rm(corr0)
        }
      }
    }
    
    if ('lassosum2' %in% methods){
      method = 'lassosum2'
      prsdir = paste0(prsdir0, method,'/')
      # Parameters
      delta = as.numeric(str_split(opt$delta,",")[[1]]) # candidate values of the shrinkage parameter in L2 regularization
      
      # --------------------- Step 2.2: Run lassosum2 ----------------------------
      set.seed(2023)
      beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = NCORES, 
                                      delta = delta, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
      if (ite == 1){
        params.lassosum2 <- attr(beta_lassosum2, "grid_param")
        write_delim(params.lassosum2, path = paste0(prsdir, 'params.lassosum2.txt'), delim='\t')
      }
      
      beta_lassosum2[is.na(beta_lassosum2)] = 0
      # Further fix potential non-convergent issues:
      n.nonconvergent = sapply(1:nrow(params.lassosum2), function(x){sum(abs(beta_lassosum2[,x])>1)})
      indx.nonconvergent = which(n.nonconvergent > 0)
      if (length(indx.nonconvergent)>0) beta_lassosum2[,indx.nonconvergent] = 0
      
      beta_lassosum2 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_lassosum2)
      colnames(beta_lassosum2) = c('chr','rsid','a1','a0', paste0('lassosum2_',1:nrow(params.lassosum2)))
      output_lassosum2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
      write_delim(beta_lassosum2,path = output_lassosum2, delim='\t')
      
      print(paste0('** Complete training ', method, ' for MCCV ite ', ite, ' **'))
    }
    
    
    if ('LDpred2' %in% methods){
      method = 'LDpred2'
      prsdir = paste0(prsdir0, method,'/')
      
      (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                      sample_size = n_eff, blocks = NULL)))
      ldsc_h2_est <- abs(ldsc[["h2"]])
      h2_seq <- round(ldsc_h2_est * h2.ratio, 5); 
      h2_seq[h2_seq == 0] = 1e-5
      h2_seq[duplicated(h2_seq)] = h2_seq[duplicated(h2_seq)] * 1.01
      n.inflated = sum(h2_seq>1)
      if (n.inflated > 0) h2_seq[h2_seq>1] = 0.95 + seq(0,0.01*(n.inflated-1), by = 0.01)
      params.ldpred2 <- expand.grid(p = p_seq, h2 = h2_seq, sparse = sparse.option)
      
      set.seed(2023)
      beta_ldpred2 <- snp_ldpred2_grid(corr, df_beta, params.ldpred2, ncores = NCORES)
      beta_ldpred2[is.na(beta_ldpred2)] = 0
      # Further fix potential non-convergent issues:
      # n.nonconvergent = sapply(1:nrow(params.ldpred2), function(x){sum(abs(beta_ldpred2[,x])>1)})
      # indx.nonconvergent = which(n.nonconvergent > 0)
      # if (length(indx.nonconvergent)>0) beta_ldpred2[,indx.nonconvergent] = 0
      
      beta_ldpred2 = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_ldpred2)
      colnames(beta_ldpred2) = c('chr','rsid','a1','a0', paste0('LDpred2_',1:nrow(params.ldpred2)))
      output_LDpred2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
      write_delim(beta_ldpred2,path = output_LDpred2, delim='\t')
      
      print(paste0('** Complete training ', method, ' for MCCV ite ', ite, ' **'))
    }
  }
}



# ---------------------------------------------------------------------------------------
# ----------------------------- Step 3: Single Model Tuning -----------------------------
# ---------------------------------------------------------------------------------------
# single_prs() in PUMA-CUBS.evaluation.R
if ( opt$verbose >= 1 ) {
  cat(paste0("\n********************************************************"))
  cat(paste0("\n******* Step 3: Parameter Tuning for Each Method *******"))
  cat(paste0("\n********************************************************\n"))
}


if ('C+T' %in% methods){
  method = 'C+T'
  prsdir = paste0(prsdir0, method,'/')
  
  # This step is also needed for C+T, these input files for evaluation will be used for calculating R2 on testing data and for training ensemble PRS
  for (ite in 1:k){
    output_ct = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
    if(file.exists(output_ct)){
      score = bigreadr::fread2(output_ct) 
      n.tuning = ncol(score) - 4
      colnames(score) = c('CHR', 'SNP', 'A1', 'A2', paste0('BETA',1:n.tuning))
    }
    
    # Match alleles with GWAS summary data:
    sumstats = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
    stateval = sumstats[, c('SNP','A1','A2')]
    colnames(stateval) = c('SNP', 'A1.ref','A2.ref')
    stateval = left_join(stateval,score,by="SNP")
    
    na.ind = which(is.na(stateval$A1))
    if (length(na.ind) > 0){
      stateval[na.ind, paste0('BETA',1:n.tuning)] = 0
      stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
      stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
    }
    flipped = which(stateval$A1.ref != stateval$A1)
    print(paste0(length(flipped), ' flipped SNPs.'))
    if (length(flipped) > 0){
      stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
      stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
      for (t in 1:n.tuning){
        stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
      }
    }
    scores = stateval[,c('CHR','SNP','A1','A2',paste0('BETA',1:n.tuning))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
    write_delim(scores,path = paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim='\t')
    # write.table(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
    rm(stateval)
  }
  
  xty_path = stats_path = output_path # the "output_path" used for storing output from pumas.subsampling.R
  R2.tuned = cbind(params.ct, matrix(0,length(R2)*length(Pvalthr),k))
  colnames(R2.tuned) = c('r2','pvalthr',paste0('ite',1:k))
  
  tunecode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                   paste0('--k ',k), 
                   paste0('--ref_path ', eval_ld_ref),
                   paste0('--trait_name ',trait_name),
                   paste0('--prs_method ', method),
                   paste0('--xty_path ', xty_path),
                   paste0('--stats_path ', stats_path),
                   paste0('--weight_path ', prsdir),
                   paste0('--output_path ', output_path_eval))
  system(tunecode)
  ##### Extract parameter tuning results
  R2.tuned = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  # Select tuning parameters
  nonzero.ind = which(colMeans(R2.tuned) != 0)
  params.tuned.ct = nonzero.ind[which.max(colMeans(R2.tuned)[nonzero.ind])]
  r2 = params.ct[params.tuned.ct,'r2']; pvalthr = params.ct[params.tuned.ct,'pvalthr']
  r2.ct = r2; pval.ct = pvalthr
  print(paste0('Tuned parameters: r2 = ', r2, ', pval = ', pvalthr))
  tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  save(r2, pvalthr, file = tuned.parameters.file)
}


# Additional filtering for lassosum2 and LDpred2 tuning parameter settings
if ('lassosum2' %in% methods){
  method = 'lassosum2'
  prsdir = paste0(prsdir0, method,'/')
  
  for (ite in 1:k){
    output_lassosum2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
    if(file.exists(output_lassosum2)){
      score = bigreadr::fread2(output_lassosum2) 
      n.tuning = ncol(score) - 4 # as.numeric(strsplit(colnames(score)[ncol(score)],split='_')[[1]][2])
      colnames(score) = c('CHR', 'SNP', 'A1', 'A2', paste0('BETA',1:n.tuning))
    }
    
    # Match alleles with GWAS summary data:
    sumstats = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
    stateval = sumstats[, c('SNP','A1','A2')]
    colnames(stateval) = c('SNP', 'A1.ref','A2.ref')
    stateval = left_join(stateval,score,by="SNP")
    
    na.ind = which(is.na(stateval$A1))
    if (length(na.ind) > 0){
      stateval[na.ind, paste0('BETA',1:n.tuning)] = 0
      stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
      stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
    }
    flipped = which(stateval$A1.ref != stateval$A1)
    print(paste0(length(flipped), ' flipped SNPs.'))
    if (length(flipped) > 0){
      stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
      stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
      for (t in 1:n.tuning){
        stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
      }
    }
    scores = stateval[,c('CHR','SNP','A1','A2',paste0('BETA',1:n.tuning))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
    write_delim(scores,path = paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim='\t')
    # write.table(scores, paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
    rm(stateval)
  }
  
  xty_path = stats_path = output_path # the "output_path" used forstoring output from pumas.subsampling.R
  
  pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                    paste0('--k ',k), 
                    paste0('--ref_path ', eval_ld_ref),
                    paste0('--trait_name ',trait_name),
                    paste0('--prs_method ',method),
                    paste0('--xty_path ', xty_path),
                    paste0('--stats_path ', stats_path),
                    paste0('--weight_path ', input_path),
                    paste0('--output_path ', output_path_eval))
  system(pumascode)
  
  
  ##### Extract parameter tuning results
  r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  r2.avg = colMeans(r2)
  r2.order = list()
  n.candidates = 20
  for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
  # Select the top parameter settings:
  # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
  params.tuned = Reduce(intersect, r2.order)
  nonzero.indx = which(r2.avg != 0)
  params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
  # print(paste0('Maximum r2 of ',method,': ', max(r2.avg)))
  tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  save(params.tuned, file = tuned.parameters.file)
  if (length(params.tuned) == 0){
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    nonzero.indx = which(r2.avg != 0)
    r2.order = list()
    n.candidates = length(nonzero.indx)
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk, nonzero.indx]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    params.tuned = nonzero.indx[params.tuned[!is.na(params.tuned)]]
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    save(params.tuned, file = tuned.parameters.file)
  }
}


if ('LDpred2' %in% methods){
  method = 'LDpred2'
  prsdir = paste0(prsdir0, method,'/')
  
  for (ite in 1:k){
    output_LDpred2 = paste0(prsdir, trait_name,'.',method,'.ite',ite,'.txt')
    if(file.exists(output_LDpred2)){
      score = bigreadr::fread2(output_LDpred2) 
      n.tuning = ncol(score) - 4 # as.numeric(strsplit(colnames(score)[ncol(score)],split='_')[[1]][2])
      colnames(score) = c('CHR', 'SNP', 'A1', 'A2', paste0('BETA',1:n.tuning))
    }
    
    # Match alleles with GWAS summary data:
    stateval = bigreadr::fread2(paste0(output_path, trait_name, '.gwas.ite', ite, '.txt'))
    stateval = stateval[, c('SNP','A1','A2')]
    colnames(stateval) = c('SNP', 'A1.ref','A2.ref')
    stateval = left_join(stateval,score,by="SNP")
    
    na.ind = which(is.na(stateval$A1))
    if (length(na.ind) > 0){
      stateval[na.ind, paste0('BETA',1:n.tuning)] = 0
      stateval[na.ind,'A1'] = stateval[na.ind,'A1.ref']
      stateval[na.ind,'A2'] = stateval[na.ind,'A2.ref']
    }
    flipped = which(stateval$A1.ref != stateval$A1) # April 14: 0 (already corrected)
    print(paste0(length(flipped), ' flipped SNPs.'))
    if (length(flipped) > 0){
      stateval[flipped,'A1'] = stateval[flipped,'A1.ref']
      stateval[flipped,'A2'] = stateval[flipped,'A2.ref']
      for (t in 1:n.tuning){
        stateval[flipped,paste0('BETA',t)] = - stateval[flipped,paste0('BETA',t)]
      }
    }
    scores = stateval[,c('CHR','SNP','A1','A2',paste0('BETA',1:n.tuning))] # other files: SNP	CHR	A1	BETA1	BETA2	A2
    write_delim(scores,path = paste0(input_path, trait_name,'.',method, '.ite',ite,'.txt'), delim='\t')
    # write.table(scores, paste0(input_path, trait_name,'.',method, '.ite', ite, '.txt'), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
    rm(stateval, scores)
  }
  
  xty_path = stats_path = output_path # the "output_path" used forstoring output from pumas.subsampling.R
  
  pumascode = paste(paste0('Rscript ', PUMAS_path, 'PUMAS.evaluation.customized.R'),
                    paste0('--k ',k), 
                    paste0('--ref_path ', eval_ld_ref),
                    paste0('--trait_name ',trait_name),
                    paste0('--prs_method ',method),
                    paste0('--xty_path ', xty_path),
                    paste0('--stats_path ', stats_path),
                    paste0('--weight_path ', input_path),
                    paste0('--output_path ', output_path_eval))
  system(pumascode)
  
  ##### Extract parameter tuning results
  r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  r2.avg = colMeans(r2)
  r2.order = list()
  n.candidates = 20
  for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk,]),decreasing = T)[1:n.candidates]
  # Select the top parameter settings:
  # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
  params.tuned = Reduce(intersect, r2.order)
  nonzero.indx = which(r2.avg != 0)
  params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
  # print(paste0('Maximum r2 of ',method,': ', max(r2.avg)))
  tuned.parameters.file = paste0(workdir,'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  save(params.tuned, file = tuned.parameters.file)
  if (length(params.tuned) == 0){
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    nonzero.indx = which(r2.avg != 0)
    r2.order = list()
    n.candidates = length(nonzero.indx)
    for (kk in 1:k) r2.order[[kk]] = order(as.numeric(r2[kk, nonzero.indx]),decreasing = T)[1:n.candidates]
    # Select the top parameter settings:
    # params.tuned = as.numeric(substr(names(sort(r2.avg, decreasing = T)[1:5]), 5, 10))
    params.tuned = Reduce(intersect, r2.order)
    params.tuned = nonzero.indx[params.tuned[!is.na(params.tuned)]]
    params.tuned = unique(params.tuned[(params.tuned <= length(r2.avg)) & (params.tuned %in% nonzero.indx)])
    save(params.tuned, file = tuned.parameters.file)
  }
}




# ---------------------------------------------------------------------------------------
# --------------------- Step 4: Generate best PRS(s) for each method --------------------
# ---------------------------------------------------------------------------------------
# Once we select the "best" PRSs for different models, train on the whole dataset
# For LDpred2 and lassosum, select the best ones that are not all 0 or NA's
# Then come back and train the best PRS model
if ( opt$verbose >= 1 ) {
  cat(paste0("\n************************************************************"))
  cat(paste0("\n******* Step 4: Generate best PRS(s) for each method *******"))
  cat(paste0("\n************************************************************\n"))
}
# ---------------------------------------------------------------------------------------
# ---------- Step 4.1: Train the best PRS(s) from each method on the whole data ---------
# ---------------------------------------------------------------------------------------

sumraw = bigreadr::fread2(paste0(output_path,trait_name,".gwas_matched.txt"))
err.CT = err.lassosum2 = err.LDpred2 = 0

# C+T: find the one best PRS 
if ('C+T' %in% methods){
  method = 'C+T'
  prsdir = paste0(prsdir0, method,'/')
  NCORES = 1
  if ( opt$verbose >= 1 ){
    print(paste0('************************************************************'))
    print(paste0('******* Start training ', method, ' on the original GWAS data *******'))
    print(paste0('************************************************************'))
  }
  # The default number of tuning parameter settings is 12, so NCORES can be 3 (default), 6, or 12. 
  # Parameters for the Clumping step
  z = qnorm(p=p.ldclump/2,lower.tail=FALSE)
  Pvalthr = as.numeric(str_split(opt$Pvalthr,",")[[1]])
  R2 = as.numeric(str_split(opt$R2,",")[[1]])
  params.ct = expand.grid(r2 = R2, pvalthr = Pvalthr)
  
  # --------------------- Step 2.1: Preparation for input files for C+T ---------------------
  sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
  sumstats$P = as.numeric(sumstats$P)
  sumstats$BETA = as.numeric(sumstats$BETA)
  sumstats$SE = as.numeric(sumstats$SE)
  sumstats$N = as.numeric(sumstats$N)
  sumstats$MAF = as.numeric(sumstats$MAF)
  colnames(sumstats) <- c("CHR", "SNP", "REF", "ALT", "BETA", "SE", "P", "N", "FRQ") # MAF or REF_FRQ
  rownames(sumstats) = sumstats$SNP
  sumstats0 = sumstats
  # REF: effect allele
  print(paste0('** Complete Loading GWAS summary data **'))
  
  tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
  load(tuned.parameters.file) # Load tuned r2 and pvalthr
  
  
  # Create base data (summary statistic) file containing the P-value information:
  sumstats_input <- sumstats[,c('SNP','P')] #,'REF_FRQ'
  fwrite(sumstats_input, file=paste0(prsdir,'snplist/',trait_name,'.txt'),row.names = F, quote = F, sep=' ')

  # Save all candidate C+T models trained on the full GWAS so downstream
  # validation can compare individual-level tuning against summary-based tuning.
  candidate.cols.ct = sapply(1:nrow(params.ct), function(x){
    paste0('r2_', params.ct[x, 'r2'], '_p_', params.ct[x, 'pvalthr'])
  })
  sumstats.full = sumstats0
  for (candidate.col in candidate.cols.ct) sumstats.full[[candidate.col]] = 0
  for (r2.full in R2){
    ldclumpcode.full <- paste0(plink_path, 'plink --bfile ', eval_ld_ref,
                               ' --clump ',prsdir,'snplist/',trait_name,'.txt',
                               ' --clump-p1 ',p.ldclump,
                               ' --clump-r2 ',r2.full,
                               ' --clump-kb ',kb,
                               ' --threads 1',
                               ' --silent',
                               ' --out ',prsdir,'clumped/',trait_name,'_r2=',r2.full)
    system(ldclumpcode.full)
    clump.file.full = paste0(prsdir,'clumped/',trait_name,'_r2=',r2.full,'.clumped')
    clumped.snp.full = character(0)
    if (file.exists(clump.file.full)){
      LD.full = bigreadr::fread2(clump.file.full)
      if (nrow(LD.full) > 0) clumped.snp.full = LD.full[,3,drop=F][,1]
      rm(LD.full)
    }
    sumstats.clumped.full = sumstats0[sumstats0$SNP %in% clumped.snp.full, ]
    for (pvalthr.full in Pvalthr){
      candidate.col = paste0('r2_', r2.full, '_p_', pvalthr.full)
      keep.SNP.full = sumstats.clumped.full[sumstats.clumped.full$P <= pvalthr.full, 'SNP']
      if (length(keep.SNP.full) > 0){
        keep.ind.full = sumstats.full$SNP %in% keep.SNP.full
        sumstats.full[keep.ind.full, candidate.col] = sumstats0[keep.ind.full, 'BETA']
      }
    }
  }
  ct.full = sumstats.full[, c('CHR','SNP','REF','ALT', candidate.cols.ct)]
  colnames(ct.full)[3:4] = c('A1','A2')
  nonzero.rows.ct = which(rowSums(abs(as.matrix(ct.full[, candidate.cols.ct, drop = FALSE])) > 0) > 0)
  if (length(nonzero.rows.ct) > 0) ct.full = ct.full[nonzero.rows.ct, ]
  write_delim(ct.full, path = paste0(prsdir, trait_name,'.',method,'.full.txt'), delim='\t')
  
  # --------------------- Step 2.2: Run C+T ----------------------------
  set.seed(2023)
  # --------------------- LD Clumping (C) ---------------------
  ldclumpcode <- paste0(plink_path, 'plink --bfile ', eval_ld_ref,
                        ' --clump ',prsdir,'snplist/',trait_name,'.txt',
                        ' --clump-p1 ',p.ldclump,
                        # ' --clump-p2 ',pc,
                        ' --clump-r2 ',r2,
                        ' --clump-kb ',kb,
                        ' --threads 1',
                        ' --silent',
                        ' --out ',prsdir,'clumped/',trait_name,'_r2=',r2)
  system(ldclumpcode)
  # --------------------- Thresholding (T) ---------------------
  LD <- bigreadr::fread2(paste0(prsdir,'clumped/',trait_name,'_r2=',r2,'.clumped'))
  clumped.snp <- LD[,3,drop=F][,1]
  sumstats.clumped <- sumstats[clumped.snp,]
  
  keep.SNP = sumstats.clumped[sumstats.clumped$P <= pvalthr,c('SNP')]
  dump.SNP = which(!sumstats$SNP %in% keep.SNP)
  sumstats[dump.SNP, 'BETA'] = 0
  
  SCORE = sumstats[,c('CHR','SNP','REF','ALT', 'BETA')]
  colnames(SCORE)[3:4] = c('A1','A2')
  # Remove SNPs with zero effects:
  beta_ct.out = data.frame(CHR = ct.full$CHR,
                           SNP = ct.full$SNP,
                           A1 = ct.full$A1,
                           A2 = ct.full$A2,
                           BETA = ct.full[, candidate.cols.ct[which((params.ct$r2 == r2) & (params.ct$pvalthr == pvalthr))]])
  beta_ct.out = beta_ct.out[beta_ct.out$BETA != 0,]
  
  ### If all effects are 0:
  R2.tuned = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
  nonzero.ind = which(colMeans(R2.tuned) != 0)
  params.tuned.ct = nonzero.ind[which.max(colMeans(R2.tuned)[nonzero.ind])]
  while((nrow(beta_ct.out) == 0) & (length(nonzero.ind) > 0)){
    nonzero.ind = nonzero.ind[-which(nonzero.ind == params.tuned.ct)]
    if (length(nonzero.ind) == 0) break
    params.tuned.ct = nonzero.ind[which.max(colMeans(R2.tuned)[nonzero.ind])]
    r2 = params.ct[params.tuned.ct,'r2']; pvalthr = params.ct[params.tuned.ct,'pvalthr']
    beta_ct.out = data.frame(CHR = ct.full$CHR,
                             SNP = ct.full$SNP,
                             A1 = ct.full$A1,
                             A2 = ct.full$A2,
                             BETA = ct.full[, candidate.cols.ct[params.tuned.ct]])
    beta_ct.out = beta_ct.out[beta_ct.out$BETA != 0,]
  }
  if (nrow(beta_ct.out) > 0){
    prs_ct_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
    write_delim(beta_ct.out, path = prs_ct_outputfile, delim='\t')
    
    ct.optimal.indx = which((params.ct$r2 == r2) & (params.ct$pvalthr == pvalthr))
    optimal.indx = ct.optimal.indx
    save(optimal.indx, file = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData'))
    rm(optimal.indx)
    write_delim(params.ct[ct.optimal.indx,], path = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'))
    if ( opt$verbose >= 1 ){
      print(paste0('*************************************'))
      print(paste0('******* Complete training ', method, ' *******'))
      print(paste0('*************************************'))
    }
    # If R2<0:
    R2.tuned = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    if (colMeans(R2.tuned[ct.optimal.indx]) <= 0){
      err.CT = 1
      if (ensemble){
        ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
        print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power and thus we will not incorporate it in the ensemble PRS. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
      } 
      if (!ensemble){
        print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
      }
    }
  }
  if (nrow(beta_ct.out) == 0){
    err.CT = 2
    if (ensemble){
      ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
      print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs and thus we will not incorporate it in the ensemble PRS. Please try other methods.'))
    } 
    if (!ensemble){
      print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. Please try other methods.'))
    }
  }
  rm(SCORE)
}




if (('lassosum2' %in% methods) | ('LDpred2' %in% methods)){
  NCORES = 17
  map_ldref <- readRDS(paste0(ld_path, 'map_',LDrefpanel,'_ldref_',race,'.rds'))
  
  sumstats = sumraw[,c('CHR','SNP','A1','A2','BETA','SE','P','N', 'MAF')]
  sumstats$P = as.numeric(sumstats$P)
  sumstats$BETA = as.numeric(sumstats$BETA)
  sumstats$SE = as.numeric(sumstats$SE)
  sumstats$N = as.numeric(sumstats$N)
  sumstats$MAF = as.numeric(sumstats$MAF)
  names(sumstats) <- c("chr", "rsid", "a1", "a0", "beta", "beta_se", "p", "n_eff", "a1_sumdata_af")
  
  info_snp <- snp_match(sumstats, map_ldref, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
  info_snp <- tidyr:: drop_na(tibble::as_tibble(info_snp))
  sd_ldref <- with(info_snp, sqrt(2 * a1_af * (1 - a1_af)))
  sd_ss <- with(info_snp, sqrt(2 * a1_sumdata_af * (1 - a1_sumdata_af)))
  is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
  df_beta <- info_snp[!is_bad, ]
  

  
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  ldsc_h2_est <- abs(ldsc[["h2"]])
  cat(paste0('Heritability estimate based on LD score regression: ', signif(ldsc[["h2"]], 3)))
  if (ldsc[["h2"]] < 0) cat(paste0('Warning: negative hertability estimate based on LD score regression.'))
  
  if ('lassosum2' %in% methods){
    method = 'lassosum2'
    prsdir = paste0(prsdir0, method,'/')
    if ( opt$verbose >= 1 ){
      print(paste0('******************************************************************'))
      print(paste0('******* Start training lassosum2 on the original GWAS data *******'))
      print(paste0('******************************************************************'))
    }
    delta = as.numeric(str_split(opt$delta,",")[[1]]) # candidate values of the shrinkage parameter in L2 regularization
    params.lassosum2 = bigreadr::fread2(paste0(prsdir, 'params.lassosum2.txt'))
    
    # Load tuned parameters (top choices from summary-data tuning) and train the
    # full candidate grid on the whole GWAS for downstream individual-level tuning.
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    tuned.parameters.file = paste0(workdir, 'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
    load(tuned.parameters.file)
    
    set.seed(2023)
    beta_lassosum2.full = snp_lassosum2(corr, df_beta, ncores = NCORES,
                                        delta = delta, nlambda = nlambda,
                                        lambda.min.ratio = lambda.min.ratio)
    params.lassosum2.full = attr(beta_lassosum2.full, "grid_param")
    lassosum2.ref.keys = paste(params.lassosum2$lambda, signif(params.lassosum2$delta, 15), sep = '|')
    lassosum2.full.keys = paste(params.lassosum2.full$lambda, signif(params.lassosum2.full$delta, 15), sep = '|')
    lassosum2.reorder = match(lassosum2.ref.keys, lassosum2.full.keys)
    if (any(is.na(lassosum2.reorder))){
      stop('Unable to align the lassosum2 full-data parameter grid with the original candidate ordering.')
    }
    beta_lassosum2.full = beta_lassosum2.full[, lassosum2.reorder, drop = FALSE]
    beta_lassosum2.full[is.na(beta_lassosum2.full)] = 0
    n.nonconvergent = sapply(seq_len(ncol(beta_lassosum2.full)), function(x){sum(abs(beta_lassosum2.full[,x]) > 1)})
    indx.nonconvergent = which(n.nonconvergent > 0)
    if (length(indx.nonconvergent) > 0) beta_lassosum2.full[, indx.nonconvergent] = 0
    
    lassosum2.full = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_lassosum2.full)
    colnames(lassosum2.full) = c('CHR','SNP','A1','A2', paste0('lassosum2_', 1:nrow(params.lassosum2)))
    nonzero.rows.lassosum2 = which(rowSums(abs(as.matrix(lassosum2.full[, -(1:4), drop = FALSE])) > 0) > 0)
    if (length(nonzero.rows.lassosum2) > 0) lassosum2.full = lassosum2.full[nonzero.rows.lassosum2, ]
    write_delim(lassosum2.full, path = paste0(prsdir, trait_name,'.',method,'.full.txt'), delim='\t')
    
    if (length(params.tuned) == 0){
      err.lassosum2 = 2
      if (ensemble){
        ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
        print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. The full candidate matrix is still saved for downstream individual-level validation, but no summary-based optimal model was selected.'))
      }
      if (!ensemble){
        print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. The full candidate matrix is still saved for downstream individual-level validation, but no summary-based optimal model was selected.'))
      }
    }
    if (length(params.tuned) > 0){
      nonzero.indx = which(sapply(params.tuned, function(x){sum(abs(beta_lassosum2.full[,x]) > 1e-7) > 0}))
      if (length(nonzero.indx) > 0){
        candidates = params.tuned[nonzero.indx]
        if (length(candidates) == 1){
          optimal.indx = candidates[1]
        }
        if (length(candidates) > 1){
          stop = 0
          ii = 0
          while ((stop == 0) & (ii < length(candidates))){
            ii = ii + 1
            if (((candidates[ii] - 1) %in% candidates) | ((candidates[ii] + 1) %in% candidates)){
              stop = 1
            }
          }
          if (stop == 0) ii = 1
          optimal.indx = candidates[ii]
        }
        save(optimal.indx, file = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData'))
        lassosum2.out = signif(params.lassosum2[optimal.indx, , drop = FALSE], 4)
        write_delim(lassosum2.out,
                    path = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'))
        beta_lassosum2 = data.frame(df_beta[,c('chr','rsid','a1','a0')],
                                    BETA = beta_lassosum2.full[, optimal.indx])
        colnames(beta_lassosum2) = c('CHR','SNP','A1','A2', 'BETA')
        
        nonzero = (sum(abs(beta_lassosum2$BETA) > 1e-7) > 0)
        if (nonzero){
          beta_lassosum2.out = beta_lassosum2[beta_lassosum2$BETA != 0,]
          prs_lassosum2_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
          write_delim(beta_lassosum2.out, path = prs_lassosum2_outputfile, delim='\t')
          print(paste0('Optimal parameter setting: '))
          print(paste0('Delta: ', signif(params.lassosum2[optimal.indx,'delta'], digits = 3)))
          print(paste0('Lambda: ', signif(params.lassosum2[optimal.indx,'lambda'], digits = 3)))
          print(paste0('Complete training ', method))
          
          if (r2.avg[optimal.indx] <= 0){
            err.lassosum2 = 1
            if (ensemble){
              ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power and thus we will not incorporate it in the ensemble PRS. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            }
            if (!ensemble){
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            }
          }
        }
        if (!nonzero){
          err.lassosum2 = 2
          if (ensemble){
            ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
          }
          if (!ensemble){
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
          }
        }
      }
      if (length(nonzero.indx) == 0){
        err.lassosum2 = 2
        if (ensemble){
          ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
        }
        if (!ensemble){
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
        }
      }
      if ( opt$verbose >= 1 ){
        print(paste0('*************************************'))
        print(paste0('******* Complete training ', method, ' *******'))
        print(paste0('*************************************'))
      }
      suppressWarnings(rm(optimal.indx))
    }
  }
  
  
  if ('LDpred2' %in% methods){
    method = 'LDpred2'
    prsdir = paste0(prsdir0, method,'/')
    
    r2 = bigreadr::fread2(paste0(output_path_eval,trait_name, '.', method, '.txt'))
    r2.avg = colMeans(r2)
    tuned.parameters.file = paste0(workdir,'PRS_model_training/',method,'/tuned_parameters_',trait_name,'.RData')
    load(tuned.parameters.file)
    
    h2_seq <- round(ldsc_h2_est * h2.ratio, 5)
    h2_seq[h2_seq == 0] = 1e-5
    h2_seq[duplicated(h2_seq)] = h2_seq[duplicated(h2_seq)] * 1.01
    n.inflated = sum(h2_seq > 1)
    if (n.inflated > 0) h2_seq[h2_seq > 1] = 0.95 + seq(0,0.01 * (n.inflated - 1), by = 0.01)
    params.ldpred2.train <- expand.grid(p = p_seq, h2 = h2_seq, sparse = sparse.option)
    
    set.seed(2023)
    beta_ldpred2.full = snp_ldpred2_grid(corr, df_beta, params.ldpred2.train, ncores = NCORES)
    beta_ldpred2.full[is.na(beta_ldpred2.full)] = 0
    ldpred2.full = data.frame(df_beta[,c('chr','rsid','a1','a0')], beta_ldpred2.full)
    colnames(ldpred2.full) = c('CHR','SNP','A1','A2', paste0('LDpred2_', 1:nrow(params.ldpred2.train)))
    nonzero.rows.ldpred2 = which(rowSums(abs(as.matrix(ldpred2.full[, -(1:4), drop = FALSE])) > 0) > 0)
    if (length(nonzero.rows.ldpred2) > 0) ldpred2.full = ldpred2.full[nonzero.rows.ldpred2, ]
    write_delim(ldpred2.full, path = paste0(prsdir, trait_name,'.',method,'.full.txt'), delim='\t')
    
    if (length(params.tuned) == 0){
      err.LDpred2 = 2
      if (ensemble){
        ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
        print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. The full candidate matrix is still saved for downstream individual-level validation, but no summary-based optimal model was selected.'))
      }
      if (!ensemble){
        print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs. The full candidate matrix is still saved for downstream individual-level validation, but no summary-based optimal model was selected.'))
      }
    }
    if (length(params.tuned) > 0){
      nonzero.indx = which(sapply(params.tuned, function(x){sum(abs(beta_ldpred2.full[,x]) > 1e-10) > 0}))
      if (length(nonzero.indx) > 0){
        candidates = params.tuned[nonzero.indx]
        if (length(candidates) == 1){
          optimal.indx = candidates[1]
        }
        if (length(candidates) > 1){
          stop = 0
          ii = 0
          while ((stop == 0) & (ii < length(candidates))){
            ii = ii + 1
            if (((candidates[ii]-1) %in% candidates) | ((candidates[ii]+1) %in% candidates)){
              stop = 1
            }
          }
          if (stop == 0) ii = 1
          optimal.indx = candidates[ii]
        }
        save(optimal.indx, file = paste0(prsdir, trait_name, '.', method, '.optimal.indx.RData'))
        ldpred2.out = params.ldpred2.train[optimal.indx, , drop = FALSE]
        write_delim(ldpred2.out,
                    path = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.optimal_params.txt'))
        beta_ldpred2 = data.frame(df_beta[,c('chr','rsid','a1','a0')],
                                  BETA = beta_ldpred2.full[, optimal.indx])
        colnames(beta_ldpred2) = c('CHR','SNP','A1','A2', 'BETA')
        
        nonzero = (sum(abs(beta_ldpred2[,'BETA']) > 1e-10) > 0)
        if (nonzero){
          beta_ldpred2.out = beta_ldpred2[beta_ldpred2$BETA != 0,]
          prs_ldpred2_outputfile = paste0(PennPRS_finalresults_path, trait_name,'.',method,'.PRS.txt')
          write_delim(beta_ldpred2.out, path = prs_ldpred2_outputfile, delim='\t')
          print(paste0('Optimal parameter setting: '))
          print(paste0('p: ', signif(params.ldpred2.train[optimal.indx,'p'], digits = 3)))
          print(paste0('h2: ', signif(params.ldpred2.train[optimal.indx,'h2'], digits = 3)))
          print(paste0('sparse: ', signif(params.ldpred2.train[optimal.indx,'sparse'], digits = 3)))
          print(paste0('Complete training ', method))
          
          if (r2.avg[optimal.indx] <= 0){
            err.LDpred2 = 1
            if (ensemble){
              ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power and thus we will not incorporate it in the ensemble PRS. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            }
            if (!ensemble){
              print(paste0('[Warning] Trained PRS model based on ', method, ' has an estimated R2 < 0, which indicates that the PRS lacks prediction power. \nPotential explanations for R2 < 0 are:\n 1. The trait is not heritable.\n 2. The GWAS have insufficient power (e.g., due to low sample size) to develop a predictive PRS.\n 3. Issues with the input GWAS summary data (e.g., problematic BETA or SE).\n 4. ', method, ' is not powerful for developing PRS for the trait, in which case other methods can be considered.'))
            }
          }
        }
        if (!nonzero){
          err.LDpred2 = 2
          if (ensemble){
            ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
          }
          if (!ensemble){
            print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
          }
        }
      }
      if (length(nonzero.indx) == 0){
        err.LDpred2 = 2
        if (ensemble){
          ensemble.methods = ensemble.methods[which(ensemble.methods != method)]
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues) and thus we will not incorporate it in the ensemble PRS. Please consider other methods.'))
        }
        if (!ensemble){
          print(paste0('Trained PRS model based on ', method, ' has zero effect estimate for all SNPs (potentially due to convergence issues). Please consider other methods.'))
        }
      }
      if (opt$verbose >= 1 ){
        print(paste0('*************************************'))
        print(paste0('******* Complete training ', method, ' *******'))
        print(paste0('*************************************'))
      }
    }
  }
  rm(corr, ld)
  # system(paste0('rm -rf ', paste0(prsdir0, 'temporary_LDpred2_lassosum2_ite1')))
}



if (ensemble) save(ensemble.methods, err.CT, err.lassosum2, err.LDpred2, file = paste0(PennPRS_finalresults_path, 'step1.RData'))
if (!ensemble) save(err.CT, err.lassosum2, err.LDpred2, file = paste0(PennPRS_finalresults_path, 'step1.RData'))

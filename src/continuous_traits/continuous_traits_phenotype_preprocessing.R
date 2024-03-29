config <- yaml::read_yaml(here::here("config.yml"))

### biochemistry ###
biochem <- read.csv(paste0(config$baseline_phenotypes_path, "clinical/biochemistry.csv"))

### body ###
body <- read.csv(paste0(config$baseline_phenotypes_path, "clinical/body.csv"))

### Blood pressure and heart rate ###
bphr <- read.csv(paste0(config$baseline_phenotypes_path, "clinical/BPHR.csv"))

### Cognitive ###
cog <- read.csv(paste0(config$baseline_phenotypes_path, "clinical/cognitive.csv"))

### Lung ###
lung <- read.csv(paste0(config$baseline_phenotypes_path, "clinical/spirometry.csv"))

### Alcohol ###
alc <- read.csv(paste0(config$baseline_phenotypes_path, "PCQ/alcohol.csv"))

### Smoking ###
smk <- read.csv(paste0(conig$baseline_phenotypes_path, "pack_years.csv"))

### merge files ### 
tmp <- merge(biochem[,1:8], body[,c(1,2,4,7,8)], by.x="ID", by.y="id", all=T)
tmp1 <- merge(tmp, bphr[,c(1,8:10)], by.x="ID", by.y="id", all=T)
tmp2 <- merge(tmp1, cog[,c(1,4,5,9,10,11)], by.x="ID", by.y="ID", all=T)
tmp3 <- merge(tmp2, lung[,c(1,11,12)], by.x="ID", by.y="id", all=T)
tmp4 <- merge(tmp3, alc, by.x="ID", by.y="id", all=T)
tmp5 <- merge(tmp4, smk, by.x="ID", by.y="Sample_Name", all=T)

### load in target file and merge ###
tar <- readRDS(paste0(config$continuous_traits_data_path, "GS20k_Targets.rds"))

tmp6 <- merge(tmp5, tar[,c(1,3,4,5)], by.x="ID", by.y="Sample_Name")


##############################################
### recode 0 to NA for each cognitive test ###
##############################################

tmp6$LM <- tmp6$logical_mem_1 + tmp6$logical_mem_2
tmp6$LM[tmp6$LM==0] <- NA
tmp6$digit_symbol[tmp6$digit_symbol==0] <- NA
tmp6$vocabulary[tmp6$vocabulary==0] <- NA
tmp6$verbal_total[tmp6$verbal_total==0] <- NA

### outlier removal - lifted from https://github.com/davebraze/FDB1/blob/master/R/outliers.R ###

outlierID <- function(x, cut=4) {
    xx <- scale(x)
    retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
    retval
}

outlierTrim <- function(x, cut=4) {
    id <- outlierID(x, cut=cut)
    retval <- ifelse(id, NA, x)
}


waves <- list('w1' = tmp6[tmp6$Set == 'wave1',],
              'w3' = tmp6[tmp6$Set == 'wave3',],
              'w4' = tmp6[tmp6$Set == 'wave4',])

waves <- lapply(waves, function(w) {
    w$Set <- NULL
    w
})

w4PostFamilyFilteringIDs <- readRDS(paste0(continuous_traits_data_path, "w4_ids_after_family_filtering.rds"))

waves$w4 <- waves$w4[waves$w4$ID %in% w4PostFamilyFilteringIDs,]

process_wave <- function(tmp7, wave) {

    ### recode outliers to NA ###
    
    tmp7$Glucose_trim <- outlierTrim(log(tmp7$Glucose))
    tmp7$Chol_trim <- outlierTrim(tmp7$Total_cholesterol)
    tmp7$HDL_trim <- outlierTrim(tmp7$HDL_cholesterol)
    tmp7$Sodium_trim <- outlierTrim(tmp7$Sodium)
    tmp7$Potassium_trim <- outlierTrim(tmp7$Potassium)
    tmp7$Urea_trim <- outlierTrim(tmp7$Urea)
    tmp7$Creatinine_trim <- outlierTrim(tmp7$Creatinine)
    
    tmp7$bmi_trim <- ifelse(tmp7$bmi<18 | tmp7$bmi>50, NA, log(tmp7$bmi))
    tmp7$whr_trim <- outlierTrim(tmp7$whr)
    tmp7$fat_trim <- outlierTrim(tmp7$body_fat)
    tmp7$height_trim <- outlierTrim(tmp7$height)
    
    tmp7$sBP_trim <- outlierTrim(tmp7$avg_sys)
    tmp7$dBP_trim <- outlierTrim(tmp7$avg_dia)
    tmp7$HR_trim <- outlierTrim(tmp7$avg_hr)
    
    tmp7$FEV_trim <- outlierTrim(tmp7$FEV)
    tmp7$FVC_trim <- outlierTrim(tmp7$FVC)
    
    tmp7$DST_trim <- outlierTrim(tmp7$digit_symbol)
    tmp7$VFT_trim <- outlierTrim(tmp7$verbal_total)
    tmp7$Vocab_trim <- outlierTrim(tmp7$vocabulary)
    tmp7$LM_trim <- outlierTrim(tmp7$LM)
    
    tmp7$Alc_trim <- outlierTrim(log(tmp7$units+1))
    tmp7$PckYrs_trim <- outlierTrim(log(tmp7$pack_years+1))
    
    #################################
    ### Create gf and g using PCA ###
    #################################
    
    library(psych)
    tmp7$g = scale(principal(tmp7[,c("DST_trim","VFT_trim","Vocab_trim","LM_trim")], factors=1, rotate="none", na.action="na.exclude")$score)
    
    tmp8 <- tmp7[,c(1,27,28,30:45,50:52)]
    
    
    # par(mfrow=c(4,5))
    # for(i in 4:22){
    # title <- names(tmp8)[i]
    # hist(tmp8[,i], breaks=100, main=title)
    # }
    
    ### create residuals - adjust for age, age^2, sex (and height for FEV and FVC) ###
    tmp9 <- tmp8
    
    for(i in c(4:13, 15:17, 20:22)){
    tmp9[,i] <- scale(resid(lm(tmp9[,i] ~ age + I(age^2) + sex, data=tmp9, na.action="na.exclude")))
    }
    
    for(i in c(18,19)){
    tmp9[,i] <- scale(resid(lm(tmp9[,i] ~ age + I(age^2) + sex + height_trim, data=tmp9, na.action="na.exclude")))
    }
    
    
    ### make sure residuals aren't wildly different from original phenotypes ###
    for(i in c(4:13,15:22)){
    print(names(tmp8)[i])
    print(cor(tmp8[,i], tmp9[,i], use="pairwise.complete.obs"))
    }
    
    ##################################
    ### final dataset for analysis ###
    ##################################
    
    d <- tmp9[,c(1,2,4:13,15:22)]
    names(d) <- c("ID", "age", "glucose", "cholest", "HDL", "sodium", "potassium", "urea", "creatinine", "bmi", "whr", "fat", "sBP", "dBP", "HR", "FEV", "FVC", "Alc", "PckYrs", "g")
    
    original_phenotypes <- tmp8[,c(1,2,4:13,15:22)]
    names(original_phenotypes) <- c("ID", "age", "glucose", "cholest", "HDL", "sodium", "potassium", "urea", "creatinine", "bmi", "whr", "fat", "sBP", "dBP", "HR", "FEV", "FVC", "Alc", "PckYrs", "g")
    
    
    saveRDS(d, file=paste0(config$continuous_traits_data_path, "GS_pheno_resids_20k_", wave, ".rds"))
    list(residuals = d, original = original_phenotypes)
}

results <- lapply(names(waves), function(x) {process_wave(waves[[x]], x)})
names(results) <- names(waves)

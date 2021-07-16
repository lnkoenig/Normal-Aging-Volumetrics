#make R environment for use in SHINY app
#Load the R packages necessary for rest of script
packages <- c("cowplot", "ggseg", "ggplot2", "ggrepel", "tidyverse", "shiny", "lm.beta")
for( i in packages ){
  #  require returns TRUE invisibly if it was able to load package
  if( ! require( i , character.only = TRUE ) ){
    #  If package was not able to be loaded then re-install
    if(i == "ggseg") {
      source("https://neuroconductor.org/neurocLite.R")
      neuro_install('ggseg')
      library('ggseg')
    }else{
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}


#Read in the data

controls <- read.csv("../AIDA/Making_Cohorts/Normals_Cohort/Normals_Cohort_OASIS_DIAN_2019-02-11_with_supra.csv")
pcad <- read.csv("../AIDA/Making_Cohorts/Preclinical_AD/oasis_nonILP_PCAD_list_by_dx_2020-06-02.csv")

#clean controls
controls$AD_status <- 0

controls$Cent_fSUVR_rsf_TOT_CORTMEAN <- NA
controls[controls$source == "OASIS" & controls$tracer == "PIB",]$Cent_fSUVR_rsf_TOT_CORTMEAN <-   45*controls[controls$source == "OASIS" & controls$tracer == "PIB",]$PET_fSUVR_rsf_TOT_CORTMEAN - 47.5 
controls[controls$source == "DIAN" & controls$tracer == "PIB",]$Cent_fSUVR_rsf_TOT_CORTMEAN <- 40.7*controls[controls$source == "DIAN" & controls$tracer == "PIB",]$PET_fSUVR_rsf_TOT_CORTMEAN - 42.9
controls[controls$tracer == "AV45",]$Cent_fSUVR_rsf_TOT_CORTMEAN <-   53.6*controls[controls$tracer == "AV45",]$PET_fSUVR_rsf_TOT_CORTMEAN - 43.2 

controls$cdr_long <- NA


#clean PCAD
pcad <- pcad[pcad$AD_status == 2,]
pcad$AD_status <- 1
pcad$visit <- NA
pcad$DIAN_ID <- NA

pcad[pcad$Sex == 1,]$Sex <- "M"
pcad[pcad$Sex == 0 ,]$Sex <- "F"

pcad <- pcad[pcad$Age>60,]


#convert SUVRs to centiloid
pcad$Cent_fSUVR_rsf_TOT_CORTMEAN <- NA
pcad[pcad$tracer == "PIB",]$Cent_fSUVR_rsf_TOT_CORTMEAN <-   45*pcad[pcad$tracer == "PIB",]$PET_fSUVR_rsf_TOT_CORTMEAN - 47.5 
pcad[pcad$tracer == "AV45",]$Cent_fSUVR_rsf_TOT_CORTMEAN <-   53.6*pcad[pcad$tracer == "AV45",]$PET_fSUVR_rsf_TOT_CORTMEAN - 43.2 



#merge controls and pcad
merged <- rbind(pcad, controls)

#remove columns that are all NA
merged <- merged[colSums(!is.na(merged)) > 0]

#remove columns that are all 0/no variance (but ignores those that aren't numeric)
inx <- sapply(merged, function(x){is.numeric(x)&& var(x,na.rm = T) == 0})
merged <- merged[,which(inx == FALSE|is.na(inx))]



merged$apoe4 <- NA
merged[merged$apoe == "44" | merged$apoe == "34",]$apoe4 <- 1
merged[merged$apoe != "44" & merged$apoe != "34",]$apoe4 <- 0

merged[merged$AD_status == 0,]$AD_status <- "Amyloid-"
merged[merged$AD_status == 1,]$AD_status <- "Amyloid+"

#separate volume and thickness measurements, for 'patient' we're only looking at those within merged age range (would include partially if also
# interested in graphing them)
merged_thickness <- merged[grep("thickness", names(merged), value = FALSE, ignore.case = T)]
merged_volume <- merged[grep("thickness", names(merged), value = FALSE, invert = TRUE, ignore.case = T)]



#remove extra columns from volumes that aren't volumes
non_regions <- c("MR_Session_Label","version","CLINICALDATA.ID","Age","Subject","cdr","mmse","apoe","FS_FSDATA.ID","PUP_PUPTIMECOURSEDATA.ID","tracer","PET_fSUVR_rsf_TOT_CORTMEAN","Scanner","Race","Ethnicity","Education","cdr_long","AD_status","Sex","source","visit","DIAN_ID","Cent_fSUVR_rsf_TOT_CORTMEAN","apoe4")
merged_demo <- merged_volume[, (colnames(merged_volume) %in% non_regions)]
merged_volume <- merged_volume[, !(colnames(merged_volume) %in% non_regions)]



#Combine Hemispheres

#####Sum the Volumes
#combine regions that are left/right in merged
comb_data <- matrix()  
for (i in colnames(merged_volume)) {
  region <- i
  if (grepl("^lh_", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "lh_", fixed = FALSE))))[,2])
    regionb <- paste("rh_", region, sep = "")
  }
  
  if (grepl("^rh_", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "rh_", fixed = FALSE))))[,2])
    regionb <- paste("lh_", region, sep = "")
  }
  if (grepl("^lh", i) == TRUE & grepl("lh_", i) == FALSE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "lh", fixed = FALSE))))[,2])
    regionb <- paste("rh", region, sep = "")
  }
  if (grepl("^rh", i) == TRUE & grepl("rh_", i) == FALSE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "rh", fixed = FALSE))))[,2])
    regionb <- paste("lh", region, sep = "")
  }
  if (grepl("^Left[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Left.", fixed = FALSE))))[,2])
    regionb <- paste("Right.", region, sep = "")
  }
  if (grepl("^Right[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Right.", fixed = FALSE))))[,2])
    regionb <- paste("Left.", region, sep = "")
  }
  if (grepl("^Left[:_:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Left_", fixed = FALSE))))[,2])
    regionb <- paste("Right_", region, sep = "")
  }
  if (grepl("^Right[:_:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Right_", fixed = FALSE))))[,2])
    regionb <- paste("Left_", region, sep = "")
  }
  if (grepl("^L[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "L.", fixed = FALSE))))[,2])
    regionb <- paste("R.", region, sep = "")
  }
  if (grepl("^R[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "R.", fixed = FALSE))))[,2])
    regionb <- paste("L.", region, sep = "")
  }
  if (region %in% colnames(merged_volume)) {
    comb_data <- c(comb_data, merged_volume[region]) #if already have a combined region included, add that
  }else{
    comb_data[region] <- merged_volume[i] + merged_volume[regionb]  #if don't already have a combined region, make one and add
  }
}
comb_data <- data.frame(comb_data, check.names = FALSE) 
comb_data <- comb_data[, !duplicated(colnames(comb_data))] #remove duplicates from comb_data
merged_volume <- comb_data[-1]



##### Average the thicknesses
comb_data <- matrix() 
for (i in colnames(merged_thickness)) {
  region <- i
  if (grepl("^lh_", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "lh_", fixed = FALSE))))[,2])
    regionb <- paste("rh_", region, sep = "")
  }
  
  if (grepl("^rh_", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "rh_", fixed = FALSE))))[,2])
    regionb <- paste("lh_", region, sep = "")
  }
  if (grepl("^lh", i) == TRUE & grepl("lh_", i) == FALSE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "lh", fixed = FALSE))))[,2])
    regionb <- paste("rh", region, sep = "")
  }
  if (grepl("^rh", i) == TRUE & grepl("rh_", i) == FALSE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "rh", fixed = FALSE))))[,2])
    regionb <- paste("lh", region, sep = "")
  }
  if (grepl("^Left[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Left.", fixed = FALSE))))[,2])
    regionb <- paste("Right.", region, sep = "")
  }
  if (grepl("^Right[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Right.", fixed = FALSE))))[,2])
    regionb <- paste("Left.", region, sep = "")
  }
  if (grepl("^Left[:_:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Left_", fixed = FALSE))))[,2])
    regionb <- paste("Right_", region, sep = "")
  }
  if (grepl("^Right[:_:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "Right_", fixed = FALSE))))[,2])
    regionb <- paste("Left_", region, sep = "")
  }
  if (grepl("^L[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "L.", fixed = FALSE))))[,2])
    regionb <- paste("R.", region, sep = "")
  }
  if (grepl("^R[:.:]", i) == TRUE){
    region <- as.character(data.frame(t(data.frame(strsplit(as.character(i), "R.", fixed = FALSE))))[,2])
    regionb <- paste("L.", region, sep = "")
  }
  
  if (region %in% colnames(merged_thickness)) {
    comb_data <- c(comb_data, merged_thickness[region]) #if already have a combined region included, add that
  }else{
    comb_data[region] <- (merged_thickness[i] + merged_thickness[regionb])/2  #if don't already have a combined region, make one and add
  }
}
comb_data <- data.frame(comb_data, check.names = FALSE) 
comb_data <- comb_data[, !duplicated(colnames(comb_data))] #remove duplicates from comb_data - since redoes L and R
merged_thickness <- comb_data[-1]


####ICV correction

#Make subset of the data we're going to correct for head size  
merged_corrected  <- merged_volume[, !grepl(("IntraCranialVol|WM.hypointensities_volume|non.WM.hypointensities_volume"), colnames(merged_volume), ignore.case = T)]

Normals_ICV_mean <- mean(merged_volume$IntraCranialVol, na.rm = TRUE)

ICV_correction_coeff <- data.frame(NA)
for (i in colnames(merged_corrected)) { #goes through each ROI in normals_volume (which contains only volumes)
  Reg <- lm(as.matrix(merged_corrected[i]) ~ as.matrix(merged$IntraCranialVol)) #linear regression determined how the region varies by ICV
  ICV_correction_coeff[i] <- as.numeric(Reg$coefficients[2])
  merged_corrected[i] <- data.frame( as.matrix(merged_corrected[i]) - as.numeric(ICV_correction_coeff[i])*(merged$IntraCranialVol - Normals_ICV_mean )) #normalizes all the normals data based on ICV
}

#add back in the volumes we didn't ICV correct
merged_corrected <- cbind(merged_corrected, merged_volume[,c("IntraCranialVol", "non.WM.hypointensities_volume", "WM.hypointensities_volume")])




#get rid of demographics don't want
merged_demo <- merged_demo[,c("Age", "Sex", 'source', 'apoe4', "AD_status")]

#### combine thicknesses and volumes
normals  <- data.frame(merged_demo, merged_corrected, merged_thickness)

#clear out environment
rm(list = ls()[!(ls() %in% c('normals','ICV_correction_coeff'))])


normals <- relocate(normals, bankssts_thickness, caudalanteriorcingulate_thickness, caudalmiddlefrontal_thickness, cuneus_thickness, entorhinal_thickness, frontalpole_thickness, fusiform_thickness, inferiorparietal_thickness, inferiortemporal_thickness, insula_thickness, isthmuscingulate_thickness, lateraloccipital_thickness, lateralorbitofrontal_thickness, lingual_thickness, medialorbitofrontal_thickness, middletemporal_thickness, paracentral_thickness, parahippocampal_thickness, parsopercularis_thickness, parsorbitalis_thickness, parstriangularis_thickness, pericalcarine_thickness, postcentral_thickness, posteriorcingulate_thickness, precentral_thickness, precuneus_thickness, rostralanteriorcingulate_thickness, rostralmiddlefrontal_thickness, superiorfrontal_thickness, superiorparietal_thickness, superiortemporal_thickness, supramarginal_thickness, temporalpole_thickness, transversetemporal_thickness, bankssts_volume, caudalanteriorcingulate_volume, caudalmiddlefrontal_volume, cuneus_volume, entorhinal_volume, frontalpole_volume, fusiform_volume, inferiorparietal_volume, inferiortemporal_volume, insula_volume, isthmuscingulate_volume, lateraloccipital_volume, lateralorbitofrontal_volume, lingual_volume, medialorbitofrontal_volume, middletemporal_volume, paracentral_volume, parahippocampal_volume, parsopercularis_volume, parsorbitalis_volume, parstriangularis_volume, pericalcarine_volume, postcentral_volume, posteriorcingulate_volume, precentral_volume, precuneus_volume, rostralanteriorcingulate_volume, rostralmiddlefrontal_volume, superiorfrontal_volume, superiorparietal_volume, superiortemporal_volume, supramarginal_volume, temporalpole_volume, transversetemporal_volume, Amygdala_volume, Caudate_volume, Hippocampus_volume, Lateral.Ventricle_volume, Pallidum_volume, Putamen_volume, Thalamus.Proper_volume, VentralDC_volume, Accumbens.area_volume, Brain.Stem_volume, CC_Anterior_volume, CC_Central_volume, CC_Mid_Anterior_volume, CC_Mid_Posterior_volume, CC_Posterior_volume, Cerebellum.Cortex_volume, Cerebellum.White.Matter_volume, choroid.plexus_volume, CortexVol, CorticalWhiteMatterVol, CSF_volume, Inf.Lat.Vent_volume, IntraCranialVol, Optic.Chiasm_volume, SubCortGrayVol, SupraTentorialVol, TotalGrayVol, vessel_volume, WM.hypointensities_volume, non.WM.hypointensities_volume, X3rd.Ventricle_volume, X4th.Ventricle_volume, X5th.Ventricle_volume)

normals <- dplyr::rename(normals, "Banks STS Thickness" = bankssts_thickness, "Caudal Anterior Cingulate Thickness" = caudalanteriorcingulate_thickness, "Caudal Middle Frontal Thickness" = caudalmiddlefrontal_thickness, "Cuneus Thickness" = cuneus_thickness, "Entorhinal Thickness" = entorhinal_thickness, "Frontal Pole Thickness" = frontalpole_thickness, "Fusiform Thickness" = fusiform_thickness, "Inferior Parietal Thickness" = inferiorparietal_thickness, "Inferior Temporal Thickness" = inferiortemporal_thickness,  "Insula Thickness" = insula_thickness, "Isthmus Cingulate Thickness" = isthmuscingulate_thickness, "Lateral Occipital Thickness" = lateraloccipital_thickness,  "Lateral Orbitofrontal Thickness" = lateralorbitofrontal_thickness, "Lingual Thickness" = lingual_thickness, "Medial Orbitofrontal Thickness" = medialorbitofrontal_thickness, "Middle Temporal Thickness" = middletemporal_thickness, "Paracentral Thickness" = paracentral_thickness, "Parahippocampal Thickness" = parahippocampal_thickness, "Pars Opercularis Thickness" = parsopercularis_thickness, "Pars Orbitalis Thickness" = parsorbitalis_thickness, "Pars Triangularis Thickness" = parstriangularis_thickness, "Pericalcarine Thickness" = pericalcarine_thickness, "Postcentral Thickness" = postcentral_thickness, "Posterior Cingulate Thickness" = posteriorcingulate_thickness, "Precentral Thickness" = precentral_thickness, "Precuneus Thickness" = precuneus_thickness, "Rostral Anterior Cingulate Thickness" = rostralanteriorcingulate_thickness, "Rostral Middle Frontal Thickness" = rostralmiddlefrontal_thickness, "Superior Frontal Thickness" = superiorfrontal_thickness, "Superior Parietal Thickness" = superiorparietal_thickness, "Superior Temporal Thickness" = superiortemporal_thickness, "Supramarginal Thickness" = supramarginal_thickness, "Temporal Pole Thickness" = temporalpole_thickness, "Transverse Temporal Thickness" = transversetemporal_thickness)


normals <- dplyr::rename(normals, "Banks STS Volume" = bankssts_volume, "Caudal Anterior Cingulate Volume" = caudalanteriorcingulate_volume, "Caudal Middle Frontal Volume" = caudalmiddlefrontal_volume, "Cuneus Volume" = cuneus_volume, "Entorhinal Volume" = entorhinal_volume, "Frontal Pole Volume" = frontalpole_volume, "Fusiform Volume" = fusiform_volume, "Inferior Parietal Volume" = inferiorparietal_volume, "Inferior Temporal Volume" = inferiortemporal_volume, "Insula Volume" = insula_volume, "Isthmus Cingulate Volume" = isthmuscingulate_volume, "Lateral Occipital Volume" = lateraloccipital_volume, "Lateral Orbitofrontal Volume" = lateralorbitofrontal_volume, "Lingual Volume" = lingual_volume, "Medial Orbitofrontal Volume" = medialorbitofrontal_volume, "Middle Temporal Volume" = middletemporal_volume, "Paracentral Volume" = paracentral_volume, "Parahippocampal Volume" = parahippocampal_volume, "Pars Opercularis Volume" = parsopercularis_volume, "Pars Orbitalis Volume" = parsorbitalis_volume, "Pars Triangularis Volume" = parstriangularis_volume, "Pericalcarine Volume" = pericalcarine_volume, "Postcentral Volume" = postcentral_volume, "Posterior Cingulate Volume" = posteriorcingulate_volume, "Precentral Volume" = precentral_volume, "Precuneus Volume" = precuneus_volume, "Rostral Anterior Cingulate Volume" = rostralanteriorcingulate_volume, "Rostral Middle Frontal Volume" = rostralmiddlefrontal_volume, "Superior Frontal Volume" = superiorfrontal_volume, "Superior Parietal Volume" = superiorparietal_volume, "Superior Temporal Volume" = superiortemporal_volume, "Supramarginal Volume" = supramarginal_volume, "Temporal Pole Volume" = temporalpole_volume, "Transverse Temporal Volume" = transversetemporal_volume, "Amygdala Volume" = Amygdala_volume, "Caudate Volume" = Caudate_volume, "Hippocampus Volume" = Hippocampus_volume, "Lateral Ventricles Volume" = Lateral.Ventricle_volume, "Pallidum Volume" = Pallidum_volume, "Putamen Volume" = Putamen_volume, "Thalamus Volume" = Thalamus.Proper_volume, "Ventral DC Volume" = VentralDC_volume, "Nucleus Accumbens Volume" = Accumbens.area_volume, "Brain Stem Volume" = Brain.Stem_volume, "Anterior Cingulate Cortex Volume" = CC_Anterior_volume, "Central Cingulate Cortex Volume" = CC_Central_volume, "Mid Anterior Cingulate Cortex Volume" = CC_Mid_Anterior_volume, "Mid-Posterior Cingulate Cortex Volume" = CC_Mid_Posterior_volume, "Posterior Cingulate Cortex Volume" = CC_Posterior_volume,  "Cerebellum Cortex Volume" = Cerebellum.Cortex_volume, "Cerebellum White Matter Volume" = Cerebellum.White.Matter_volume, "Choroid Plexus Volume" = choroid.plexus_volume, "Cortex Volume" = CortexVol, "Cortical White Matter Volume" = CorticalWhiteMatterVol, "CSF Volume" = CSF_volume, "Inferior Lateral Ventricles Volume" = Inf.Lat.Vent_volume, "Intracranial Volume" = IntraCranialVol, "Optic Chiasm Volume" = Optic.Chiasm_volume, "Subcortical Grey Matter Volume" = SubCortGrayVol, "Supratentorial Volume" = SupraTentorialVol, "Total Gray Matter Volume" = TotalGrayVol, "Vessel Volume" = vessel_volume, "White Matter Hyperintensities Volume" = WM.hypointensities_volume, "Non-WM Hyperintensities Volume" = non.WM.hypointensities_volume, "3rd Ventricle Volume" = X3rd.Ventricle_volume, "4th Ventricle Volume" = X4th.Ventricle_volume, "5th Ventricle Volume" = X5th.Ventricle_volume)


# for(col in colnames(select(normals, "Banks STS Thickness":"5th Ventricle Volume"))){
#   c <- strsplit(gsub("\\.|_", " ", col, ignore.case = FALSE, perl = FALSE,
#                      fixed = FALSE, useBytes = FALSE), " ")[[1]] 
#   region_title <- paste(toupper(substring(c, 1,1)), substring(c, 2),
#                         sep = "", collapse = " ")
#   colnames(normals)[colnames(normals) == col] <- region_title
#   colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff) == col] <- region_title
# }


ICV_correction_coeff <- dplyr::rename(ICV_correction_coeff, "Banks STS Volume" = bankssts_volume, "Caudal Anterior Cingulate Volume" = caudalanteriorcingulate_volume, "Caudal Middle Frontal Volume" = caudalmiddlefrontal_volume, "Cuneus Volume" = cuneus_volume, "Entorhinal Volume" = entorhinal_volume, "Frontal Pole Volume" = frontalpole_volume, "Fusiform Volume" = fusiform_volume, "Inferior Parietal Volume" = inferiorparietal_volume, "Inferior Temporal Volume" = inferiortemporal_volume, "Insula Volume" = insula_volume, "Isthmus Cingulate Volume" = isthmuscingulate_volume, "Lateral Occipital Volume" = lateraloccipital_volume, "Lateral Orbitofrontal Volume" = lateralorbitofrontal_volume, "Lingual Volume" = lingual_volume, "Medial Orbitofrontal Volume" = medialorbitofrontal_volume, "Middle Temporal Volume" = middletemporal_volume, "Paracentral Volume" = paracentral_volume, "Parahippocampal Volume" = parahippocampal_volume, "Pars Opercularis Volume" = parsopercularis_volume, "Pars Orbitalis Volume" = parsorbitalis_volume, "Pars Triangularis Volume" = parstriangularis_volume, "Pericalcarine Volume" = pericalcarine_volume, "Postcentral Volume" = postcentral_volume, "Posterior Cingulate Volume" = posteriorcingulate_volume, "Precentral Volume" = precentral_volume, "Precuneus Volume" = precuneus_volume, "Rostral Anterior Cingulate Volume" = rostralanteriorcingulate_volume, "Rostral Middle Frontal Volume" = rostralmiddlefrontal_volume, "Superior Frontal Volume" = superiorfrontal_volume, "Superior Parietal Volume" = superiorparietal_volume, "Superior Temporal Volume" = superiortemporal_volume, "Supramarginal Volume" = supramarginal_volume, "Temporal Pole Volume" = temporalpole_volume, "Transverse Temporal Volume" = transversetemporal_volume, "Amygdala Volume" = Amygdala_volume, "Caudate Volume" = Caudate_volume, "Hippocampus Volume" = Hippocampus_volume, "Lateral Ventricles Volume" = Lateral.Ventricle_volume, "Pallidum Volume" = Pallidum_volume, "Putamen Volume" = Putamen_volume, "Thalamus Volume" = Thalamus.Proper_volume, "Ventral DC Volume" = VentralDC_volume, "Nucleus Accumbens Volume" = Accumbens.area_volume, "Brain Stem Volume" = Brain.Stem_volume, "Anterior Cingulate Cortex Volume" = CC_Anterior_volume, "Central Cingulate Cortex Volume" = CC_Central_volume, "Mid Anterior Cingulate Cortex Volume" = CC_Mid_Anterior_volume, "Mid-Posterior Cingulate Cortex Volume" = CC_Mid_Posterior_volume, "Posterior Cingulate Cortex Volume" = CC_Posterior_volume,  "Cerebellum Cortex Volume" = Cerebellum.Cortex_volume, "Cerebellum White Matter Volume" = Cerebellum.White.Matter_volume, "Choroid Plexus Volume" = choroid.plexus_volume, "Cortex Volume" = CortexVol, "Cortical White Matter Volume" = CorticalWhiteMatterVol, "CSF Volume" = CSF_volume, "Inferior Lateral Ventricles Volume" = Inf.Lat.Vent_volume, "Optic Chiasm Volume" = Optic.Chiasm_volume, "Subcortical Grey Matter Volume" = SubCortGrayVol, "Supratentorial Volume" = SupraTentorialVol, "Total Gray Matter Volume" = TotalGrayVol, "Vessel Volume" = vessel_volume, "3rd Ventricle Volume" = X3rd.Ventricle_volume, "4th Ventricle Volume" = X4th.Ventricle_volume, "5th Ventricle Volume" = X5th.Ventricle_volume)


#save data for shiny app
filename = paste("Normal_Aging_Input_Shiny_", Sys.Date(), ".RData", sep = "") #making csv to write data into
save("normals", "ICV_correction_coeff",file = filename) # save environment





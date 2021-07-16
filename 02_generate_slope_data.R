#Estimates Normal Aging Slopes and correlates slopes with age to find acceleration

#Read in data
merged_new <- read.csv("data-clean/clean_Normal_Aging_and_PCAD_cohorts_2020-06-02.csv")

#separate out normal aging cohort
normals <- merged_new[merged_new$AD_status==0,]


region_list <- colnames(subset(merged_new, select = -c(AD_status, CLINICALDATA.ID, Subject, cdr, mmse, apoe, apoe4, FS_FSDATA.ID, PUP_PUPTIMECOURSEDATA.ID, tracer, Cent_fSUVR_rsf_TOT_CORTMEAN, PET_fSUVR_rsf_TOT_CORTMEAN, MR_Session_Label, Scanner, Age, version, Race, Ethnicity, Education, Sex, source, visit, DIAN_ID, PET_fSUVR_rsf_TOT_CORTMEAN, cdr_long) ) )


#Set normals data's min and max ages
n_min <- ceiling(min(normals$Age))
n_max <- floor(max(normals$Age))

age_range <- 1 #range at which slopes are calculated (1 = 1 year bin)


#make dataframe to hold results
slopes <- data.frame(Age=(n_min+1):n_max)
slopes_m <- data.frame(Age=(n_min+1):n_max)

for (ROI_name in region_list) {      #Loop that goes through each region in subset matrix defined above
  
  #Smooth population data with loess
  smoothFit <- loess( as.matrix(normals[ROI_name]) ~ normals$Age, degree=1, span=0.7 ) #locally fits a polynomial surface
  smoothPred <- predict( smoothFit, newdata=seq(n_min, n_max, by=age_range), se=TRUE ) #uses loess fit to predict data for each age
  
  stand_dev2 <- matrix()
  
  
  for (i in 1:length(n_min:n_max)) {
    stand_dev2[i] <- sd(as.matrix(normals[normals$Age > seq(n_min,n_max, by=1)[i]-1 & normals$Age < seq(n_min,n_max, by=1)[i]+1,][ROI_name]))
  }
  sd_fit <- loess(stand_dev2 ~ c(n_min:n_max), degree=1, span=0.7)
  sd_pred <- predict(sd_fit, newdata=seq(n_min, n_max, by=age_range))
  
  
  plotFrame <- data.frame( Age=c(seq(n_min, n_max, by=age_range)), fit=smoothPred$fit, sd=sd_pred ) #combines all the model's data into single dataframe

  slope <- diff(plotFrame$fit, differences=1)/age_range #normalize to range fit the data (if using an age increase other than 1 year it switches it to per year)
  slope <- (slope) / plotFrame$fit[-length(plotFrame$fit)] #make it percent change

  #save data to dataframe
  slopes[,ROI_name] <- slope
  
  
  #REPEAT WITH MERGED DATA
  
  #Smooth population data with loess
  smoothFit <- loess( as.matrix(merged_new[ROI_name]) ~ merged_new$Age, degree=1, span=0.7 ) #locally fits a polynomial surface
  smoothPred <- predict( smoothFit, newdata=seq(n_min, n_max, by=age_range), se=TRUE ) #uses loess fit to predict data for each age
  
  stand_dev2 <- matrix()
  
  
  for (i in 1:length(n_min:n_max)) {
    stand_dev2[i] <- sd(as.matrix(merged_new[merged_new$Age > seq(n_min,n_max, by=1)[i]-1 & merged_new$Age < seq(n_min,n_max, by=1)[i]+1,][ROI_name]))
  }
  sd_fit <- loess(stand_dev2 ~ c(n_min:n_max), degree=1, span=0.7)
  sd_pred <- predict(sd_fit, newdata=seq(n_min, n_max, by=age_range))
  
  
  plotFrame <- data.frame( Age=c(seq(n_min, n_max, by=age_range)), fit=smoothPred$fit, sd=sd_pred ) #combines all the model's data into single dataframe
  
  slope <- diff(plotFrame$fit, differences=1)/age_range #normalize to range fit the data (age increased by 0.1, switching to per year)
  slope <- (slope - mean(slope)) / sd(slope) #normalized to the size of region at that age

  #save data to dataframe
  slopes_m[,ROI_name] <- slope
  
  
}


write.csv(slopes, paste("data-clean/NormalAgingCohort_Slopes_", Sys.Date(), ".csv", sep=""), row.names = F)
write.csv(slopes_m, paste("data-clean/MergedCohort_Slopes_", Sys.Date(), ".csv", sep=""), row.names = F)

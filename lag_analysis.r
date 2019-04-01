library(astsa)
library(neurobase)
library(parallel)
library(zoo)
library(RNifti)
library(oro.nifti)

#generic setup
BIDS_proc_folder="/Volumes/RAID_NP/Nikos/LOCUS1/Locus1_proc/BIDS_nii_structure/BIDS_proc"
subj_list<-as.matrix(read.table("/Users/nikos/Locus1/BIDS_nii_structure/subject_list_complete.txt",sep='\t', header=F, na.strings=c("NA","NaN", " ")))
task="rs1"
fmri_proce="fixed_slicecorr_cropped_smoothed"

#assign max lag
max_lag=5

for (subject in 19:length(subj_list)){
  
  fmri_name<-paste(BIDS_proc_folder,"/",subj_list[subject],"/func/",subj_list[subject],"_",task,"_",fmri_proce,".nii.gz", collapse = NULL, sep="");
  output<-paste(BIDS_proc_folder,"/",subj_list[subject],"/func/", collapse = NULL, sep="");
  output_lags<-paste(output,subj_list[subject],"_",task,"_",fmri_proce,"_lags", collapse = NULL, sep="");
  output_lagged_fmri<-paste(output,subj_list[subject],"_",task,"_",fmri_proce,"_lagged", collapse = NULL, sep="");
  reference_ts_file<-paste(BIDS_proc_folder,"/",subj_list[subject],"/func/",subj_list[subject],"_",task,"_roi/epi_vol_pve_2_",fmri_proce,".txt", collapse = NULL, sep="")
  mask_for_nii_path<-paste("/Users/nikos/Locus1/BIDS_nii_structure/BIDS_proc/",subj_list[subject],"/func/",subj_list[subject],"_",task,"_roi","/hip_amy_ec_pons_whole.nii.gz" , collapse = NULL, sep="")
  
  #assign reference to average GM timeseries
  reference_ts<-as.matrix(read.table(reference_ts_file,sep='\t', header=F, na.strings=c("NA","NaN", " ")))
  #turn it into a timeseries
  reference_ts<-ts(as.vector(reference_ts))
  
  #read fmri and get dimension
  fmri<-readNifti(fmri_name,internal=FALSE)
  #fmri<-readnii(fmri_name)
  
  dim_fmri<-dim(fmri)[1:3]

  #load mask 
  mask = readnii(mask_for_nii_path)
  mask_2d<-t(c(mask))
  mask_2d[mask_2d==0]<-NA
  #fmri<-fmri[1:5,1:5,1:5,1:5]
  
  #turn 4d timeseries to a 2D object (turn to voxelxtimepoints)
  vox_ts <- t(apply(fmri, 4, c))
  #make NA the ones that are of no interest
  vox_ts[,is.na(mask_2d)]<-NA
  
  #define function that finds optimal lag
  find_max_lag <- function(timeseries_obj){
    #for one timeseries
    if (!is.na(timeseries_obj[1])){
      ccf_obj=ccf(ts(as.vector(timeseries_obj)), reference_ts, max_lag,ylab="CCF",plot=FALSE)
    optimal_lag_for_voxel<-ccf_obj$lag[which.max(abs(ccf_obj$acf))]
    return(optimal_lag_for_voxel)
    } else {return(NA)}
  }
  
  no_cores <- detectCores() - 1
  #correct lag function
  cl <- makeCluster(no_cores, type="FORK")
  optimal_lags<-parApply(cl,vox_ts, 2,find_max_lag)
  stopCluster(cl)
  
  #get rid of NAs if NAs
  optimal_lags<-lapply(optimal_lags, function(x) x[1])
  optimal_lags<-unlist(optimal_lags)
  optimal_lags[is.na(optimal_lags)]<-0
  
  #construct lags to 3D array
  optimal_lags_3D<-array(optimal_lags, dim_fmri)
  
  #roll fMRI timeseries according to lag
  vox_ts<-rbind(vox_ts,optimal_lags)
  
  lag_timeseries <- function(timeseries_obj){
    tmp1=length(timeseries_obj)
    if (!is.na(timeseries_obj[1])){
    return(matrix(lag(zoo(timeseries_obj[1:(tmp1-1)]),timeseries_obj[tmp1],na.pad = TRUE)))
    } else { return(matrix(0,tmp1-1))}
  }
  

  cl <- makeCluster(no_cores, type="FORK")
  lagged_fmri_ts_2D<-parApply(cl,vox_ts,2,lag_timeseries)
  stopCluster(cl)
  
  # crop lagged fmri ts
  ts_length<-dim(lagged_fmri_ts_2D)[1]-max_lag
  lagged_fmri_ts_2D<-lagged_fmri_ts_2D[max_lag:(ts_length-max_lag),]
  # and put them in 4D
  lagged_fmri_ts_4D<-array(t(lagged_fmri_ts_2D), c(dim_fmri,(ts_length-max_lag)))
  
  ## SAVE into nii and you are done
    #here for the lags, for some weird reason it tends to output high intensity values 
    #i.e. 255 in the nifti format when using other writing functions than writenii. 
    #bit disturbing
  lags.nifti <- nifti(optimal_lags_3D)
  nim = copyNIfTIHeader(img = mask, arr = optimal_lags_3D)
  writenii(nim,output_lags)

  fmri_lags.nifti <- nifti(lagged_fmri_ts_4D,datatype=8)
  #nimfmri = copyNIfTIHeader(img = fmri, arr = lagged_fmri_ts_4D)
  #writenii(fmri_lags.nifti, output_lagged_fmri)  
  writeNIfTI(fmri_lags.nifti, output_lagged_fmri, verbose=TRUE)
  rm(fmri_lags.nifti,lagged_fmri_ts_4D,fmri,vox_ts, lagged_fmri_ts_2D)
  rm(optimal_lags, optimal_lags_3D,mask,nim)
  
}
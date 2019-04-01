#!/bin/bash

#thsi script offers a basic fmri preprocessing
#we do
	#motion correction with mc flirt towards the last volume (closest to PA scan for subsequent optimal topup correction)
	#motion outliers with fsl_motion_outliers. 
		#current measure: FD , task threshold: 0.9, resting state threshold: 0.4
		#motion confound mstrix NOT created if no volume exceed threshold!
	#topup: is done assuming 5 volumes for PA and that the PA scan is done after the encoding scan 
		#(you can change the order from the function and the number of volumes from the acq file)
	#melodic: is done by loading up the gui (you need to press "Go")
		#includes BBR, FNIRT to 1mm MNI, 100Hz highpass filtering, no smoothing and melodic
	#fix: is done with 30 (found optimal in our testing)
	#pnm:

PROGNAME=$(basename $0)

function error_exit
{
#	----------------------------------------------------------------
#	Function for error control
#		Accepts 1 argument:
#			string containing descriptive error message
#	----------------------------------------------------------------
	echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}


motion_correction () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4

	mkdir -p ${output}/${subject}/func/tmp
	mkdir -p ${output}/${subject}/func/tmp/mc


	num_vol=`fslval ${input}/${subject}/func/${subject}_${task}_bold.nii.gz dim4`
	num_vol=`expr $num_vol - 1`


	in_nii=${input}/${subject}/func/${subject}_${task}_bold.nii.gz
	out_nii="${output}/${subject}/func/tmp/mc/${subject}_${task}_bold_moco.nii.gz"

	(mcflirt -in $in_nii -out $out_nii -refvol $num_vol -mats -plots) || error_exit "cannot motion correct $in_nii"
	mv $out_nii ${output}/${subject}/func/${subject}_${task}_bold_corrected.nii.gz
}

motion_outliers () {

	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local threshold=$5

	mkdir -p ${output}/${subject}/func/tmp
	mkdir -p ${output}/${subject}/func/tmp/motion_outliers

	in_nii=${input}/${subject}/func/${subject}_${task}_bold.nii.gz
	out_nii="${output}/${subject}/func/tmp/motion_outliers/${subject}_${task}_bold"

	(fsl_motion_outliers -i $in_nii -o ${out_nii}_motion_confound.txt \
	-p ${out_nii}_FDplot.png -s ${out_nii}_FDvalues.txt --fd --thresh=$threshold) || error_exit "cannot do fsl_motion_outliers for $in_nii"

}

topup_apply () {

	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local acq_param=$5

	mkdir -p ${output}/${subject}/func/tmp
	mkdir -p ${output}/${subject}/func/tmp/topup

	cp ${output}/${subject}/func/${subject}_${task}_bold_corrected.nii.gz ${output}/${subject}/func/tmp/topup/${subject}_${task}_bold_corrected.nii.gz || error_exit "topup input fail 1"
	cp ${input}/${subject}/fmap/${subject}_${task}_PA.nii.gz ${output}/${subject}/func/tmp/topup/${subject}_${task}_PA.nii.gz || error_exit "topup input fail 2"

	dirpath=${output}/${subject}/func/tmp/topup
	echo $dirpath
	cd $dirpath

	file_to_slice=${subject}_${task}_bold_corrected.nii.gz;
	echo $file_to_slice
	#file to correct with (blip down)
	file_to_correct_with=${subject}_${task}_PA.nii.gz;
	echo $file_to_correct_with
	#order of scans (PA then scan or opposite)
	order_variable=2;

	#get dim4 dimensions of image
	a=`fslinfo ${file_to_slice} | grep dim4 | awk '{print $2;}'`;
	a=`echo $a|awk '{print $1;}'`;
	#if you have PA volume number=10, then you need 10 volumes
	#modify otherwise
	#and modify acq_param too
	let a=a-5;


	if   [ "${order_variable}" -eq "1" ]; then
		pos_slice=0;
	elif [ "${order_variable}" -eq "2" ]; then
		pos_slice=${a};
	else
		echo "wrong input";
	fi

	#cut away relevant part of the file for correction
	fslroi ${file_to_slice} AP_${file_to_slice} 0 -1 0 -1 0 -1 ${pos_slice} 5 || error_exit "topup fslroi fail";

	fname=`$FSLDIR/bin/remove_ext ${file_to_slice}`;

	#merge blip up and blip down to one
	fslmerge -t input_topup_${file_to_slice} AP_${file_to_slice} ${file_to_correct_with} || error_exit "topup merge fail";

	#see if we need to cut away top slice (topup subsampling prerequisite)
	z_slice_num=`fslinfo input_topup_${file_to_slice} | grep dim3 | awk '{print $2;}'`;
	z_slice_num=`echo $z_slice_num|awk '{print $1;}'`;

	if   [ $((z_slice_num%2)) -ne 0 ]; then
		let z_slice_num=z_slice_num-1;
		fslroi input_topup_${file_to_slice} input_topup_${file_to_slice} 0 -1 0 -1 0 ${z_slice_num} 0 -1
		fslroi ${file_to_slice} ${file_to_slice} 0 -1 0 -1 0 ${z_slice_num} 0 -1
	fi

	#calculate epi distortion with topup
	topup --imain=input_topup_${file_to_slice} --datain=$acq_param \
	--config=$FSLDIR/etc/flirtsch/b02b0.cnf --out=${fname}_topup  --iout=${fname}_hifi_topup || error_exit "topup distortion estimation fail";

	#apply distortion correction to file
	applytopup --imain=${file_to_slice} --inindex=1 --method=jac --datain=$acq_param \
	--topup=${fname}_topup --out=${fname}_corrected || error_exit "topup applytopup fail";

	echo "done with" 
	echo ${subject}_${task};

	mv ${fname}_corrected.nii.gz ${output}/${subject}/func/${subject}_${task}_bold_topupcorrected.nii.gz

}


run_melodic () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local melodic_template=$5

	mkdir -p ${output}/${subject}/func/tmp
	cp -rf $melodic_template ${output}/${subject}/func/tmp/${task}design_melodic.fsf

	echo $subject $task "melodic"

	input_file=${output}/${subject}/func/${subject}_${task}_bold_topupcorrected.nii.gz
	nvol_fmri=`fslval $input_file dim4`
	output_melodic=${output}/${subject}/func/${subject}_${task}_melodic
	input_brain=${input}/${subject}/anat/${subject}_T1w_brain.nii.gz

	echo $input_file
	echo $nvol_fmri
	echo $output_melodic
	echo $input_brain
	sed -e "s@OUTPUT_FILE@$output_melodic@g; s@NVOL@$nvol_fmri@g; s@INPUT_FMRI_FILE@$input_file@g; s@INPUT_T1_BRAIN@$input_brain@g" ${output}/${subject}/func/tmp/${task}design_melodic.fsf\
	 > ${output}/${subject}/func/tmp/${subject}_${task}_design_melodic.fsf \
	|| error_exit "cannot do melodic template for $subject $task"
	wait ${!}
	feat ${output}/${subject}/func/tmp/${subject}_${task}_design_melodic.fsf
}

run_fix () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4

	echo $subject $task "fix"
	#get motion parameters to melodic file
	#then do fix and cleanup
	mel_file=${output}/${subject}/func/${subject}_${task}_melodic.ica
	mkdir -p ${mel_file}/mc
	echo $mel_file
	#echo ${mel_file}/mc/prefiltered_func_data_mcf.par
	cp ${output}/${subject}/func/tmp/mc/${subject}_${task}_bold_moco.nii.gz.par ${mel_file}/mc/prefiltered_func_data_mcf.par \
	|| error_exit "cannot copy motion file for fix $subject $task"
	fix $mel_file /usr/local/fix/training_files/HCP_hp2000.RData 30 || error_exit "cannot do fix $subject $task"

	cp ${output}/${subject}/func/${subject}_${task}_melodic.ica/filtered_func_data_clean.nii.gz ${output}/${subject}/func/${subject}_${task}_fixed.nii.gz

}


pnm_regressors () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4

	echo $subject $task "make phys regressors"
	phys_file=${output}/${subject}/func/physio
	input_pnm_file=${phys_file}/${subject}_${task}_physio.txt

	pnm_file=${output}/${subject}/func/physio/${subject}_${task}_pnm/pnm
	mkdir -p ${output}/${subject}/func/physio/${subject}_${task}_pnm

	/usr/local/fsl/bin/fslFixText $input_pnm_file ${pnm_file}_input.txt
	/usr/local/fsl/bin/pnm_stage1 -i ${pnm_file}_input.txt -o ${pnm_file} -s 50 --tr=2.0 --smoothcard=0.1 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3

	/usr/local/fsl/bin/fslFixText $input_pnm_file ${pnm_file}_input.txt


	/usr/local/fsl/bin/pnm_stage1 -i ${pnm_file}_input.txt -o ${pnm_file} -s 50 --tr=2.0 --smoothcard=0.1 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3
	/usr/local/fsl/bin/popp -i ${pnm_file}_input.txt -o ${pnm_file} -s 50 --tr=2.0 --smoothcard=0.1 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3
	
	open ${pnm_file}_pnm1.html 
	wait ${!}

	cd $pwd
	obase=$pnm_file
	# if [ $# -gt 0 ] ; then 
	# 	/usr/local/fsl/bin/popp -i ${pnm_file}_input.txt \
	# 	-o ${pnm_file} -s 50 --tr=2.0 --smoothcard=0.1 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3 $@ ; 
	# fi
		/usr/local/fsl/bin/pnm_evs -i ${output}/${subject}/func/${subject}_${task}_bold_topupcorrected.nii.gz \
		-c ${pnm_file}_card.txt -r ${pnm_file}_resp.txt \
		-o $pnm_file --tr=2.0 --oc=3 --or=4 --multc=1 --multr=1 \
		--slicetiming=${input}/slicetimes_unix3.txt
	ls -1 `/usr/local/fsl/bin/imglob -extensions ${obase}ev0*` > ${pnm_file}_evlist.txt

	#get respiration phase regressor (1 and -1 depending on respiration/exhalation; don't remember which one is which)
	phase=${obase}ev001.nii.gz
	fslmaths $phase -thr 0 -mul 500000 -bin ${obase}phase.nii.gz
	fslmaths $phase -uthr 0 -mul -500000 -bin -mul -1 -add ${obase}phase.nii.gz ${obase}phase.nii.gz
	#and add it to the ev list
	ls ${obase}phase.nii.gz >> ${pnm_file}_evlist.txt

}

pnm_regressors_respiration_only () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4

	echo $subject $task "make phys regressors"
	phys_file=${output}/${subject}/func/physio
	input_pnm_file=${phys_file}/${subject}_${task}_physio.txt

	pnm_file=${output}/${subject}/func/physio/${subject}_${task}_pnm/pnm_respiration
	mkdir -p ${output}/${subject}/func/physio/${subject}_${task}_pnm/pnm_respiration

	/usr/local/fsl/bin/fslFixText $input_pnm_file ${pnm_file}_input.txt
	/usr/local/fsl/bin/pnm_stage1 -i ${pnm_file}_input.txt -o ${pnm_file} -s 50 --tr=2.0 --smoothcard=0.1 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3

	/usr/local/fsl/bin/popp -i ${pnm_file}_input.txt -o ${pnm_file} -s 50 --tr=2.0 --smoothcard=0.1 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3
	
	open ${pnm_file}_pnm1.html 
	wait ${!}

	cd $pwd
	obase=$pnm_file
	# if [ $# -gt 0 ] ; then 
	# 	/usr/local/fsl/bin/popp -i ${pnm_file}_input.txt \
	# 	-o ${pnm_file} -s 50 --tr=2.0 --smoothcard=0.1 --smoothresp=0.1 --resp=1 --cardiac=2 --trigger=3 $@ ; 
	# fi
		/usr/local/fsl/bin/pnm_evs -i ${output}/${subject}/func/${subject}_${task}_bold_topupcorrected.nii.gz \
		-c ${pnm_file}_card.txt -r ${pnm_file}_resp.txt \
		-o $pnm_file --tr=2.0 --oc=0 --or=4 --multc=0 --multr=0 \
		--slicetiming=${input}/slicetimes_unix3.txt

	
	ls -1 `/usr/local/fsl/bin/imglob -extensions ${obase}ev0*` > ${pnm_file}_evlist.txt

	#get respiration phase regressor (1 and -1 depending on respiration/exhalation; don't remember which one is which)
	phase=${obase}ev001.nii.gz
	fslmaths $phase -thr 0 -mul 500000 -bin ${obase}phase.nii.gz
	fslmaths $phase -uthr 0 -mul -500000 -bin -mul -1 -add ${obase}phase.nii.gz ${obase}phase.nii.gz
	#and add it to the ev list
	ls ${obase}phase.nii.gz >> ${pnm_file}_evlist.txt

}


run_pnm_feat () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local melodic_template=$5

	mkdir -p ${output}/${subject}/func/tmp
	cp -rf $melodic_template ${output}/${subject}/func/tmp/${task}design_pnmfeat.fsf

	echo $subject $task "feat PNM"
	#maybe add a check that FIX has run already correctly!
	input_file=${output}/${subject}/func/${subject}_${task}_fixed.nii.gz
	nvol_fmri=`fslval $input_file dim4`
	output_melodic=${output}/${subject}/func/${subject}_${task}_feat_pnm_respphase
	input_brain=${input}/${subject}/anat/${subject}_T1w_brain.nii.gz

	pnm_evlist=${output}/${subject}/func/physio/${subject}_${task}_pnm/pnm_evlist.txt

	conf_list=${output}/${subject}/func/${subject}_${task}_confound_M_PR_WM_CSF.txt

	sed -e "s@OUTPUT_DIR@$output_melodic@g; s@NVAL@$nvol_fmri@g; s@INPUT_BOLD@$input_file@g; s@VOXELWISE_REGRESSORS@$pnm_evlist@g; s@CONF_EVS@$conf_list@g; s@INPUT_T1_BRAIN@$input_brain@g" ${output}/${subject}/func/tmp/${task}design_pnmfeat.fsf\
	 > ${output}/${subject}/func/tmp/${subject}_${task}_design_pnm_feat.fsf \
	|| error_exit "cannot do feat PNM template for $subject $task"

	#open -e ${output}/${subject}/func/tmp/${subject}_${task}_design_pnm_feat.fsf
	feat ${output}/${subject}/func/tmp/${subject}_${task}_design_pnm_feat.fsf || error_exit "feat PNM failed for $subject $task"
}

run_pnm_feat_respiration_only () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local melodic_template=$5

	mkdir -p ${output}/${subject}/func/tmp
	cp -rf $melodic_template ${output}/${subject}/func/tmp/${task}design_pnmfeat.fsf


	 echo $subject $task "feat PNM"
	# #maybe add a check that FIX has run already correctly!
	 input_file=${output}/${subject}/func/${subject}_${task}_fixed.nii.gz
	 nvol_fmri=`fslval $input_file dim4`
	 output_melodic=${output}/${subject}/func/${subject}_${task}_feat_pnm_respiration
	 input_brain=${input}/${subject}/anat/${subject}_T1w_brain.nii.gz

	 pnm_evlist=${output}/${subject}/func/physio/${subject}_${task}_pnm/pnm_respiration_evlist.txt

	 conf_list=${output}/${subject}/func/${subject}_${task}_confound_M_PR_WM_CSF.txt

	 sed -e "s@OUTPUT_DIR@$output_melodic@g; s@NVAL@$nvol_fmri@g; s@INPUT_BOLD@$input_file@g; s@VOXELWISE_REGRESSORS@$pnm_evlist@g; s@CONF_EVS@$conf_list@g; s@INPUT_T1_BRAIN@$input_brain@g" ${output}/${subject}/func/tmp/${task}design_pnmfeat.fsf\
	  > ${output}/${subject}/func/tmp/${subject}_${task}_design_pnm_respiration_feat.fsf 
	# # \|| error_exit "cannot do feat PNM template for $subject $task"


	# #open -e ${output}/${subject}/func/tmp/${subject}_${task}_design_pnm_feat.fsf
	 feat ${output}/${subject}/func/tmp/${subject}_${task}_design_pnm_respiration_feat.fsf 
	 #|| error_exit "feat PNM failed for $subject $task"
}

pnm_noise_component () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4


	#save fixed fMRI data
	input_file=${output}/${subject}/func/${subject}_${task}_fixed.nii.gz
	#mv ${output}/${subject}/func/${subject}_${task}_melodic.ica/filtered_func_data_clean.nii.gz ${output}/${subject}/func/${subject}_${task}_fixed.nii.gz

	#get pnm signal and demean it
	pnm_residuals=${output}/${subject}/func/${subject}_${task}_feat_pnm_respphase.feat/stats/res4d.nii.gz

	pnm_signal=${output}/${subject}/func/${subject}_${task}_feat_pnm_respphase.feat/pnm_signal.nii.gz
	fslmaths $input_file -sub $pnm_residuals $pnm_signal
	fslmaths $pnm_signal -Tmean -mul -1 -add $pnm_signal $(dirname $pnm_signal)/$(basename $pnm_signal .nii.gz)_demeaned.nii.gz

	#clean pnm signal with FIX components
	mkdir -p ${output}/${subject}/func/tmp/${subject}_${task}_pnm.ica
	cp -rf ${output}/${subject}/func/${subject}_${task}_melodic.ica/ ${output}/${subject}/func/tmp/${subject}_${task}_pnm.ica
	cp $(dirname $pnm_signal)/$(basename $pnm_signal .nii.gz)_demeaned.nii.gz ${output}/${subject}/func/tmp/${subject}_${task}_pnm.ica/filtered_func_data.nii.gz

	/usr/local/fix/fix -a ${output}/${subject}/func/tmp/${subject}_${task}_pnm.ica/fix4melview_HCP_hp2000_thr30.txt

	fixed_pnm_signal=${output}/${subject}/func/tmp/${subject}_${task}_pnm.ica/filtered_func_data_clean.nii.gz
	fslmaths ${output}/${subject}/func/${subject}_${task}_fixed.nii.gz -sub $fixed_pnm_signal ${output}/${subject}/func/${subject}_${task}_fixedANDpnmed.nii.gz

}

pnm_noise_component_respiration () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4


	#save fixed fMRI data
	input_file=${output}/${subject}/func/${subject}_${task}_fixed.nii.gz
	#mv ${output}/${subject}/func/${subject}_${task}_melodic.ica/filtered_func_data_clean.nii.gz ${output}/${subject}/func/${subject}_${task}_fixed.nii.gz

	#get pnm signal and demean it
	pnm_residuals=${output}/${subject}/func/${subject}_${task}_feat_pnm_respiration.feat/stats/res4d.nii.gz

	pnm_signal=${output}/${subject}/func/${subject}_${task}_feat_pnm_respiration.feat/pnm_signal.nii.gz
	fslmaths $input_file -sub $pnm_residuals $pnm_signal
	fslmaths $pnm_signal -Tmean -mul -1 -add $pnm_signal $(dirname $pnm_signal)/$(basename $pnm_signal .nii.gz)_demeaned.nii.gz

	#clean pnm signal with FIX components
	mkdir -p ${output}/${subject}/func/tmp/${subject}_${task}_pnm_respiration.ica
	cp -rf ${output}/${subject}/func/${subject}_${task}_melodic.ica/ ${output}/${subject}/func/tmp/${subject}_${task}_pnm_respiration.ica
	cp $(dirname $pnm_signal)/$(basename $pnm_signal .nii.gz)_demeaned.nii.gz ${output}/${subject}/func/tmp/${subject}_${task}_pnm_respiration.ica/filtered_func_data.nii.gz

	/usr/local/fix/fix -a ${output}/${subject}/func/tmp/${subject}_${task}_pnm_respiration.ica/fix4melview_HCP_hp2000_thr30.txt

	fixed_pnm_signal=${output}/${subject}/func/tmp/${subject}_${task}_pnm_respiration.ica/filtered_func_data_clean.nii.gz
	fslmaths ${output}/${subject}/func/${subject}_${task}_fixed.nii.gz -sub $fixed_pnm_signal ${output}/${subject}/func/${subject}_${task}_fixedANDpnmed_respiration.nii.gz

}


slice_timing_correction () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4

	output_disk=$output
	output_disk=/Volumes/RAID_NP/Nikos/LOCUS1/Locus1_proc/BIDS_nii_structure/BIDS_proc

	input_file=${output_disk}/${subject}/func/${subject}_${task}_fixedANDpnmed_respiration.nii.gz

	slicetimer \
	--in=${input_file} \
	--out=$(dirname $input_file)/$(basename $input_file .nii.gz)_slicecorr.nii.gz \
	--tcustom=/Users/nikos/Locus1/BIDS_nii_structure/slicetimes_in_TR_fractions.txt \
	--repeat=2 
	#\|| error_exit "cannot slice correct $subject $task"

}



regressout_pnm_regressors () {

	local task=$1
	local subject=$2
	local input=$3
	local output=$4


	pnm_file=${output}/${subject}/func/physio/${subject}_${task}_pnm/pnm

	#read how many fmri volumes are there
	input_file=${output}/${subject}/func/${subject}_${task}_bold_topupcorrected.nii.gz
	vol=`fslval $input_file dim4`
	#make text file with 0 for each volume called blank.mat
	#rm ${pnm_file}_blank.mat
	#for (( num=1; num<=$vol; num++ )); do  echo 0; done >> ${pnm_file}_blank.mat
	#read evlist, replace /n with ,
	a=`cat ${pnm_file}_evlist.txt`
	conf_regr=`echo ${a// /,} | sed 's/ /,/g'`

	echo ${input_file}
	echo ${pnm_file}_blank.mat
	echo $(dirname $input_file)/$(basename $input_file .nii.gz)_bold_PNM_topupcorrected.nii.gz
	fsl_glm -i ${input_file} --out_res=$(dirname $input_file)/$(basename $input_file .nii.gz)_bold_PNM_topupcorrected.nii.gz -d ${pnm_file}_blank.mat

}

spatial_smoothing () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local sigma=$5

	input_file=${output}/${subject}/func/${subject}_${task}_fixedANDpnmed_respiration_slicecorr_cropped.nii.gz
	fslmaths $input_file -kernel gauss $sigma -fmean $(dirname $input_file)/$(basename $input_file .nii.gz)_smoothed.nii.gz

}

spatial_smoothing_with_susan () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local sigma=$5

	#find background value
	fslmaths ${output}/${subject}/func/${subject}_${task}_roi/epi_vol_pveseg.nii.gz -binv ${output}/${subject}/func/${subject}_${task}_roi/background.nii.gz
	background_intensity=`fslstats ${output}/${subject}/func/${subject}_${task}_roi/background.nii.gz -m`

	#find median value within brain
	fslmaths ${output}/${subject}/func/${subject}_${task}_roi/epi_vol_pveseg.nii.gz -bin ${output}/${subject}/func/${subject}_${task}_roi/brain_all.nii.gz
	brain_intensity=`fslstats ${output}/${subject}/func/${subject}_${task}_roi/brain_all.nii.gz -m`

	#calculate 0.75 of contrast as a threshold
	thresh=0.75
	ddd=`bc <<< "$brain_intensity-$background_intensity"`
	thresh_up=`bc <<< "$ddd*$thresh"`

	#call susan
	output_disk=/Volumes/RAID_NP/Nikos/LOCUS1/Locus1_proc/BIDS_nii_structure/BIDS_proc
	input_file=${output_disk}/${subject}/func/${subject}_${task}_fixedANDpnmed_slicecorr_cropped.nii.gz

	fslmaths $input_file -Tmean ${output}/${subject}/func/${subject}_${task}_roi/mean_func.nii.gz
	mean_func=${output}/${subject}/func/${subject}_${task}_roi/mean_func.nii.gz
	/usr/local/fsl/bin/susan $input_file $thresh_up $sigma 3 1 1 $mean_func $thresh_up $(dirname $input_file)/$(basename $input_file .nii.gz)_smoothedsusan_15mm.nii.gz

	
}


make_GMandWM_mask () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4

	output_folder=${output}/${subject}/func/${subject}_${task}_roi
	
	epi_vol=${output_folder}/epi_vol.nii.gz
	fast $epi_vol 
	fslmaths ${output_folder}/epi_vol_pve_1.nii.gz -add ${output_folder}/epi_vol_pve_2.nii.gz -bin ${output_folder}/GMWM_epi_mask.nii.gz
	gunzip -k ${output_folder}/GMWM_epi_mask.nii.gz
}

make_The_Plot () {
	local task=$1
	local subject=$2
	local input=$3
	local output=$4
	local FS_fil=$5

	#this is The Plot, as described by Jonathan Power. Made a few modifications to work with BIDs and include brainstem, etc 

	output_folder=${output}/${subject}/func/${subject}_${task}_roi
	T1_image=${input}/${subject}/anat/${subject}_T1w.nii.gz
	directory_demo=${output_folder}/plot
	fmri_dataset=${output}/${subject}/func/${subject}_${task}_fixedANDpnmed_slicecorr.nii.gz

	t1_2_epi_mat=${output_folder}/struct2epi.mat
	masterimg=${directory_demo}/data/mprage_at_EPI.nii.gz


	mkdir -p $directory_demo
	mkdir -p ${directory_demo}/data


	mri_convert ${FS_fil}/${subject}/mri/aparc+aseg.mgz ${directory_demo}/data/aparc+aseg.atlas.nii.gz
	cp $fmri_dataset ${directory_demo}/data/bold1_at_EPI.nii.gz
	flirt -interp sinc -in $T1_image -ref $fmri_dataset -applyxfm -init $t1_2_epi_mat -out ${directory_demo}/data/mprage_at_EPI.nii.gz



	echo ----------------------------------------------------------
	echo ----------------------------------------------------------
	echo CLEARING OUT OLD MASKS IF PRESENT
	echo ----------------------------------------------------------
	echo ----------------------------------------------------------

	rm -rf `ls ${directory_demo}/data/*_ero*.nii.gz`

	echo ----------------------------------------------------------
	echo ----------------------------------------------------------
	echo CREATE MASKS FROM FREESURFER SEGMENTATION
	echo ----------------------------------------------------------
	echo ----------------------------------------------------------

	cd ${directory_demo}/data
	cp ${output}/${subject}/func/tmp/motion_outliers/${subject}_${task}_bold_FDvalues.txt fd.txt




	# extract the GM, WM, CSF, and WB compartments
	# everything labeled in FS, followed by resampling to the BOLD resolution
	3dcalc -a aparc+aseg.atlas.nii.gz \
	-expr 'not(equals(a,0))' \
	-prefix aparc+aseg.atlas.INBRAINMASK_ero0.nii.gz

	# the major WM compartments, with 4 erosions at the T1 resolution followed by resampling to the BOLD resolution
	3dcalc -a aparc+aseg.atlas.nii.gz \
	-expr 'equals(a,2)+equals(a,7)+equals(a,41)+equals(a,46)+equals(a,251)+equals(a,252)+equals(a,253)+equals(a,254)+equals(a,255)' \
	-prefix aparc+aseg.atlas.WMMASK_ero0.nii.gz

	3dcalc -a aparc+aseg.atlas.WMMASK_ero0.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.WMMASK_ero1.nii.gz

	3dcalc -a aparc+aseg.atlas.WMMASK_ero1.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.WMMASK_ero2.nii.gz

	3dcalc -a aparc+aseg.atlas.WMMASK_ero2.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.WMMASK_ero3.nii.gz

	3dcalc -a aparc+aseg.atlas.WMMASK_ero3.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.WMMASK_ero4.nii.gz

	# the major CSF compartments, with 4 erosions at the T1 resolution followed by resampling to the BOLD resolution
	3dcalc -a aparc+aseg.atlas.nii.gz \
	-expr 'equals(a,4)+equals(a,43)+equals(a,14)' \
	-prefix aparc+aseg.atlas.CSFMASK_ero0.nii.gz

	3dcalc -a aparc+aseg.atlas.CSFMASK_ero0.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.CSFMASK_ero1.nii.gz

	3dcalc -a aparc+aseg.atlas.CSFMASK_ero1.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.CSFMASK_ero2.nii.gz

	3dcalc -a aparc+aseg.atlas.CSFMASK_ero2.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.CSFMASK_ero3.nii.gz

	3dcalc -a aparc+aseg.atlas.CSFMASK_ero3.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.CSFMASK_ero4.nii.gz

	# the gray matter ribbon (amygdala and hippocampus need to be added - 17 18 53 54
	3dcalc -a aparc+aseg.atlas.nii.gz \
	-expr 'within(a,1000,3000)+equals(a,17)+equals(a,18)+equals(a,53)+equals(a,54)' \
	-prefix aparc+aseg.atlas.GM_RIBBONMASK_ero0.nii.gz

	# the cerebellum
	3dcalc -a aparc+aseg.atlas.nii.gz \
	-expr 'equals(a,47)+equals(a,8)' \
	-prefix aparc+aseg.atlas.GM_CBLMMASK_ero0.nii.gz

	3dcalc -a aparc+aseg.atlas.GM_CBLMMASK_ero0.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.GM_CBLMMASK_ero1.nii.gz

	3dcalc -a aparc+aseg.atlas.GM_CBLMMASK_ero1.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.GM_CBLMMASK_ero2.nii.gz

	# the subcortical nuclei
	3dcalc -a aparc+aseg.atlas.nii.gz \
	-expr 'equals(a,11)+equals(a,12)+equals(a,10)+equals(a,49)+equals(a,50)+equals(a,51)' \
	-prefix aparc+aseg.atlas.GM_SCMASK_ero0.nii.gz

	3dcalc -a aparc+aseg.atlas.GM_SCMASK_ero0.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.GM_SCMASK_ero1.nii.gz

	3dcalc -a aparc+aseg.atlas.GM_SCMASK_ero1.nii.gz -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
	-prefix aparc+aseg.atlas.GM_SCMASK_ero2.nii.gz

	# all gray matter
	3dcalc -a aparc+aseg.atlas.nii.gz \
	-expr 'within(a,1000,3000)+equals(a,17)+equals(a,18)+equals(a,53)+equals(a,54)+equals(a,47)+equals(a,8)+equals(a,11)+equals(a,12)+equals(a,10)+equals(a,49)+equals(a,50)+equals(a,51)' \
	-prefix aparc+aseg.atlas.GM_ALLMASK_ero0.nii.gz


	#transform masks to EPI space and binarize
	for i in aparc+aseg.atlas.?*.nii.gz; do 
		flirt -interp sinc -in $i -ref $fmri_dataset -applyxfm -init $t1_2_epi_mat -out $(basename $i .nii.gz)_EPI.nii.gz
		fslmaths $i -thr 0.5 -bin $i
	done	

	echo ----------------------------------------------------------
	echo ----------------------------------------------------------
	echo ALL DONE MAKING MASKS
	echo ----------------------------------------------------------
	echo ----------------------------------------------------------

	cd .. 
	matlab -nodesktop -nosplash -r "try plotdemo(); catch; end; quit"

	cd ~


}

#	----------------------------------------------------------------
#	parameters  setup
#		we assume that within the dir_path there are a bunch of subjects files 
#		and within them another folder with unpacked dcm data
# 		
#	----------------------------------------------------------------

BIDS_nii_folder=/Users/nikos/Locus1/BIDS_nii_structure
BIDS_nii_folder=/Volumes/RAID_NP/Nikos/LOCUS1/Locus1_proc/BIDS_nii_structure
FS_files=/Users/nikos/Locus1/BIDS_nii_structure/FS/FreeSurferOutputv6b

subject_list=/Users/nikos/Locus1/BIDS_nii_structure/subject_list_complete.txt

acq_param=${BIDS_nii_folder}/acq_params.txt
melodic_template=${BIDS_nii_folder}/template_melodic.fsf
output_file=${BIDS_nii_folder}/BIDS_proc
feat_pnm_enc_template=${BIDS_nii_folder}/template_pnmANDconf_regression.fsf

fmri_array=( ret enc rs1 rs2 )
task_array=( ret enc )
rs_array=( rs1 rs2 )

#	----------------------------------------------------------------
#	commands
#	----------------------------------------------------------------

#make processing file
mkdir -p $output_file || error_exit "cannot make output file"

echo "start from:"
read lower_limit

echo "stop at:"
read upper_limit
n=1;

#go through all subjects in list
cat $subject_list | while read line
do
	n=$((n+1))
	if [[ $n -lt $lower_limit ]]; then
		continue;
	elif [[ $n -gt $upper_limit ]]; then
		break;
	fi
	echo $line
	echo $$
   #motion correction
    echo "Motion correction"
    for task in "${fmri_array[@]}"; do motion_correction $task $line $BIDS_nii_folder $output_file & done || error_exit "cannot do motion correction for $line"
    wait ${!}
   echo "motion outliers"
   #motion outliers for task
   for task in "${task_array[@]}"; do motion_outliers $task $line $BIDS_nii_folder $output_file 0.9 & done || error_exit "cannot find motion outliers for $line"
   wait ${!}
   echo "topup"
   for task in "${fmri_array[@]}"; do topup_apply $task $line $BIDS_nii_folder $output_file $acq_param & done || error_exit "cannot do topup for $line"
   wait ${!}
  
   echo "melodic"
   for task in "${rs_array[@]}"; do run_melodic $task $line $BIDS_nii_folder $output_file $melodic_template & done || error_exit "cannot do melodic for $line"  
   wait ${!}
   echo $$
   for task in "${rs_array[@]}"; do run_fix $task $line $BIDS_nii_folder $output_file & done || error_exit "cannot do fix for $line"  
   wait ${!}


   # #make physiological regressors
   wait ${!}
   for task in "${fmri_array[@]}"; do pnm_regressors $task $line $BIDS_nii_folder $output_file & done 
   || error_exit "cannot make pnm regressors for $line"  

	#run feat to find PNM residuals
    for task in "${fmri_array[@]}"; do run_pnm_feat $task $line $BIDS_nii_folder $output_file $feat_pnm_enc_template & done 
    wait ${!} 
 	for task in "${fmri_array[@]}"; do run_pnm_feat_respiration_only $task $line $BIDS_nii_folder $output_file $feat_pnm_enc_template & done 
    wait ${!} 
    for task in "${fmri_array[@]}"; do pnm_noise_component $task $line $BIDS_nii_folder $output_file & done 
    wait ${!} 
   	for task in "${fmri_array[@]}"; do pnm_noise_component_respiration $task $line $BIDS_nii_folder $output_file & done 
    wait ${!} 
    
   #echo "slice time correction $line"
   for task in "${fmri_array[@]}"; do slice_timing_correction $task $line $BIDS_nii_folder $output_file & done 
   wait ${!}
 
   #divide your FWHM (in mm) by 2.3548 to get sigma
   #here set at smoothing of 1.5mm FWHM
   sigma_smoothing=0.6369
   for task in "${fmri_array[@]}"; do spatial_smoothing $task $line $BIDS_nii_folder $output_file $sigma_smoothing & done 
   #for task in "${fmri_array[@]}"; do spatial_smoothing_with_susan $task $line $BIDS_nii_folder $output_file $sigma_smoothing & done 

   wait ${!}    
   #make a GM and WM mask. This can be useful 1. as a preset mask to ensure that SPM 1st level analysis doesn't throw out relevant areas 
   #2. to constrain analysis to plausible brain regions
   for task in "${fmri_array[@]}"; do  make_GMandWM_mask $task $line $BIDS_nii_folder $output_file & done
 	wait ${!} 
 	#task="rs1" 
 	for task in "${fmri_array[@]}"; do make_The_Plot $task $line $BIDS_nii_folder $output_file $FS_files & done
  	wait ${!}  
done



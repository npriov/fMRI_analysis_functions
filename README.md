A set of functions for fMRI analysis. Randomly written in shell, R and matlab. Will expand in due time

The BIDS_fmri_preproc.sh is the preprocessing pipeline I had set up for our 7T BOLD fMRI datasets. It assumes a BIDS filing setup (more about the BIDS specifications here http://bids.neuroimaging.io/). The pipeline is based around FSL and FreeSurfer functions. Kudos also to Jonathan Power for sharing his Plot with us; I have added it here, since its quite useful to evaluate pipelines.

I am quite interested in characterising the physiological noise in fMRI datasets, so a big part of the pipeline is focused on this. Within the pipeline we
1. Apply a 6dof motion correction
2. Extract motion outliers and save it in regressor form (this can probably be combined with the first step; haven’t bothered to change it yet)
3. Distort correct, assuming a reversed phase encoding acquisition
4. Run Melodic and FIX
5. Apply PNM (in two flavours, respiration only and cardiac and respiration together)
6. Slice timing correct
7. Smooth
8. Make a GM/WM mask from the EPI data themselves
9. Make The Plot

Crucially, the PNM regression is combined with the FIX regression so that we do not regress the same noise time series twice; this is important to not reintroduce noise to our data (kudos to the FSL people in Oxford for coming up with this one). The PNM also includes a regressor of the respiration phase only (i.e. binarised expiration/inspiration; calculated in get_resp_phase.m with a simple Hilbert transform).
The slice timing correction is performed assuming a multi band acquisition (acq table extracted from the DICOM and scaled to the TR).
As pointed out by @srirangakashyap, we can also reduce the 2-interpolations done in motion correction and distortion correction to a single one, with the awesome ANTs package (highly recommended for quality registrations as well).



I also recently wanted to take a look at the cross-correlation variance along different lags (to examine the need to adapt this for different populations/regions). I played with parallelising this in R, since R has a whole bunch of time series analysis packages and analysis with R seems to be gaining traction in neuroimaging circles. Loading the niftis themselves takes an obnoxious length of time; though I think there are some newer nifti reading functions that are faster. Anyway, take a look if interested in the lag_analysis.r 

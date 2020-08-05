Preparation
-----------

+sort_MIDIFSUA
  -put *all* MIDI/FSUA files in one directory, e.g. /media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_P90RATZKA/
  -MIDI files can be also Lab calibrations and will be in a tmp folder which can be deleted

OBSOLETE

  -create directory /home/amueller/work/MIDI_FSUA/MIDI_FSUA_COMM<XX>
  -create directory /media/disk_MIDIFSU/MIDI_FSUA/<XX>
  -create directory /media/disk_MIDIFSU/MIDI_FSUA/COMM<XX>/FSUAdata  +FringeScan  +LabCalibration  +not_useable  +SkyCalibration
  -create directory /media/disk_MIDIFSU/MIDI_FSUA/COMM<XX>/MIDIdata  +Acquisition_Images  +not_useable  +Photometry
  -sort data with respect to the directories
  -make sure that MIDI...99.fits are moved in not_useable

dfits MIDI.2014-04-13T0*01.fits | fitsort dpr.type obs.targ.name dpr.catg date
dfits Photometry/MIDI.2014-04-13T0*01.fits | fitsort dpr.type obs.targ.name dpr.catg date
dfits FSUA_OBS*_00* | fitsort dpr.type obs.targ.name dpr.catg date
dfits SkyCalibration/FSUA_SKY_DARK*_00* | fitsort dpr.type obs.targ.name dpr.catg date
dfits SkyCalibration/FSUA_SKY_FLAT*_00* | fitsort dpr.type obs.targ.name dpr.catg date

+observation.txt
  -enter the necessary data

+observation_FSUA.txt
  -enter the necessary data

+MIDI_NULL_and_replace_OPD_values_ZERO
  -!!!!!   IF A FILE HAS MORE THAN ONE TIME STAMP REDUCTION PROBABLY FAILS. MAYBE USE ORIGINAL DATA IF FSU REDUCTION NOT PERFORMED
  -looks for MIDI NULL entries and replaces them
  -OPD values are replaced with 0
  -backup is created

+in Calibration: calibrate_FSUA_LAB, 'B'
		  calibrate_FSUA_SKY...., /nodark
		OR use calibrate_all_scans and call writeresult.pro


+FSUA_replace_GD_PD_values
  -recompute GD and PD values with LabCal and replace the original values
  -backup is created

+merge_MIDI_FSU_observations
  -MERGE ONLY AFTER NULL TIME STAMP SEARCH
  -works only for two observations for the moment and only for GEN_RECORD files
  -OBSERVATION.TXT and OBSERVATION_FSUA.TXT have to be modified!
  -in case there are N>1 files recorded merge them into one
  -files to merge have to be set up manually

+merge_MIDI_Photometry: if it was decided to record another photometry

+MIDI_MakeCustomMask
  -only the first 3 files are used, otherwise problems with RAM

  copy checkmasks.pro in the corresponding run directory

OBSOLETE

  +create mask
    start mia in the data folder, select fringe track file
    f = midigui()
    midiMakeMask, 'star', f
    copy files in /opt/MIDI_FSUA/MIDI_FSUA_COMM<XX>

  +MIA
    files = midigui()
    mask=miamask(files)
    mask->gui
    mask->write_mask
    maskname=mask->get_mask_name()

  +create mask from calibrator and copy it to /opt/MIDI_FSUA/MIDI_FSUA_COMM'<commnum>
    -reduceMIDI_photometry


Conditions
----------

+get_V_projbaseline
  -create stars_COMMXX.txt containing
    id,ra,dec,vmag,jmag,hmag,kmag,f10,diam
  -creates Vis_projBaseline_COMXX.txt and projBaseline_COMMXX.txt

+get_atmospheric_paramameters
  -creates conditions_COMMXX.txt

+FSU_FringeTrackPerformance / FSU_FringeTrackPerformance_wrappedPD
  -!!! can only be executes if phase unwrapping was performed!!! calcMIDI_Koresko
  -creates FSUA_lockratio_COMMXX.txt

+compute_vsibility_binary
  -computes the visibility of a binray
  -values (coordinates, files, fluxes) have to be provided manually

Reduction
---------

+calcMIDI_Koresko (for Custom- and StandardMask)
  -predicts N band GD, PD and Dispersion using PRIMA FSUA

+reduceMIDI_EWS_CustomMask / reduceMIDI_EWS_StandardMask
  reduces data in standard EWS fashion but using the individual masks

+reduceMIDI_CohInt_CustomMask / reduceMIDI_CohInt_StandardMask

+reduceMIDI_Koresko_CustomMask / reduceMIDI_Koresko_StandardMask -- run the calibrators first!

+reduceMIDI_faintEWS_CustomMask

+MIDI_number_of_good_frames_CustomMask / MIDI_number_of_good_frames_StandardMask
  get number of used frames for MIDI reduction


Calibration
-----------

+MIDI_Calibrator_TF_Fcorr_CustomMask / (MIDI_Calibrate_Fcorr_StandardMask does not exist yet)
  -has to be executed from MIA+EWS
  -measures the TF of each night
  -Only caliibrator stars of course
  -make modifications similar to MIDI_Calibrator_TF_Fcorr_CustomMask_P92CHOQUET w.r.t. e.g. errors and plots

OR in case of WISE photometry
+MIDI_Calibrator_TF_Fcorr_CustomMask_P92CHOQUET
  -in observation.txt MANUALLY SELECT/deselect observation
  -diameters, used baseline have to be provided manually

+MIDI_Calibrator_TF_Fcorr_CustomMask_Photometry
  -create TF from Calibrator observations if good photometry is available
  -diameters, used baseline have to be provided manually
  -in observation.txt MANUALLY SELECT/deselect observation

---

CHOQUET
+calibrateVis_Koresko_CustomMask_P92CHOQUET
  claibrate visibilities using binned data (corr, photometry)
+calibrate_Fcorr_EWS_CustomMask_P92Choquet
  using WISE photometry to get calibrated correlated fluxes

---

also for CHOQUET
+MIDI_calibrate_Fcorr_CustomMask
  -in observation.txt MANUALLY SELECT/deselect observation for night and baseline
  -calibrates rawFcorr, Calibrator and Science

+MIDI_Calibrator_calibrate_Vis_CustomMask
  -computes a calibrated Visibility of Calibrator Stars using the total photometry of the vanBoekel database
  -HAS TO BE EXECUTED FROM MIA+EWS!


Photometry
----------

  +reduceMIDI_KoreskoPhotometry_CustomMask

CHOQUET
  +average_faint_photometry_CustomMask
    DON'T FORGET TO COPY ORIGINAL DATA BACK!
    BINNING
    average the data of same target observations in wavelength and time
    plots photometry of the same target acquired over several exposures of the same run

  +in case photometry is useless, try WISE and 2MASS photometry
    -was done with P92CHOQUET
    -go there, execute totflux.pro after editing the data_***.txt files with the corresponding magnitudes
    -inside the routine proivde parallax and stellar diameter
    -MIDIreduce_EWS_CustomMask has to exist to get the wavelength


Don't care about photometry use vanBoekel Database
  +reduceMIDI_EWSPhotometry_StandardMask
    -extract photometry if present and get instrumental visibility (.redcal)

  +reduceMIDI_KoreskoPhotometry_StandardMask
    -extract photometry if present and get instrumental visibility (.redcal)



Data Quality
------------

+MIDI_plot_OPD_Delay_CustomMask / MIDI_plot_OPD_Delay_StandardMask
  NO TWOPASS
  ORIGversion plots things a bit differently and has the twopass option

+check_MIDI_data_qual
  -has to be modified for each individual file

+check_MIDI_AcqImage


Plots
-----

(under DataQuality)

+plot_fcorr_phi
  -plots of (un)calibrated diff phi and Fcorr

+plot_calFcorr_EWS_Koresko_CustomMask
  -create plots with EWS,Koresko calFcorr

+plot_MIDI_FSUA_MeasuredPredicted_CustomMask
  plot MIDI GD found by EWS and computed by Koresko
  makes only sense if target was bright enough for EWS



----------------------------------------------------------------------------------------------


+DO NOT USE - reduceMIDI_usingFSUAGD
  -reduces MIDI+FSUA data
  -FSUA GD used and written in MIDI file instead of MIDI GD
  -requires analyze_FSUAGD_for_MIDI_reduction
  -requires the previous results from the reduction by EWS

+DO NOT USE - MIDI_NULL_and_replace_OPD_values
  -looks for MIDI NULL entries and replaces them
  -OPD values are used from the FSUA data
  -original MIDI FT files are copied for backup


plot_Fcorr_IRAS_allstars_CohInt
analyze_FSUAGD_for_MIDI_reduction_GDwin
analyze_get_fsu_gd_orig
analyze_get_fsu_gd_orig_GDwin
analyze_MIDI_FSUA_GD_COMM9
Fcorr_plot
get_fsu_gd_orig
get_fsu_gd_orig_GDwin
N_Vis_Com9	;computes Visibility
plot_Fcorr_allstars_COMM9_GRISM
plot_Fcorr_allstars_COMM9_PRISM
plot_FSU_FUOFFSET_COMM9_GRISM
plot_FSU_FUOFFSET_COMM9_PRISM
plot_FSU_RTOFFSET_COMM9_GRISM
plot_FSU_RTOFFSET_COMM9_PRISM
plot_FSU_UWPhase_COMM9_GRISM
plot_FSU_UWPhase_COMM9_PRISM
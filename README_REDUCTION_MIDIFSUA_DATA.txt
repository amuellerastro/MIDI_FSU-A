mask=miamask(files)
mask->gui
mask->write_mask
maskname=mask->get_mask_name()

valid for files that were taken while FSU-A fringe tracking

1.) steps to perform reduction of bright (>10Jy) MIDI targets that were taken using FSU-A as fringe tracker
can be executed with reduceMIDI_EWS_COMMX.pro

oir1dCompressData 'exposures' <Mask> <out_compressed>
oirFormFringes <compressed> <out_firnges>  / for sources with F<7Jy '-removeAverage' option is applied
oirRotateInsOpd.noOPD <fringes> <out_insopd>
oirGrouDelay.noOPD <insopd> <out_groupdelay>		there is no affect if the line with OPD is comment out or not
oirRotateGroupDelay.noOPD <fringes> <groupdelay> <out_ungroupdelay>
oirAutoFlag.noOPD <ungroupdelay> <groupdelay> <out_flag>
oirAverageVis <ungroupdelay> <flag> <out_corr>

all files are moved to 'MIDIreduced_EWS'


2.) steps to perform reduction of faint (<10Jy) MIDI targets that were taken using FSU-A as fringe tracker, coherent integration
can be executed with reduceMIDI_cohint_COMMX.pro

oir1dCompressData 'exposures' <Mask> <out_compressed>
oirFormFringes <compressed> <out_firnges>  / for sources with F<7Jy '-removeAverage' option is applied
oirRotateInsOpd.noOPD <fringes> <out_insopd>
analyze_FSUAGD_for_MIDI_reduction.pro - produces <.time> and <.flags> and <.midigd_calc_fsugd> and <.fsugd> / all ascii files
oirAverageVis.cohint <insopd> <.fsuflags> <out_corr>

all files are moved to 'MIDIreduced_CohInt'


3.) steps to perform reduction of MIDI targets that were taken using FSU-A as fringe tracker, FSU-A GD only is used
can be executed with reduceMIDI_usingFSUAGD_COMMX.pro

oir1dCompressData 'exposures' <Mask> <out_compressed>
oirFormFringes <compressed> <out_firnges>  / for sources with F<7Jy '-removeAverage' option is applied
oirRotateInsOpd.noOPD <fringes> <out_insopd>
analyze_FSUAGD_for_MIDI_reduction.pro - produces <.time> and <.flags> and <.fsugd> and <.fsugd> / all ascii files
oirRotateGroupDelay.cohint <fringes> <.fsugd> <out_ungroupdelay>
oirAverageVis.cohint <ungroupdelay> <.fsuflags> <out_corr>

all files are moved to 'MIDIreduced_usingFSUAGD'


4.) steps to perform reduction of MIDI targets that were taken using FSU-A as fringe tracker, FSU-A GD is used to compute MIDI GD
can be executed with reduceMIDI_calcMIDIGD_COMMX.pro

oir1dCompressData 'exposures' <Mask> <out_compressed>
oirFormFringes <compressed> <out_firnges>  / for sources with F<7Jy '-removeAverage' option is applied
oirRotateInsOpd.noOPD <fringes> <out_insopd>
analyze_FSUAGD_for_MIDI_reduction.pro - produces <.time> and <.flags> and <.calc_midigd> and <.fsugd> / all ascii files
oirRotateGroupDelay.cohint <fringes> <.calc_midigd> <out_ungroupdelay>
oirAverageVis.cohint <ungroupdelay> <.fsuflags> <out_corr>

all files are moved to 'MIDIreduced_calcMIDIGD'


5.) Misc. useable for both

get_atmospheric_paramameters: retrieve t_c, seeing, airmass saved in conditions.txt
get_fsu_lockratio: retrieve estimated lockratio for FSUA
analyze_MIDI_FSUA_GD_COMMX: for bright sources where GD by EWS could be measured, gives theoretical TF, MIDI-GD vs FSUA-GD and linear fits and CCF
plot_theoTF: plots theoretical TF in one graph, output from analyze_MIDI_FSUA_GD_COMMX
plot_Fcorr_allstars/COMX: produces a plot Fcorr vs wavelength for all used reduction methods (superimposed)
plot_Fcorr_IRAS_allstars: produces a plot Fcorr vs IRAS flux for all different used reduction methods (single plots)
get_V_projbaseline: computes from MIDI files the projectes baseline and V^2 if diameter given, else V^2 set to 1.00
check_MIDI_data_qual: quality check on MIDI data, i.e., delay functions
N_Vis_ComX: computes photometry and instrumental visibility
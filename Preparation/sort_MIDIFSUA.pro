;put *all* MIDI/FSUA files in one directory, e.g. /media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_P90RATZKA/
;MIDI files can be also Lab calibrations, will be deleted

pro sort_MIDIFSUA

base = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/'

obsnum = ''
read, 'Enter Observing Number (e.g. P90RATZKA): ', obsnum
dir = base+'MIDI_FSUA_'+obsnum+'/'
;dir = base+'MIDI_FSUA_'+obsnum+'/'

;create FSUA directories
fsudir1 = dir+'FSUAdata/'
fsudir2 = fsudir1+'FringeSearch/'
fsudir3 = fsudir1+'FringeScan/'
fsudir4 = fsudir1+'LabCalibration/'
fsudir5 = fsudir1+'SkyCalibration/'
fsudir6 = fsudir1+'not_useable/'

file_mkdir, fsudir1
file_mkdir, fsudir2
file_mkdir, fsudir3
file_mkdir, fsudir4
file_mkdir, fsudir5
file_mkdir, fsudir6

;create MIDI directories
mididir1 = dir+'MIDIdata/'
mididir2 = mididir1+'Acquisition_Images/'
mididir3 = mididir1+'Photometry/'
mididir4 = mididir1+'not_useable/'
mididir5 = mididir1+'tmpcalib/'

file_mkdir, mididir1
file_mkdir, mididir2
file_mkdir, mididir3
file_mkdir, mididir4
file_mkdir, mididir5


;sort FSUA files
file = ''
file = file_search(dir+'FSUA_FIRST_FRINGE_SCAN_*.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+fsudir2+'.'
file = ''
file = file_search(dir+'FSUA_SECOND_FRINGE_SCAN_*.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+fsudir3+'.'
file = ''
file = file_search(dir+'FSUA_OBS_*.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+fsudir1+'.'
file = ''
file = file_search(dir+'FSUA_SKY_*.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+fsudir5+'.'
file = ''
file = file_search(dir+'PACMAN_LAB_*.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+fsudir4+'.'


;sort MIDI files

file = ''
file = file_search(dir+'MIDI*_99.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+mididir5+'.'

file = ''
file = file_search(dir+'MIDI*.fits', count=nfiles)

for i=0,nfiles-1 do begin

  dum = mrdfits(file[i], 0, hdr, /silent)
  obsname = get_eso_keyword(hdr, 'HIERARCH ESO OBS NAME')
  dprtype = get_eso_keyword(hdr, 'HIERARCH ESO DPR TYPE')
  tplid = get_eso_keyword(hdr, 'HIERARCH ESO TPL ID')

;Calib files
  if (obsname eq 'Start-up-test' and tplid eq 'MIDI_autotest_tec_startup') then $
    spawn, 'mv '+file[i]+' '+mididir5+'.'
  if (obsname eq 'Autotest-complete' and tplid eq 'MIDI_autotest_tec_detron') then $
    spawn, 'mv '+file[i]+' '+mididir5+'.'
  if (obsname eq 'Autotest-complete' and tplid eq 'MIDI_autotest_tec_detlin') then $
    spawn, 'mv '+file[i]+' '+mididir5+'.'
  if (obsname eq 'Autotest-complete' and tplid eq 'MIDI_autotest_tec_wavecal') then $
    spawn, 'mv '+file[i]+' '+mididir5+'.'
  if (obsname eq 'Autotest-complete' and tplid eq 'MIDI_autotest_tec_dsptrn') then $
    spawn, 'mv '+file[i]+' '+mididir5+'.'
  if (obsname eq 'Autotest-complete' and tplid eq 'MIDI_autotest_tec_refpix') then $
    spawn, 'mv '+file[i]+' '+mididir5+'.'

;Aquisition
  if (dprtype eq 'COARSE,OBJECT' and tplid eq 'MIDI_midi+fsu_acq') then $
    spawn, 'mv '+file[i]+' '+mididir2+'.'

;Observations
  if (dprtype eq 'TRACK,OBJECT,DISPERSED' and tplid eq 'MIDI_midi+fsu_obs') then $
    spawn, 'mv '+file[i]+' '+mididir1+'.'

;Photometry
  if (dprtype eq 'PHOTOMETRY,OBJECT' and tplid eq 'MIDI_midi+fsu_obs') then $
    spawn, 'mv '+file[i]+' '+mididir3+'.'

endfor



stop
end


; OBJECT  = 'COARSE,OBJECT'      / Original target
; HIERARCH ESO OBS NAME        = 'CAL_HD218594' / OB name
; ERARCH ESO TPL NAME        = 'MIDI MIDI+FSU ACQ' / Template name
; HIERARCH ESO TPL ID          = 'MIDI_midi+fsu_acq' / Template signature ID
; HIERARCH ESO DPR TYPE        = 'COARSE,OBJECT' / Observation type

; FTK
; HIERARCH ESO TPL ID          = 'MIDI_midi+fsu_obs' / Template signature ID
; HIERARCH ESO TPL NAME        = 'MIDI MIDI+FSU OBS' / Template name
; HIERARCH ESO DPR TYPE        = 'TRACK,OBJECT,DISPERSED' / Observation type
; 
; PHOT
; HIERARCH ESO TPL ID          = 'MIDI_midi+fsu_obs' / Template signature ID
; HIERARCH ESO TPL NAME        = 'MIDI MIDI+FSU OBS' / Template name
; HIERARCH ESO DPR TYPE        = 'PHOTOMETRY,OBJECT' / Observation type

; Calib1
; HIERARCH ESO OBS NAME        = 'Start-up-test' / OB name
; HIERARCH ESO TPL ID          = 'MIDI_autotest_tec_startup' / Template signature
; HIERARCH ESO TPL NAME        = 'MIDI AUTOTEST TEC STARTUP' / Template name
; HIERARCH ESO DPR TYPE        = 'IMAGING '   / Observation type

; Calib2
; HIERARCH ESO OBS NAME        = 'Autotest-complete' / OB name
; HIERARCH ESO DPR TYPE        = 'BIAS    '   / Observation type
; HIERARCH ESO TPL ID          = 'MIDI_autotest_tec_detron' / Template signature I
; HIERARCH ESO TPL NAME        = 'MIDI AUTOTEST TEC DETRON' / Template name

; Calib3
; HIERARCH ESO TPL NAME        = 'MIDI AUTOTEST TEC DETLIN' / Template name
; HIERARCH ESO TPL ID          = 'MIDI_autotest_tec_detlin' / Template signat
; HIERARCH ESO DPR TYPE        = 'FLAT    '   / Observation type
; HIERARCH ESO OBS NAME        = 'Autotest-complete' / OB name

; Calib4
; HIERARCH ESO OBS NAME        = 'Autotest-complete' / OB name
; HIERARCH ESO DPR TYPE        = 'WAVE,SPECTEMPL' / Observation type
; HIERARCH ESO TPL NAME        = 'MIDI AUTOTEST TEC WAVECAL' / Templ
; AcqHIERARCH ESO TPL ID          = 'MIDI_autotest_tec_wavecal' / Template si

; calib5
; HIERARCH ESO DPR TYPE        = 'WAVE    '   / Observation type
; HIERARCH ESO OBS NAME        = 'Autotest-complete' / OB name
; HIERARCH ESO TPL NAME        = 'MIDI AUTOTEST TEC DSPTRN' / Template name
; HIERARCH ESO TPL ID          = 'MIDI_autotest_tec_dsptrn' / Template signature

; calib6
; HIERARCH ESO OBS NAME        = 'Autotest-complete' / OB name
; HIERARCH ESO DPR TYPE        = 'FMTCHCK '   / Observation type
; HIERARCH ESO TPL ID          = 'MIDI_autotest_tec_refpix' / Template sign
; HIERARCH ESO TPL NAME        = 'MIDI AUTOTEST TEC REFPIX' / Template name
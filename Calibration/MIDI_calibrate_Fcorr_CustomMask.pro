@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@julday.pro
@caldat.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@interpol.pro


pro MIDI_calibrate_Fcorr_CustomMask

;**********************************************************************************************


obsnum = ''
; read, 'Enter Observation ID: ', obsnum
; obsnum = 'P92CHOQUET'
obsnum = 'P92CHOQUETPOTT'

;==========================================

; if (obsnum eq 'P92CHOQUET') then begin
if (obsnum eq 'P92CHOQUETPOTT') then begin

  base = ['U1U3','U1U4']	;used baselines for this program

endif

;==========================================

basequest = ''
for i=0,n_elements(base)-1 do print, strcompress(i,/rem), '  ', base[i]
read, 'Which baseline?: ', basequest

readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', skipline=1, /silent

idx = where(comm eq obsnum)
if (idx[0] ne -1) then begin

  comm = comm[idx]
  id = id[idx]
  scical = scical[idx]
  nfl = nfl[idx]
  midi = midi[idx]
  fsu = fsu[idx]
  maskname = maskname[idx]
  photA = photA[idx]
  photB = photB[idx]

endif else begin

  print, ''
  print, 'Commissioning number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

time = strmid(midi,16,24)
nobs = n_elements(midi)

;select reduction method
redmethod = ''
print, ''
print, '     1 EWS'
print, '     2 CohInt'
print, '     3 Koresko'
read, 'For which reduction method?: ', redmethod
if (redmethod eq '1') then $
  resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_EWS_CustomMask/'
if (redmethod eq '2') then $
  resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_CohInt_CustomMask/'
if (redmethod eq '3') then $
  resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_Koresko_CustomMask/'
print, ''



;midipath = '/media/disk/MIDI_FSUA/COMM'+obsnum+'/MIDIdata/'
pcrdir = resultsdir+'Calibrator_PhotonCountRate/'
; workingdir = '/opt/MIDI_FSUA/temp_WorkingDir/'
; file_mkdir, workingdir

;**********************************************************************************************

;extract to which night the files belong

year = double(strmid(midi, 5, 4))
month = double(strmid(midi, 10, 2))
day = double(strmid(midi, 13, 2))
hour = double(strmid(midi, 16, 2))
minute = double(strmid(midi, 19, 2))
second = double(strmid(midi, 22, 2))
jd = julday(month, day, year, hour, minute, second)
jdtmp = jd-0.5d0

night = strarr(nobs, 3)
for i=0,nobs-1 do begin
  if ((jdtmp[i]-floor(jdtmp[i])) gt 0.5d0) then begin
    caldat, jdtmp[i]-1.d0, month1, day1, year1, hour1, minute1, second1
    night[i,*] = strcompress([day1, month1, year1], /rem)
    if (uint(night[i,1] lt 10)) then night[i,1] = '0'+strcompress(night[i,1],/rem)
    if (uint(night[i,0] lt 10)) then night[i,0] = '0'+strcompress(night[i,0],/rem)
  endif else begin
    caldat, jdtmp[i], month1, day1, year1, hour1, minute1, second1
    night[i,*] = strcompress([day1, month1, year1], /rem)
    if (uint(night[i,1] lt 10)) then night[i,1] = '0'+strcompress(night[i,1],/rem)
    if (uint(night[i,0] lt 10)) then night[i,0] = '0'+strcompress(night[i,0],/rem)
  endelse
endfor



;read in the averaged TF for the corresponding night and multiply it with the raw correlated flux
for i=0,nobs-1 do begin
;for i=3,3 do begin

  ;wavelength
  st1 = mrdfits(resultsdir+id[i]+'_'+time[i]+'.corr.fits', 1, /silent)
  wl = 1.d6*st1.eff_wave	;micron

  ;raw correlatd flux
  st3 = mrdfits(resultsdir+id[i]+'_'+time[i]+'.corr.fits', 3, /silent)
  rawFcorr = st3.visamp
  rawFcorr_err = st3.visamperr

  ;mean TF of the night
  readcol, pcrdir+night[i,2]+night[i,1]+night[i,0]+'_'+base[basequest]+'_meanTF_Fcorr_Calibrator.txt', wltf_orig, tf_orig, tferr_orig, format='d,d,d', /silent	;wltf is in micron


  ;sort data with respect to ascending wavelength because for longer wavelengths order is not correct
  idxsort = sort(wl)
  wl = wl[idxsort]
  rawFcorr = rawFcorr[idxsort]
  rawFcorr_err = rawFcorr_err[idxsort]

  idxsorttf = sort(wltf_orig)
  wltf_orig = wltf_orig[idxsort]
  tf_orig = tf_orig[idxsort]
  tferr = tferr_orig[idxsort]

  ;to be sure interpolate TF on wavelength from file
  tf = interpol(tf_orig, wltf_orig, wl, /spline)

  ;compute calibrated correlated flux for each target
  calFcorr = rawFcorr/tf
  calFcorr_err = calFcorr*sqrt((rawFcorr_err/rawFcorr)^2. + (tferr/tf)^2.)


  ;output, wl, fcorr, err
  openw, lun, resultsdir+id[i]+'_'+time[i]+'.calcorr', width=1400, /get_lun

    printf, lun, '## wavelength in micrometer, Fcorr and its error in Jy'
    for j=0,n_elements(calFcorr)-1 do $
      printf, lun, wl[j], calFcorr[j], calFcorr_err[j], format='(f10.5, 2f20.10)'
  close, lun
  free_lun, lun


endfor


;stop
; spawn, 'rm -r '+workingdir

end
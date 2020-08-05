@readcol.pro
@numlines.pro
@remchar.pro
@gettok.pro
@repchr.pro
@strsplit.pro
@strnumber.pro
@valid_num.pro
@fxpar.pro
@uniq.pro
@fxbtform.pro
@fxbfind.pro
@ieee_to_host.pro
@interpol.pro


;computes a calibrated Visibility of Calibrator Stars using the total photometry of the vanBoekel database

;HAS TO BE EXECUTED FROM MIA+EWS!

pro MIDI_Calibrator_calibrate_Vis_CustomMask

;***************************************************************************************

;ewsversion = 'MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013Feb02'
; ewsversion = 'MIA+EWS-2014Feb16'
; dbpath = '/home/amueller/src/'+ewsversion+'/idl/calibrators/databases/'
; db = dbpath+'vBoekelDatabase.fits'

print, ''
print, 'HAS TO BE EXECUTED FROM MIA+EWS!'
print, ''

commnum = ''
read, 'Enter number of Commissioning: ', commnum

readcol, 'observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', skipline=1, /silent

idx = where(comm eq commnum)
if (idx[0] ne -1) then begin

  comm = comm[idx]
  id = id[idx]
  scical = scical[idx]
  nfl = nfl[idx]
  midi = midi[idx]
  fsu = fsu[idx]
  maskname = maskname[idx]
;   photA = photA[idx]
;   photB = photB[idx]

endif else begin

  print, ''
  print, 'Commissioning number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

time = strmid(midi,16,24)

idxsc = where(scical eq 'Cal')
 if (idxsc[0] ne -1) then begin
  comm = comm[idxsc]
  id = id[idxsc]
  scical = scical[idxsc]
  nfl = nfl[idxsc]
  midi = midi[idxsc]
  fsu = fsu[idxsc]
  maskname = maskname[idxsc]
  time = time[idxsc]
  n = n_elements(idxsc)
endif else begin
  print, 'No calibrator targets observed in that run.'
  stop
endelse

;select reduction method
redmethod = ''
print, ''
print, '     1 EWS'
print, '     2 CohInt'
print, '     3 Koresko'
read, 'For which reduction method?: ', redmethod
if (redmethod eq '1') then $
  resultsdir = '/media/disk/MIDI_FSUA/COMM'+commnum+'/MIDIreduced_EWS_CustomMask/'
if (redmethod eq '2') then $
  resultsdir = '/media/disk/MIDI_FSUA/COMM'+commnum+'/MIDIreduced_CohInt_CustomMask/'
if (redmethod eq '3') then $
  resultsdir = '/media/disk/MIDI_FSUA/COMM'+commnum+'/MIDIreduced_Koresko_CustomMask/'
print, ''


n = n_elements(idxsc)

midipath = '/media/disk/MIDI_FSUA/COMM'+commnum+'/MIDIdata/'
pcrdir = resultsdir+'Calibrator_PhotonCountRate/'
file_mkdir, pcrdir
workingdir = '/opt/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir

;***************************************************************************************

cdb = vboekelbase()	;EWS shortcut to the database
for i=0,n-1 do begin

;   counter, i+1, n, 'Processing observation '


;extract data from database
  nfound = cdb->sourceByName(id[i], sourceData, diamData, photData, specData)
  if (total(nfound) eq 0.) then begin
    print, ''
    print, 'Calibrator is not in vanBoekel Database. Stop.'
    print, ''
    print, stop
  endif


  model_Ftot_orig = transpose((*specData[0]).flux)	;Jy
  model_wl_orig = 1.d6*(transpose((*specData[0]).wavelength))	;micrometer
  ;plot, (*specData[0]).wavelength ,(*specData[0]).flux

;get calibrated correlated flux
  readcol, resultsdir+id[i]+'_'+time[i]+'.calcorr', wl, calFcorr, calFcorr_err, format='d,d,d', /silent

  Ftot = interpol(model_Ftot_orig, model_wl_orig, wl, /spline)

  V = calFcorr/Ftot
  IS THE ERROR COMPUTATION CORRECT?? no error in roys data base?
  Verr = V*calFcorr_err/calFcorr

  openw, lun, resultsdir+id[i]+'_'+time[i]+'.calvis', width=1400, /get_lun
    wlout = wl*1.d6
    printf, lun, 'calibrated Visibilities'
    for j=0,n_elements(wl)-1 do begin
      printf, lun, wl[j], V[j], Verr[j], format="(f10.5,2f20.10)"
    endfor

  close, lun
  free_lun, lun

;   window, 0
;   ploterror, wl, v, verr, xr=[8,13], xst=1, yr=[0.6, 1.6], yst=1
;   plots, !x.crange, [1,1]
;   wait, 1.5

endfor


stop
end
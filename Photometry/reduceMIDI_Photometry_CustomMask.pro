@remchar.pro
@counter.pro
@fifteenb.pro
@readcol.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@get_eso_keyword.pro

;reduces MIDI photometry data if present

pro reduceMIDI_Photometry_CustomMask

;**********************************************************************************************

;ewspath = '/opt/MIA+EWS-2011Feb02/c/bin/'
;ewspath = '/opt/MIA+EWS-2009Dec02/c/bin/'

; ewsversion = 'MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013Feb02'
; ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'


obsnum = ''
read, 'Enter Observation ID: ', obsnum


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
  photo = photo[idx]
  photA = photA[idx]
  photB = photB[idx]

endif else begin

  print, ''
  print, 'Commissioning number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

;path = '/opt/MIDI_FSUA/MIDI_FSUA_COMM'+obsnum+'/'	;masks are there
photpath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/Photometry/'
maskpath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDI_CustomMasks/'
; reddir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_Koresko_CustomMask/'
reddir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_faintEWS_CustomMask/'
resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/Photometry_CustomMask/'
file_mkdir, resultsdir
workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir


;**********************************************************************************************

; idxphot = where(photo eq 'y' and nfl ge 10.)	;if ATs used
; if (idxphot[0] ne -1) then begin
; 
;   comm = comm[idxphot]
;   id = id[idxphot]
;   nfl = nfl[idxphot]
;   midi = midi[idxphot]
;   fsu = fsu[idxphot]
;   maskname = maskname[idxphot]
;   photo = photo[idxphot]
; 
; endif

idxphot = where(photo eq 'y')	;if ATs used
if (idxphot[0] ne -1) then begin

  comm = comm[idxphot]
  id = id[idxphot]
  nfl = nfl[idxphot]
  midi = midi[idxphot]
  fsu = fsu[idxphot]
  maskname = maskname[idxphot]
  photo = photo[idxphot]
  photA = photA[idxphot]
  photB = photB[idxphot]

endif

prime = "'"

time = strmid(midi,16,24)

n = n_elements(id)

for i=0,n-1 do begin
;for i=3,3 do begin

  counter, i+1, n, 'Processing observation '
;   print, ''
;   print, 'Photometry '+strcompress(i+1, /rem)+' / '+strcompress(n, /rem)
;   print, ''


  ;Mask selection
  srcmask = maskpath+maskname[i]+'.srcmask.fits'
  skymask = maskpath+maskname[i]+'.skymask.fits'
  ;skymask = srcmask
;   srcmask = '/home/amueller/src/MIA+EWS-2014Feb16/maskfiles/minrtsMask_FIELD_PRISM_HIGH_SENS_SLIT_2011.fits'


;search for photometry files
  date = strmid(midi[i], 0, 16)
  photAfile = file_search(photpath+date+photA[i]+'*.fits', count=nfilesA)
  photBfile = file_search(photpath+date+photb[i]+'*.fits', count=nfilesB)

if (obsnum eq 'HDP91PANIC') then begin
  if (i eq 3) then begin
    tmp = 'MIDI.2013-05-11T'
  ;   photAfile = file_search(photpath+tmp+photA[i]+'*.fits', count=nfilesA)
  photBfile = file_search(photpath+tmp+photb[i]+'*.fits', count=nfilesB)
  endif
  if (i eq 7) then begin
    tmp = 'MIDI.2013-05-12T'
    photAfile = file_search(photpath+tmp+photA[i]+'*.fits', count=nfilesA)
    photBfile = file_search(photpath+tmp+photb[i]+'*.fits', count=nfilesB)
  endif
endif

  ;reduce data and extract spectrum

    if (n_elements(photAfile) eq 1) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photAfile[0]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Aphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photAfile) eq 2) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photAfile[0]+' '+photAfile[1]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Aphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photAfile) eq 3) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photAfile[0]+' '+photAfile[1]+' '+photAfile[2]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Aphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photAfile) eq 4) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photAfile[0]+' '+photAfile[1]+' '+photAfile[2]+' '+photAfile[3]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Aphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photAfile) eq 5) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photAfile[0]+' '+photAfile[1]+' '+photAfile[2]+' '+photAfile[3]+' '+photAfile[4]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Aphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photAfile) eq 6) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photAfile[0]+' '+photAfile[1]+' '+photAfile[2]+' '+photAfile[3]+' '+photAfile[4]+' '+photAfile[5]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Aphotometry.fits -ref '+srcmask;+' -skymask '+skymask

    if (n_elements(photBfile) eq 1) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photBfile[0]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Bphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photBfile) eq 2) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photBfile[0]+' '+photBfile[1]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Bphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photBfile) eq 3) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photBfile[0]+' '+photBfile[1]+' '+photBfile[2]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Bphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photBfile) eq 4) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photBfile[0]+' '+photBfile[1]+' '+photBfile[2]+' '+photBfile[3]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Bphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photBfile) eq 5) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photBfile[0]+' '+photBfile[1]+' '+photBfile[2]+' '+photBfile[3]+' '+photBfile[4]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Bphotometry.fits -ref '+srcmask;+' -skymask '+skymask
    if (n_elements(photBfile) eq 6) then spawn, ewspath+'oirChopPhotoImages -in '+prime+photBfile[0]+' '+photBfile[1]+' '+photBfile[2]+' '+photBfile[3]+' '+photBfile[4]+' '+photBfile[5]+prime+' -out '+workingdir+id[i]+'_'+time[i]+'.Bphotometry.fits -ref '+srcmask;+' -skymask '+skymask

    spawn, ewspath+'oirMakePhotoSpectra -A '+workingdir+id[i]+'_'+time[i]+'.Aphotometry.fits'+' -B '+workingdir+id[i]+'_'+time[i]+'.Bphotometry.fits'+' -mask '+srcmask+' -out '+workingdir+id[i]+'_'+time[i]+'.photometry.fits'; -autoshift'


    spawn, 'cp '+reddir+id[i]+'_'+time[i]+'.corr.fits '+workingdir	;temporary copy, needed by oirRedCal

    spawn, ewspath+'oirRedCal '+workingdir+id[i]+'_'+time[i]

    ;DON'T CHANGE ORDER
    spawn, 'rm '+workingdir+id[i]+'_'+time[i]+'.corr.fits'	;delete temporary copy
    spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'.redcal.fits '+reddir+'/'


    spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'*.fits '+resultsdir+'/'


endfor

spawn, 'rm -r '+workingdir

stop
end
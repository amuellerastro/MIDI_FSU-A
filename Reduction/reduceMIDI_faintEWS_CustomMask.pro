@counter.pro
@fifteenb.pro
@get_eso_keyword.pro

;reduces MIDI data taken with FSU-A as fringe sensor

pro reduceMIDI_faintEWS_CustomMask

print, ''
print, 'Execute this script from MIA+EWS!'
print, ''



;**********************************************************************************************
;ewsversion = 'MIA+EWS-2009Dec02'
;ewsversion = 'MIA+EWS-2011Sep08'
;ewsversion = 'MIA+EWS-2011Oct22'
;ewsversion = 'MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013Feb02'
; ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'


;**********************************************************************************************

obsnum = ''
read, 'Enter Observation Number: ', obsnum


sm = ''
read, 'Smooth (Default: 1.0) : ', sm
gsm = ''
read, 'Gsmooth (Default: 0.1): ', gsm
msm = ''
read, 'Msmooth (Default: 1.0): ', msm
; sm = '0.2'
; gsm = '0.6'
; msm = strcompress(1, /rem)

; mxopd = ''
; read, 'Max OPD (Default: 100): ', mxopd
mxopd = strcompress(100, /rem)

; sm = strcompress(round(double(sm/1.8d-2)), /rem)	;based on midiVisPipe
; strsm = 'smooth='+sm
; strgsm = 'smooth='+gsm
; strmsm = 'msmooth='+msm
; strmxopd = ' -maxopd '+mxopd



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


;path = '/opt/MIDI_FSUA/MIDI_FSUA_COMM'+obsnum+'/'
maskpath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDI_CustomMasks/'
spawn, 'mkdir '+maskpath

midipath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/'
fsupath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'
maskpath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDI_CustomMasks/'
resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_faintEWS_CustomMask'
file_mkdir, resultsdir
; workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
; file_mkdir, workingdir

;**********************************************************************************************

prime = "'"

time = strmid(midi,16,24)

n = n_elements(id)
for i=0,n-1 do begin
; for i=1,1 do begin

;   counter, i+1, n, 'Processing observation '
  print, ''
  print, 'Observation '+strcompress(i+1,/rem)+' / '+strcompress(n,/rem)
  print, ''

  file = file_search(midipath+midi[i]+'*',count=nfiles)

  check99 = strmid(file[nfiles-1],strlen(file[nfiles-1])-7,2)
  if (check99 eq '99') then begin
    file = file[0:nfiles-2]
    nfiles = n_elements(file)
  endif


  ;get DIT
  dum = readfits(file[0], hdr, /silent)
  dit = get_eso_keyword(hdr, 'HIERARCH ESO DET DIT')
;**********************************************************************

  ;DATA BASED ON DEFAULT VALUES OF MIDIUTILITIES OF EWS FAINTVISPIPE.PRO
; 
;   ;default values
;   strSm = ' -smooth '+string(round(1./dit))+' '
;   gsmooth = 0.1
;   strGsm = ' -smooth '+string(gsmooth)+' '
;   psmooth = gsmooth
;   ngrad = 2
;   strNg = ' -ngrad '+string(ngrad)+' '
;   maxopd = 100.
;   strMxopd = '-maxopd '+strtrim(maxopd,2)+' ' 
;   strMnopd = " "
;   strGsf = ' -smooth '+string(psmooth)+' '
;   strMsm = ' -medsmooth '+string(gsmooth)
; ;   noave = 0	;not if data were not recorded with offset tracking

;**********************************************************************

;Mask selection
  mask = maskpath+maskname[i]+'.srcmask.fits'
;   if (scical[i] eq 'Cal') then mask = maskpath+id[i]+'_'+time[i]+'.srcmask.fits'
;   if (scical[i] eq 'Sci') then begin
;     filemask = file_search(maskpath+'*.srcmask.fits', count=nmasks)
;     for j=0,nmasks-1 do print, j+1, '  ', filemask[j]
;     print, ''
;     read, 'Choose a mask from a calibrator for '+id[i]+'_'+time[i]+': ', questmask
;     print, ''
;     mask = filemask[questmask-1]
;   endif

;   if (nfiles eq 1) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 2) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 3) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 4) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 5) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 6) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 7) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 8) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 9) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 10) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 11) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 12) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
;   if (nfiles eq 13) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+' '+file[12]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
; 
;   spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits'+strSm
; ;   if (nfl[i] le 7.) then begin
; ;     spawn, ewspath+'oirFormFringes '+id[i]+'_'+time[i]+'.compressed.fits '+id[i]+'_'+time[i]+'.fringes.fits -removeAverage'
; ;   endif else begin
; ;     spawn, ewspath+'oirFormFringes '+id[i]+'_'+time[i]+'.compressed.fits '+id[i]+'_'+time[i]+'.fringes.fits'
; ;   endelse
; 
;   spawn, ewspath+'oirRotateInsOpd '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.insopd.fits'
; 
;   ;necessary in order to have GD values for each frame
;   spawn, ewspath+'oirGroupDelay '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.ORIGINALgroupdelay.fits'
; 
;   ;twopass groupdelay
;   spawn, ewspath+'oirFGroupDelay '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay1.fits'+ strGsm + strMxopd + strMsm + strNg + ' -exinterval 0.6'
;   spawn, ewspath+'oirPowerDelay '+workingdir+id[i]+'_'+time[i]+'.groupdelay1.fits '+workingdir+id[i]+'_'+time[i]+'.powerdelay.fits'+ strGsm +' -ampsmooth 8.'
;   spawn, ewspath+'oirFGroupDelay '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits -first '+workingdir+id[i]+'_'+time[i]+'.powerdelay.fits' + strNg + strGsm + ' -maxopd 40 '+ strMsm
; 
; 
;   spawn, ewspath+'oirRotateGroupDelay '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits'+strGsf+'  -ngrad 0 -debias 1.5 '
;   spawn, ewspath+'oirAutoFlag '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.flag.fits -maxopd 100 -jitteropd 1.5 ' + strMnopd
;   spawn, ewspath+'oirAverageVis '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.flag.fits '+workingdir+id[i]+'_'+time[i]+'.corr.fits'

  tmp = ''
  for j=0,nfiles-1 do tmp = tmp+' '+file[j]
  file = tmp

  midivispipe, id[i]+'_'+time[i], file, mask=mask, /twopass, smooth=sm, gsmooth=gsm, msmooth=msm, /noVLTIdelay


  spawn, 'mv '+id[i]+'_'+time[i]+'* '+resultsdir+'/ '
;stop
endfor

; spawn, 'rm -r '+workingdir

stop
end

@fxpar.pro
@fxbtform.pro
@fxbfind.pro
@fxbtdim.pro
@ieee_to_host.pro
@fxparpos.pro
@detabify.pro
@interpol.pro
@readfits.pro
@sxpar.pro
@counter.pro
@fifteenb.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@get_eso_keyword.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro

;reduces MIDI data taken with FSU-A as fringe sensor

pro reduceMIDI_EWS_CustomMask


;**********************************************************************************************


;ewsversion = 'MIA+EWS-2009Dec02'
;ewsversion = 'MIA+EWS-2011Sep08'
;ewsversion = 'MIA+EWS-2011Dec13'
;ewsversion = 'MIA+EWS-2013Feb02'
;ewsversion = 'MIA+EWS-2013May18'
ewsversion = 'MIA+EWS-2014Feb16'

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'

obsnum = ''
read, 'Enter Observation Number: ', obsnum
; obsnum = 'P92CHOQUET'

;Choquet: sm=0.2, gsm=0.6
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

sm = strcompress(round(double(sm/1.8d-2)), /rem)	;based on midiVisPipe
strsm = ' -smooth '+sm
strgsm = ' -gsmooth '+gsm
strmsm = ' -medsmooth '+msm
strmxopd = ' -maxopd '+mxopd
ngrad = 2
strNg = ' -ngrad '+strcompress(ngrad, /rem)


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

;path = '/opt/MIDI_FSUA/MIDI_FSUA_COMM'+obsnum+'/'	;masks are there
midipath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/'
fsupath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'
maskpath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDI_CustomMasks/'
resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_EWS_CustomMask'
file_mkdir, resultsdir
workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir

;**********************************************************************************************


prime = "'"

time = strmid(midi,16,24)

n = n_elements(id)
for i=0,n-1 do begin
;for i=1,n-1 do begin
; for i=16,16 do begin

;   counter, i+1, n, 'Processing observation '
  print, ''
  print, 'Observation '+strcompress(i+1,/rem)+' / '+strcompress(n,/rem)
  print, ''


;Mask selection
  mask = maskpath+maskname[i]+'.srcmask.fits'
;    mask = '/home/amueller/src/MIA+EWS-2014Feb16/maskfiles/minrtsMask_FIELD_PRISM_HIGH_SENS_SLIT_2011.fits'
;   if (scical[i] eq 'Cal') then mask = maskpath+id[i]+'_'+time[i]+'.srcmask.fits'
;   if (scical[i] eq 'Sci') then begin
;     filemask = file_search(maskpath+'*.srcmask.fits', count=nmasks)
;     for j=0,nmasks-1 do print, j+1, '  ', filemask[j]
;     print, ''
;     read, 'Choose a mask from a calibrator for '+id[i]+'_'+time[i]+': ', questmask
;     print, ''
;     mask = filemask[questmask-1]
;   endif


  file = file_search(midipath+midi[i]+'*',count=nfiles)

  check99 = strmid(file[nfiles-1],strlen(file[nfiles-1])-7,2)
  if (check99 eq '99') then begin
    file = file[0:nfiles-2]
    nfiles = n_elements(file)
  endif


  ;check if FSUA was used for the noVLTIDelay option
  dum = mrdfits(file[0], 0, hdr, /silent)
  mode = get_eso_keyword(hdr, 'HIERARCH ESO TPL ID')
  if (mode eq 'MIDI_starintf_obs_fringe') then strdelay = ' -noVLTIDelay' else strdelay = ' '


  if (nfiles eq 2) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 3) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 4) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 5) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 6) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 7) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 8) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 9) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 10) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 11) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 12) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 13) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+' '+file[12]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 21) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+' '+file[12]+' '+file[13]+' '+file[14]+' '+file[15]+' '+file[16]+' '+file[17]+' '+file[18]+' '+file[19]+' '+file[20]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 27) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+' '+file[12]+' '+file[13]+' '+file[14]+' '+file[15]+' '+file[16]+' '+file[17]+' '+file[18]+' '+file[19]+' '+file[20]+' '+file[21]+' '+file[22]+' '+file[23]+' '+file[24]+' '+file[25]+' '+file[26]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'


  spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits'+strsm+' -averageOrder 0'+strdelay

  spawn, ewspath+'oirRotateInsOpd '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.insopd.fits'

  spawn, ewspath+'oirFGroupDelay '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay1.fits'+strgsm+strmxopd+strmsm+strNg+' -exinterval 0.6 -sigout .125'

  spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes1.fits '+strsm+strdelay

;   if (nfl[i] le 7.) then begin
;     spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits -removeAverage'
;   endif else begin
;     spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits'
;   endelse
;   if (midi[i] eq 'MIDI.2012-12-15T02:11:26' or midi[i] eq 'MIDI.2012-12-15T02:41:34' or midi[i] eq 'MIDI.2012-12-15T03:30:16') then begin
;     spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits -smooth 200 -averageOrder 1'
;   endif else begin
;     spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits'
;   endelse


;   spawn, ewspath+'oirRotateInsOpd.noOPD '+id[i]+'_'+time[i]+'.fringes.fits '+id[i]+'_'+time[i]+'.insopd.fits'
;   spawn, ewspath+'oirGroupDelay.noOPD '+id[i]+'_'+time[i]+'.insopd.fits '+id[i]+'_'+time[i]+'.groupdelay.fits'
;   spawn, ewspath+'oirRotateGroupDelay.noOPD '+id[i]+'_'+time[i]+'.fringes.fits '+id[i]+'_'+time[i]+'.groupdelay.fits '+id[i]+'_'+time[i]+'.ungroupdelay.fits'
;   spawn, ewspath+'oirAutoFlag.noOPD '+id[i]+'_'+time[i]+'.ungroupdelay.fits '+id[i]+'_'+time[i]+'.groupdelay.fits '+id[i]+'_'+time[i]+'.flag.fits'
;   spawn, ewspath+'oirAverageVis '+id[i]+'_'+time[i]+'.ungroupdelay.fits '+id[i]+'_'+time[i]+'.flag.fits '+id[i]+'_'+time[i]+'.corr.fits'


;   spawn, ewspath+'oirGroupDelay '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits'

  spawn, ewspath+'oirRotateGroupDelay '+workingdir+id[i]+'_'+time[i]+'.fringes1.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay1.fits '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits'+strgsm+'  -ngrad 0 -debias 1.5 '+strdelay

  spawn, ewspath+'oirAutoFlag '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay1.fits '+workingdir+id[i]+'_'+time[i]+'.flag.fits -maxopd 100 -jitteropd 1.5'

  spawn, ewspath+'oirAverageVis '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.flag.fits '+workingdir+id[i]+'_'+time[i]+'.corr.fits'

  spawn, 'mv '+workingdir+id[i]+'* '+resultsdir+'/'

endfor

spawn, 'rm -r '+workingdir
;stop
end

@counter.pro
@fifteenb.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@get_eso_keyword.pro
@fxpar.pro
@fxbtform.pro
@fxbfind.pro
@fxbtdim.pro
@ieee_to_host.pro
@fxparpos.pro
@detabify.pro
@interpol.pro

;reduces MIDI data taken with FSU-A as fringe sensor

pro reduceMIDI_EWS_CustomMask_midiVisPipe


;**********************************************************************************************


;ewsversion = 'MIA+EWS-2009Dec02'
;ewsversion = 'MIA+EWS-2011Sep08'
;ewsversion = 'MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013Feb02'
; ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
; ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'

obsnum = ''
read, 'Enter Observation Number: ', obsnum
; obsnum = 'P92CHOQUET'
; obsnum = 'P93MONNIERA1B2'
; obsnum = 'P93MONNIERA1C1'
; obsnum = 'P93MONNIERB2C1'

;Choquet: sm=0.2, gsm=0.6
sm = ''
read, 'Smooth (Default: 1.0) : ', sm
gsm = ''
read, 'Gsmooth (Default: 0.1): ', gsm
msm = ''
read, 'Msmooth (Default: 1.0): ', msm
; sm = '0.2'
; gsm = '0.6'


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
; fsupath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'
maskpath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDI_CustomMasks/'
resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_EWS_CustomMask'
file_mkdir, resultsdir
; workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
; file_mkdir, workingdir

;**********************************************************************************************


prime = "'"

time = strmid(midi,16,24)

n = n_elements(id)
for i=0,n-1 do begin
;for i=5,5 do begin

;   counter, i+1, n, 'Processing observation '
  print, ''
  print, 'Observation '+strcompress(i+1,/rem)+' / '+strcompress(n,/rem)
  print, ''


;Mask selection
  mask = maskpath+maskname[i]+'.srcmask.fits'


  file = file_search(midipath+midi[i]+'*',count=nfiles)

  check99 = strmid(file[nfiles-1],strlen(file[nfiles-1])-7,2)
  if (check99 eq '99') then begin
    file = file[0:nfiles-2]
    nfiles = n_elements(file)
  endif

  ;f = [file[0]+' '+file[1]+' '+file[2]]
  f = ''
  for x=0,n_elements(file)-1 do f = f+' '+file[x]

  ;for P93MONNIERA1B2 with MIDI stand alone track/notrack data
  ;but with or without /noVLTIDelay produces the same result for notrack data
;   if (id[i] eq 'HD187642') then midivispipe, id[i]+'_'+time[i], f, mask=mask, smooth=sm, ave=0, gsmooth=gsm, msmooth=msm
;   if (id[i] eq 'HD187929' or id[i] eq 'V1295Aql') then midivispipe, id[i]+'_'+time[i], f, mask=mask, smooth=sm, ave=0, gsmooth=gsm, msmooth=msm;, /noVLTIDelay

  midivispipe, id[i]+'_'+time[i], f, mask=mask, smooth=sm, ave=0, gsmooth=gsm, msmooth=msm, /noVLTIDelay

  resfile = file_search(id[i]+'_'+time[i]+'*.fits', count=nresfiles)
  for j=0,nresfiles-1 do spawn, 'mv '+resfile[j]+' '+resultsdir+'/'

endfor

; spawn, 'rm -r '+workingdir
;stop
end

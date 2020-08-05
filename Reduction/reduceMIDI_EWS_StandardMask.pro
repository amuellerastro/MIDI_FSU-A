@counter.pro
@fifteenb.pro
@get_eso_keyword.pro

;reduces MIDI data taken with FSU-A as fringe sensor

pro reduceMIDI_EWS_StandardMask


;**********************************************************************************************


;ewsversion = 'MIA+EWS-2009Dec02'
;ewsversion = 'MIA+EWS-2011Sep08'
; ewsversion = 'MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013Feb02'
; ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'
mask = '/home/amueller/src/'+ewsversion+'/maskfiles/minrtsMask_FIELD_PRISM_HIGH_SENS_SLIT_2011.fits'

commnum = ''
read, 'Enter Observation ID: ', commnum


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
  photA = photA[idx]
  photB = photB[idx]

endif else begin

  print, ''
  print, 'Commissioning number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

midipath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/MIDIdata/'
fsupath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/FSUAdata/'
resultsdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/MIDIreduced_EWS_StandardMask'
file_mkdir, resultsdir
workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir


;**********************************************************************************************


prime = "'"

time = strmid(midi,16,24)

n = n_elements(id)
for i=0,n-1 do begin
;for i=42,42 do begin

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

  spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits'
;   if (nfl[i] le 7.) then begin
;     spawn, ewspath+'oirFormFringes '+id[i]+'_'+time[i]+'.compressed.fits '+id[i]+'_'+time[i]+'.fringes.fits -removeAverage'
;   endif else begin
;     spawn, ewspath+'oirFormFringes '+id[i]+'_'+time[i]+'.compressed.fits '+id[i]+'_'+time[i]+'.fringes.fits'
;   endelse

;   spawn, ewspath+'oirRotateInsOpd.noOPD '+id[i]+'_'+time[i]+'.fringes.fits '+id[i]+'_'+time[i]+'.insopd.fits'
;   spawn, ewspath+'oirGroupDelay.noOPD '+id[i]+'_'+time[i]+'.insopd.fits '+id[i]+'_'+time[i]+'.groupdelay.fits'
;   spawn, ewspath+'oirRotateGroupDelay.noOPD '+id[i]+'_'+time[i]+'.fringes.fits '+id[i]+'_'+time[i]+'.groupdelay.fits '+id[i]+'_'+time[i]+'.ungroupdelay.fits'
;   spawn, ewspath+'oirAutoFlag.noOPD '+id[i]+'_'+time[i]+'.ungroupdelay.fits '+id[i]+'_'+time[i]+'.groupdelay.fits '+id[i]+'_'+time[i]+'.flag.fits'
;   spawn, ewspath+'oirAverageVis '+id[i]+'_'+time[i]+'.ungroupdelay.fits '+id[i]+'_'+time[i]+'.flag.fits '+id[i]+'_'+time[i]+'.corr.fits'

  spawn, ewspath+'oirRotateInsOpd '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.insopd.fits'
  spawn, ewspath+'oirGroupDelay '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits'
  spawn, ewspath+'oirRotateGroupDelay '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits'
  spawn, ewspath+'oirAutoFlag '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.flag.fits'
  spawn, ewspath+'oirAverageVis '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.flag.fits '+workingdir+id[i]+'_'+time[i]+'.corr.fits'


  spawn, 'mv '+workingdir+id[i]+'* '+resultsdir
;stop
endfor

spawn, 'rm -r '+workingdir
;stop
end

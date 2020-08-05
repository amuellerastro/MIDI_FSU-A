;reduces MIDI data taken with FSU-A as fringe sensor
;calculates correlated flux using flaging by analyze_FSUA_GD_MIDICohInt.pro
;is not using GD from MIDI or FSU-A at all

@counter.pro
@fifteenb.pro
@get_eso_keyword.pro
@sync_FSUAMIDI.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@readfits.pro
@sxpar.pro
@mean.pro
@moment.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@mrd_hread.pro
@ts_diff.pro
@numlines.pro
@repchr.pro


pro reduceMIDI_CohInt_CustomMask

;**********************************************************************************************

;ewspath = '/opt/MIA+EWS-2009Dec02/c/bin/'
;ewspath = '/opt/MIA+EWS-2011Sep08/c/bin/'
;ewspath = '/opt/MIA+EWS-2011Oct22/c/bin/'

; ewsversion = 'MIA+EWS_OLD/MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'

obsnum = ''
read, 'Enter Observing Number: ', obsnum


sm = ''
read, 'Smooth (Default: 1.0) : ', sm
sm = strcompress(round(double(sm/1.8d-2)), /rem)	;based on midiVisPipe
strsm = ' -smooth '+sm

readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', /silent;, skipline=1

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
  print, 'Observing number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

;path = '/opt/MIDI_FSUA/MIDI_FSUA_COMM'+obsnum+'/'	;masks are there
midipath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/'
fsupath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'
maskpath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDI_CustomMasks/'
resultsdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_CohInt_CustomMask/'
midiewspath  = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_EWS_CustomMask/'
file_mkdir, resultsdir
workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir

;**********************************************************************************************

prime = "'"
time = strmid(midi,16,24)

n = n_elements(id)
for i=0,n-1 do begin
; for i=11,n-1 do begin

;   counter, i+1, n, 'Processing observation '
  print, ''
  print, 'Observation '+strcompress(i+1,/rem)+' / '+strcompress(n,/rem)
  print, ''


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

  if (nfiles eq 1) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
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
;  spawn, ewspath+'oir1dCompressData '+prime+midipath+midi[i]+'.000_01.fits '+midipath+midi[i]+'.000_02.fits '+midipath+midi[i]+'.000_03.fits '+midipath+midi[i]+'.000_04.fits '+midipath+midi[i]+'.000_05.fits '+midipath+midi[i]+'.000_06.fits '+midipath+midi[i]+'.000_07.fits '+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'_compressed.fits'

;   par_smooth = '10'	;smooth value/dit, value = 10*0.018
  spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits'+strsm+' -averageOrder 0'+strdelay; -removeAverage'; -smooth '+par_smooth; -removeAverage'
;   par_smooth = '10'	;smooth value/dit, value = 10*0.018
;   if (nfl[i] le 7.) then begin
;     spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits  -smooth '+par_smooth+' -removeAverage'
;   endif else begin
;     spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits -smooth '+par_smooth
;   endelse

;    spawn, ewspath+'oirRotateInsOpd.noOPD '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.insopd.fits'
;  spawn, ewspath+'oirRotateInsOpdmod '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.insopd.fits'
  spawn, ewspath+'oirRotateInsOpd '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.insopd.fits'


  sync_FSUAMIDI, id[i]+'_'+time[i],workingdir+id[i]+'_'+time[i]+'.compressed.fits', fsupath+fsu[i], obsnum, workingdir

  spawn, ewspath+'oirAverageVis.cohint '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.fsuflags '+workingdir+id[i]+'_'+time[i]+'.corr.fits'

  spawn, 'mv '+workingdir+id[i]+'* '+resultsdir+'/'

endfor

spawn, 'rm -r '+workingdir

;stop
end
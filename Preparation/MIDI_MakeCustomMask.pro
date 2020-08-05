@counter.pro
@fifteenb.pro
@fxpar.pro
@fxbtform.pro
@fxbfind.pro
@fxbtdim.pro
@ieee_to_host.pro
@fxparpos.pro
@detabify.pro
@interpol.pro
@curvefit.pro
@poly_fit.pro
@poly.pro
@fxbhmake.pro
@fxaddpar.pro
@fxhmake.pro
@get_date.pro
@daycnv.pro
@host_to_ieee.pro
@counter.pro
@fifteenb.pro
@fxpar.pro
@fxbtform.pro
@fxbfind.pro
@fxbtdim.pro
@ieee_to_host.pro
@fxparpos.pro
@detabify.pro
; @ls2fit.pro
@interpol.pro
@curvefit.pro
@poly_fit.pro
@poly.pro
@fxbhmake.pro
@fxaddpar.pro
@fxhmake.pro
@get_date.pro
@daycnv.pro
@host_to_ieee.pro
@readcol.pro
@numlines.pro
@remchar.pro
@gettok.pro
@repchr.pro
@strsplit.pro
@strnumber.pro
@valid_num.pro

;USE THIS SCRIPT WITH MIA+EWS

pro MIDI_MakeCustomMask

;**********************************************************************************************

;ewsversion = 'MIA+EWS-2011Dec13'
;ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'


print, ''
print, 'USE THIS SCRIPT WITH MIA+EWS!'
print, ''

obsnum = ''
read, 'Enter Observing Run: ', obsnum


readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', /silent;, skipline=1

;select observing run
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


;select only Calibrators
idxcal = where(scical eq 'Cal')
if (idxcal[0] ne -1) then begin
  comm = comm[idxcal]
  id = id[idxcal]
  scical = scical[idxcal]
  nfl = nfl[idxcal]
  midi = midi[idxcal]
  fsu = fsu[idxcal]
  maskname = maskname[idxcal]
  photA = photA[idxcal]
  photB = photB[idxcal]

endif

prime = "'"
time = strmid(midi,16,24)

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'

midipath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/'

workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir

path = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/
maskpath = path+'MIDI_CustomMasks/'
file_mkdir, maskpath


;**********************************************************************************************

n = n_elements(id)
for i=0,n-1 do begin
;for i=1,1 do begin

  counter, i+1, n, 'Processing observation '
;   print, ''
;   print, 'Observation '+strcompress(i+1,/rem)+' / '+strcompress(n,/rem)
;   print, ''

  file = file_search(midipath+midi[i]+'*',count=nfiles)

  check99 = strmid(file[nfiles-1],strlen(file[nfiles-1])-7,2)
  if (check99 eq '99') then begin
    file = file[0:nfiles-2]
    nfiles = n_elements(file)
  endif


  if (nfiles eq 1) then filetmp = file[0]
  if (nfiles eq 2) then filetmp = file[0]+' '+file[1]
  if (nfiles eq 3) then filetmp = file[0]+' '+file[1]+' '+file[2]
  ;we have to limit the number of used files here because of limit memory (RAM)
  if (nfiles gt 3) then filetmp = file[0]+' '+file[1]+' '+file[2]
;   if (nfiles eq 4) then filetmp = file[0]+' '+file[1]+' '+file[2]+' '+file[3]
;   if (nfiles eq 5) then filetmp = file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]
;   if (nfiles eq 6) then filetmp = file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]
;   if (nfiles eq 7) then filetmp = file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]

;   midiMakeMask, workingdir+id[i]+'_'+time[i], filetmp, smooth=10., gsmooth=5., /noVLTIDelay ;,factor=1.5,  smooth=double(par_smooth), gsmooth=double(par_gsmooth)
;   midiMakeMask, workingdir+id[i]+'_'+time[i], filetmp, smooth=0.2, gsmooth=0.6, /noVLTIDelay, factor=1.2
; if (id[i] eq 'HD187929') then begin

;   midiMakeMask, workingdir+id[i]+'_'+time[i], filetmp, smooth=0.2, gsmooth=0.6, /noVLTIDelay, factor=1.0
  midiMakeMask, workingdir+id[i]+'_'+time[i], filetmp, smooth=0.2, gsmooth=0.2;, /noVLTIDelay;, factor=1.0
;   midiMakeMask, workingdir+id[i]+'_'+time[i], filetmp, smooth=0.2, gsmooth=0.2, factor=1.0

  spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'*.fits '+maskpath

endfor


spawn, 'rm -r '+workingdir

;stop
end
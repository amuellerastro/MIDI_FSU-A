@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@readfits.pro
@sxpar.pro
@get_eso_keyword.pro

pro MIDI_number_of_good_frames_StandardMask

commnum = ''
read, 'Enter number of Commissioning: ', commnum

;*************************************************************
main = '/media/disk/MIDI_FSUA/COMM'+commnum+'/'
resultpath = '/opt/MIDI_FSUA/MIDI_FSUA_COMM'+commnum+'/'
pathews = 'MIDIreduced_EWS_StandardMask/'
pathcohint = 'MIDIreduced_CohInt_StandardMask/'
pathkoresko = 'MIDIreduced_Koresko_StandardMask/'

; pathfaintews = 'MIDIreduced_faintEWS/'


;*************************************************************

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


path = [main+pathews, main+pathcohint, main+pathkoresko]	;main+pathfaintews

ngood = intarr(n_elements(path), n_elements(midi)) & nframes = intarr(n_elements(midi))

for i=0,n_elements(path)-1 do begin

  for j=0,n_elements(midi)-1 do begin

    time = strmid(midi[j], 16, 8)
    file = path[i]+id[j]+'_'+time+'.corr.fits'

    ;read in fsu flag file to get total number of recorded MIDI frames
    flagfile = main+pathcohint+id[j]+'_'+time+'.fsuflags'
    readcol, flagfile, ntmp, format = 'i', /silent

    openr, lun, flagfile, /get_lun
    rows = file_lines(flagfile)
    data = dblarr(1, rows)  
    readf, lun, data
      tmp = data[0,*]
      ntmp = reform(tmp)
    close, lun
    free_lun, lun
    if (i eq 0) then nframes[j] = n_elements(ntmp)

    dum = readfits(file, hdr, /silent)
    ngood[i,j] = get_eso_keyword(hdr, 'NGOOD')

  endfor

endfor

ngoodews = ngood[0,*]
ngoodcohint = ngood[1,*]
ngoodkoresko = ngood[2,*]

openw, lun, resultpath+'MIDI_good_frames_StandardMask.txt', width=1400, /get_lun

  printf, lun, '                  ID                          File Total Frames    EWS    ChoInt    Koresko'

  for i=0,n_elements(midi)-1 do begin

    printf, lun, id[i], midi[i], nframes[i,0], ngoodews[i], ngoodcohint[i], ngoodkoresko[i], format='(a20, a30, i10, i10, i10, i10)'

  endfor



close, lun
free_lun, lun


stop
end
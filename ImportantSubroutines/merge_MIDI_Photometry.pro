;merge two photometric observations

pro merge_MIDI_Photometry

;ASSUMES 2 seperate files only, i.e. _01.fits and _02.fits

indmidi = ['MIDI.2013-03-03T08:36:53', 'MIDI.2013-03-03T08:44:20']
finalmidi = indmidi[0]
;phot B OPEN
indmidi1 = ['MIDI.2013-03-03T08:39:47', 'MIDI.2013-03-03T08:47:16']
finalmidi1 = indmidi1[0]

; ;phot A OPEN
; indmidi = ['MIDI.2014-03-17T08:31:19', 'MIDI.2014-03-17T08:41:18']
; finalmidi = 'MIDI.2014-03-17T08:31:19'
; ;phot B OPEN
; indmidi1 = ['MIDI.2014-03-17T08:34:13', 'MIDI.2014-03-17T08:44:12']
; finalmidi1 = 'MIDI.2014-03-17T08:34:13'

; ;phot A OPEN
; indmidi = ['MIDI.2013-05-28T07:59:02', 'MIDI.2013-05-28T08:12:27']
; finalmidi = 'MIDI.2013-05-28T07:59:02'
; ;phot B OPEN
; indmidi1 = ['MIDI.2013-05-28T08:01:18', 'MIDI.2013-05-28T08:14:43']
; finalmidi1 = 'MIDI.2013-05-28T08:01:18'

;phot A OPEN
; indmidi = ['MIDI.2013-05-27T08:32:54', 'MIDI.2013-05-27T08:44:53']
; finalmidi = 'MIDI.2013-05-27T08:32:54'
; ;phot B OPEN
; indmidi1 = ['MIDI.2013-05-27T08:35:10', 'MIDI.2013-05-27T08:47:09']
; finalmidi1 = 'MIDI.2013-05-27T08:35:10'


backupmidi = 'not_merged/'
spawn, 'mkdir '+backupmidi

;==========================================================================

;PHOTOMOTRY A

;backup the individual files before merging
nobs = n_elements(indmidi)
for i=0,nobs-1 do begin

  if (i eq 0) then spawn, 'cp '+indmidi[i]+'* '+backupmidi $
    else spawn, 'mv '+indmidi[i]+'* '+backupmidi

endfor


;extract the last number
reffile = file_search(backupmidi+indmidi[0]+'*', count=ln)

;extract last frame number of reference file
refst9 = mrdfits(reffile[ln-1], 7, /silent)
refframe = refst9[n_elements(refst9.frame)-1].frame


file = file_search(backupmidi+indmidi[1]+'*', count=nfiles)
nframes = uintarr(nfiles)
for i=0,nfiles-1 do begin

  st9 = mrdfits(file[i], 7, /silent)
  nframes[i] = n_elements(st9.frame)

  if (i eq 0) then framenum = lindgen(nframes[i])+long(1)+refframe $
    else framenum = lindgen(nframes[i])+long(1)+refframe+long(total(nframes[0:i-1]))
  ;print, framenum[0], framenum[n_elements(framenum)-1]

  ;laufender index
  num = ln+i+1
  if (num lt 10) then num = '0'+strcompress(num, /rem) $
    else num = strcompress(num, /rem)

  ;write MIDI files with updated framenumbers
  for j=0,7 do begin

    data = mrdfits(file[i], j, hdr, /silent)

    if (j ne 7) then begin 

      mwrfits, data, finalmidi+'.000_'+num+'.fits', hdr, /silent

    endif else begin

      data.frame = framenum
      mwrfits, data, finalmidi+'.000_'+num+'.fits', hdr, /silent

    endelse

  endfor

endfor

;==========================================================================

;PHOTOMOTRY B

;backup the individual files before merging
nobs = n_elements(indmidi1)
for i=0,nobs-1 do begin

  if (i eq 0) then spawn, 'cp '+indmidi1[i]+'* '+backupmidi $
    else spawn, 'mv '+indmidi1[i]+'* '+backupmidi

endfor


;extract the last number
reffile = file_search(backupmidi+indmidi1[0]+'*', count=ln)

;extract last frame number of reference file
refst9 = mrdfits(reffile[ln-1], 7, /silent)
refframe = refst9[n_elements(refst9.frame)-1].frame


file = file_search(backupmidi+indmidi1[1]+'*', count=nfiles)
nframes = uintarr(nfiles)
for i=0,nfiles-1 do begin

  st9 = mrdfits(file[i], 7, /silent)
  nframes[i] = n_elements(st9.frame)

  if (i eq 0) then framenum = lindgen(nframes[i])+long(1)+refframe $
    else framenum = lindgen(nframes[i])+long(1)+refframe+long(total(nframes[0:i-1]))
  ;print, framenum[0], framenum[n_elements(framenum)-1]

  ;laufender index
  num = ln+i+1
  if (num lt 10) then num = '0'+strcompress(num, /rem) $
    else num = strcompress(num, /rem)

  ;write MIDI files with updated framenumbers
  for j=0,7 do begin

    data = mrdfits(file[i], j, hdr, /silent)

    if (j ne 7) then begin 

      mwrfits, data, finalmidi1+'.000_'+num+'.fits', hdr, /silent

    endif else begin

      data.frame = framenum
      mwrfits, data, finalmidi1+'.000_'+num+'.fits', hdr, /silent

    endelse

  endfor

endfor

stop
end
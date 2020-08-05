@counter.pro
@fifteenb.pro
@arrdelete.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@readfits.pro
@sxpar.pro
@get_eso_keyword.pro
@mean.pro
@moment.pro
@mwrfits.pro
@fxaddpar.pro
@fxparpos.pro
@detabify.pro 
@sixlin.pro
@ts_diff.pro
@sxaddhist.pro

pro MIDI_NULL_and_replace_OPD_values_ZERO



;VAARIABLES
;**********************************************************************************************
obsnum = ''
read, 'Enter Observing Run: ', obsnum

midipath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/'
; fsupath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'
copypath = midipath+'originalMIDIFT/'
nullpath = '/home/amueller/work/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/'
file_mkdir, nullpath
file_mkdir, copypath

readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', skipline=1, /silent

;**********************************************************************************************

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

;**********************************************************************************************

for i=0,n_elements(midi)-1 do begin
; for i=18,n_elements(midi)-1 do begin
; for i=0,0 do begin

  counter, i, n_elements(midi), 'Processing observation '

  ;copy original MIDI fringe track files
  print, 'backup original MIDI FT data '+strcompress(i+1,/rem)+' / '+strcompress(n_elements(midi),/rem)
  midi_file_orig = file_search(midipath+midi[i]+'*.fits', count=nfiles)
  for k=0,nfiles-1 do  spawn, 'mv '+midi_file_orig[k]+' '+copypath

  midi_file = file_search(copypath+midi[i]+'*.fits', count=nfiles)

  ;read in header of all files to get total number of subintegrations
  ndits = intarr(nfiles)
  name = strarr(nfiles)
  for k=0,nfiles-1 do begin

    pos = strpos(midi_file[k], '/', /reverse_search)
    name[k] = strmid(midi_file[k], pos+1, 36)

    dum = readfits(midi_file[k], hdr, /silent)
    ndits[k] = get_eso_keyword(hdr, 'HIERARCH ESO DET NDIT')

    if (k eq 0) then begin

      midi_inttime = double(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'))	;MIDI integration time in seconds
      midi_cycle = double(get_eso_keyword(hdr,'HIERARCH ESO INS PIEZ CYCLT'))	;MIDI cycle time in seconds

    endif

  endfor
  totframes = total(ndits)	;with that number we can create arrays of correct size

  ;read in time of all files into a single array
  timeall = dblarr(totframes)
  for k=0,nfiles-1 do begin

    st9 = mrdfits(midi_file[k], 9, /silent)
    idx0 = where(timeall eq 0.d0)
    timeall[idx0[0]:idx0[0]+n_elements(st9.time)-1] = st9.time

  endfor


  ;check if where are NULL time stamps
  print, 'checking for NULL time stamps'
  NULLindex = where(finite(timeall) eq 0)

  if (NULLindex[0] ne -1) then begin

    print, ' '
    print, '***********************************'
    print, 'NULL TIMESTAMP AT FRAME(S) '+strcompress(NULLindex+1,/rem)
    print, '***********************************'
    print, ' '

    ndata = total(ndits)	;total number of frames

    ;replace all NULL time stamps of timeall
    timeallclean = dblarr(ndata)
    idx = where(finite(timeall) ne 0)
    timeallclean[0:n_elements(idx)-1] = timeall[idx]
    ;extrapolate the missing values by comuting average 
    idx1 = where(timeallclean ne 0.d0)
    diff = abs(median(ts_diff(timeallclean[idx1], 1)))	;days
    idx0 = where(timeallclean eq 0.d0)
    for m=0,n_elements(idx0)-1 do timeallclean[idx0[m]] = timeallclean[idx0[m]-1]+diff


    ;read in all corresponding MIDI files
    for k=0,nfiles-1 do begin

      d9 = mrdfits(midi_file[k], 9, hdr, /silent)
      if (k eq 0) then d9all = d9 else d9all = [d9all,d9]

    endfor

 
    ;write frame with NULL time stamp into a file
    for l=0L,n_elements(NULLindex)-1 do begin

      openw, lun, nullpath+'MIDI_NULLvalues.txt', width=1400, /get_lun, /append
	printf, lun, id[i], midi[i], strcompress(NULLindex[l]+1,/rem), format='(a20, a25, a10)'
      close, lun
      free_lun, lun

    endfor

;       endif

  ;replace NULL time stamp with next time stamp
    midi_time = timeallclean;[0:ndits[k]-1]

    ;replace OPD table with zeros
    d9all.opd[1,*] = 0.d0

    ;fix NULL entries
;     for l=0,n_elements(NULLindex)-1 do begin
; 
;       ;insert exptime at NULL
;       d9all[NULLindex[l]].exptime = midi_inttime
; 
;       ;1st column of OPD, insert 0 at NULL
;       d9all[NULLindex[l]].opd[0,*] = 0.d0
; 
;       ;1st column of OPD, insert 0 at NULL
;       d9all[NULLindex[l]].localopd[1] = 0.d0
;       lo = d9all.localopd[0,*]
; 
;       ;interpolate NULL in 1st column of localopd
;       step_orig = d9all.stepping_phase
;       max_step = max(step_orig)
;       null_step = step_orig[NULLindex[l]]
;       null_time = midi_time[NULLindex[l]]
; 
;       t = midi_time[NULLindex[l]-null_step+1:NULLindex[l]+abs(null_step-max_step)]
;       lo = lo[NULLindex[l]-null_step+1:NULLindex[l]+abs(null_step-max_step)]
;       idx = where(finite(lo) ne 0)
;       t = t[idx]
;       lo = lo[idx]
;       sixlin, t, lo, intercept, sigi, slope, sigs
;       lo_new = slope[0]*null_time+intercept[0]
;       d9all[NULLindex[l]].localopd[0] = lo_new
; 
;       ;subtract ~2.1mu from localopd
;       diff = abs(ts_diff(lo, 1))
;       diff = median(diff[0:n_elements(diff)-2])
; 
;       d9all[NULLindex[l]:ndata-1].localopd[0] = d9all[NULLindex[l]:ndata-1].localopd[0]-diff
; 
;     endfor

    mtimes = median((ts_diff(timeallclean,1))[0:n_elements(timeallclean)-2])
    nfram = n_elements(d9all.frame)

    tmp = ts_diff(d9all.localopd[0], 1)
    idx = where(tmp le -4.d-5)
;     mlopds = (median(tmp[idx]))
    mlopds = median(d9all.localopd[0,*] - shift(d9all.localopd[0],1))

    for l=0,n_elements(NULLindex)-1 do begin

      ; we simply anticipate that the NULL entries are save one frame later, since
      ;    a raw data set reduced seem to have a jump in delay by one scanning piezo-step, so
      ;    data and localopd are shifted wrt each other

      d9all[NULLindex[l]:nfram-2].time           = d9all[NULLindex[l]+1:*].time
      d9all[NULLindex[l]:nfram-2].localopd[0]    = d9all[NULLindex[l]+1:*].localopd[0]
      d9all[NULLindex[l]:nfram-2].stepping_phase = d9all[NULLindex[l]+1:*].stepping_phase
      d9all[NULLindex[l]].exptime = midi_inttime
      d9all[NULLindex[l]].localopd[1,*] = 0.d0
      d9all[NULLindex[l]].opd[0,*] = 0.d0

      ; fill up the last entry
      d9all[nfram-1].time = d9all[nfram-2].time + mtimes

      if d9all[nfram-2].stepping_phase eq 40 then begin

        d9all[nfram-1].stepping_phase = 1
        d9all[nfram-1].localopd[0] = d9all[nfram-2].localopd[0]-39.d0*mlopds

      endif else begin

        d9all[nfram-1].stepping_phase = d9all[nfram-2].stepping_phase + 1
        d9all[nfram-1].localopd[0]    = d9all[nfram-2].localopd[0]    + mlopds

      endelse

      NULLindex = NULLindex-1

    endfor 

    d9all.time = midi_time

    for j=0,9 do begin

      for k=0,nfiles-1 do begin

	struc = mrdfits(midi_file[k], j, hdr, status=status, /silent)
	if (j ne 9) then begin

	  mwrfits, struc, midipath+name[k], hdr, /silent

	endif else begin

	  if (k eq 0) then mwrfits, d9all[0:ndits[0]-1], midipath+name[k], hdr, /silent else $;, /create
	  mwrfits, d9all[total(ndits[0:k-1]):total(ndits[0:k])-1], midipath+name[k], hdr, /silent

	endelse

      endfor

    endfor


  endif else begin

;if there is no NULL just write out new files with OPD=0

    for k=0,nfiles-1 do begin	;

      for j=0,9 do begin

	struc = mrdfits(midi_file[k], j, hdr, status=status, /silent)

	if (j ne 9) then begin

	  mwrfits, struc, midipath+name[k], hdr, /silent;, /create

	endif else begin

	  struc.opd[1,*] = 0.d0
	  mwrfits, struc, midipath+name[k], hdr, /silent;, /create

	endelse

      endfor

    endfor

  endelse

endfor	;main loop

stop
end


@get_eso_keyword.pro
; @showsym.pro

function closest,vector,find,SORT=s
  nf=n_elements(find)
  sort=keyword_set(s) || arg_present(s)
  if sort && n_elements(s) ne n_elements(vector) then s=sort(vector)
  j=value_locate(sort?vector[s]:vector,find)
  b=[[j>0],[(j+1)<(n_elements(vector)-1)]]
  mn=min(abs((sort?vector[s[b]]:vector[b])- $
             rebin([find],nf,2)),DIMENSION=2,pos)
  pos=j>0+pos/nf
  return,sort?s[pos]:pos
end 


pro sync_FSUAMIDI, star, midi_file, fsu_file, commnum, workingdir

print, ''
print, 'synchronize MIDI with FSU data...takes some time'

;MIDI
d1 = mrdfits(midi_file, 2, midi_hdr, /silent)
d2 = readfits(midi_file, midi_hdr, /silent)

midi_time = d1.time	;time in table is the end of the single frame!


;midi_time_orig = midi_time
midi_inttime = double(get_eso_keyword(midi_hdr,'HIERARCH ESO DET DIT'))*1.D6	;MIDI integration time in usec


;FSU
; if (commnum eq '6' or commnum eq '8' or commnum eq '9' or commnum eq '14') then begin
;   d3 = mrdfits(fsu_file,13)
; endif else begin
;   d3 = mrdfits(fsu_file,11, /silent)
; endelse
d3 = mrdfits(fsu_file,'OPDC', /silent)

fsu_time_double = d3.time	;usec, double frequency of fsu_time2
fsu_state_orig = d3.state	;OPDC state


d4 = readfits(fsu_file,fsu_hdr, /silent)
; if (commnum eq '6' or commnum eq '8' or commnum eq '9' or commnum eq '14') then begin
;   d5 = mrdfits(fsu_file,9)
; endif else begin
;   if (commnum eq 'P90LOPEZ') then d5 = mrdfits(fsu_file,7, /silent) else d5 = mrdfits(fsu_file,6, /silent)
; endelse
d5 = mrdfits(fsu_file,'IMAGING_DATA_FSUA', /silent)

fsu_gd = d5.gd	;meter, FSU GD
fsu_time = d5.time	;usec, this should be the reference

;starting time MJD
fsu_mjd = double(sxpar(fsu_hdr,'MJD-OBS'))

;extract starting time in UTC
fsu_utc_dum = get_eso_keyword(fsu_hdr,'DATE-OBS')
fsu_utc = dblarr(3)
fsu_utc[0] = double(strmid(fsu_utc_dum,11,2))	;UTC hrs
fsu_utc[1] = double(strmid(fsu_utc_dum,14,2))	;UTC min
fsu_utc[2] = double(strmid(fsu_utc_dum,17,strlen(fsu_utc_dum)-17))	;UTC sec

;recalculate the start MJD more accurate
mjd_start = floor(fsu_mjd) + fsu_utc[0]/24.d0 + fsu_utc[1]/(24.d0*60.d0) + fsu_utc[2]/(24.d0*3600.d0)

;calculate time of MIDI frames in usec with respect of begin of FSU fringe recording
mo = midi_time
midi_time = (midi_time-mjd_start)*86400.d0*1.D6

midi_time_orig = midi_time

openw, lun, workingdir+strcompress(star,/rem)+'.time', /get_lun
for i=0L,n_elements(midi_time)-1 do begin

  printf, lun, midi_time[i], format='(f21.11)'

endfor
close, lun
free_lun, lun


;remove MIDI frames that were taken before FSU fringe recording started
index = where(midi_time ge 0.)
if (index[0] ne -1) then begin
  midi_time = midi_time[index]
endif

;remove MIDI frames that were taken after FSU fringe recording ended
index = where(midi_time le max(fsu_time))
if (index[0] ne -1) then begin
  midi_time = midi_time[index]
endif

if (n_elements(midi_time) lt 2) then begin

  print, ' '
  print, '************************************************************'
  print, 'MIDI data were not recorded with FSUA data at the same time.'
  print, '************************************************************'
  print, ' '
endif

;to be safe remove the 1st MIDI frame because FSU fringe recording could have started right in the middle of the integration
;midi_time = midi_time[1:*]

;time difference between two consecutive MIDI frames with respect to the previous one, i.e. time span for every frame
diff = ts_diff(midi_time,1)*(-1.d0)
diff[1:n_elements(diff)-1] = diff[0:n_elements(diff)-2]
diff[0] = median(diff[1:*])


;assuming that 1st integrating of e.g. 18ms, than overhead due to read-out of e.g. 2ms -> time span 20ms
overhead = diff - midi_inttime


;filter out FSU frames that are different in time from MIDI frames
;check which FSU frames are taken during one MIDI frame, i.e. in 18ms of MIDI integration time
;1st calculate begin and end of MIDI exposure
midi_end = midi_time - overhead
midi_start = midi_time - diff


;get OPDC state for fsu_time because time of OPDC state is twice the frequency
index_fsu_time = closest(fsu_time_double,fsu_time)
fsu_state = fsu_state_orig[index_fsu_time]

fsu_gd_ave = dblarr(n_elements(midi_time))
fsu_gd_orig = dblarr(n_elements(midi_time),50)

for i=0L,n_elements(midi_time)-1 do begin

  ;check which FSU frames belong to one MIDI frame in time domain
  index1 = 0
  index1 = where(fsu_time ge midi_start[i] and fsu_time le midi_end[i])

  ;check if ALL FSU frames have OPDC state 7 inside one MIDI frame
  ;1st selection criteria
  index2 = 0
  index2 = where(fsu_state[index1] eq 7 or fsu_state[index1] eq 5)
  if (n_elements(index2) ne n_elements(index1)) then begin	;FSU lost fringes during one MIDI integration
 
    midi_time[i] = -99999 

  endif else begin

    fsu_gd_ave[i] = mean(fsu_gd[index1])
    fsu_gd_orig[i,0:n_elements(index1)-1] = fsu_gd[index1]

  endelse

endfor

index = 0
index = where(midi_time ne -99999)
if (index[0] ne -1) then begin

  midi_time = midi_time[index]
  fsu_gd_ave = fsu_gd_ave[index]
  fsu_gd_orig = fsu_gd_orig[index,*]

endif

;don't know why but deleting the last element, 'last' GD value is close to zero
midi_time = midi_time[0:n_elements(midi_time)-2]
fsu_gd_ave = fsu_gd_ave[0:n_elements(fsu_gd_ave)-2]
fsu_gd_orig = fsu_gd_orig[0:n_elements(fsu_gd_orig[*,0])-2,*]


flag = uintarr(n_elements(midi_time_orig))
midi_gd_fsu = dblarr(n_elements(midi_time_orig))
;do flagging and write file
for i=0L,n_elements(midi_time)-1 do begin

  idxflag = where(midi_time[i] eq midi_time_orig)
  if (idxflag[0] ne -1) then begin
    flag[idxflag] = 1
    midi_gd_fsu[idxflag] = fsu_gd_ave[i]
  endif

endfor


; idxval = where(fsu_gd_orig ne 0.)
; openw, lun, strcompress(star,/rem)+'.fsugd_orig',/get_lun
; for i=0L,n_elements(idxval)-1 do begin
;   printf, lun, fsu_gd_orig[idxval[i]]
; endfor
; close, lun
; free_lun, lun


openw, lun, workingdir+strcompress(star,/rem)+'.fsuflags', /get_lun
for i=0L,n_elements(midi_time_orig)-1 do begin
  printf, lun, flag[i], format='(I1)'
endfor
close, lun
free_lun, lun

; ;pure GD from FSU-A
; openw, lun, workingdir+strcompress(star,/rem)+'.fsugd', /get_lun
; for i=0L,n_elements(midi_time_orig)-1 do begin
;   printf, lun, midi_gd_fsu[i]
; endfor
; close, lun
; free_lun, lun


end


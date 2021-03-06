@counter.pro
@fifteenb.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@readfits.pro
@sxpar.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@get_eso_keyword.pro
@stddev.pro
@moment.pro
@unwrap_disp.pro
@uniq.pro
@ts_diff.pro
@reverse.pro
@setdifference.pro
@setintersection.pro


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

pro FSU_FringeTrackPerformance_wrappedPD

;**********************************************************************************************

obsnum = ''
read, 'Enter Observation ID: ', obsnum

calquest = ''
read, 'Skycal [1] or Labcal [2]: ', calquest

path = '/home/amueller/work/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/'
file_mkdir, path

fsupath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'

;file = file_search(fsupath+'PACMAN_OBS_GENERIC*.fits',count=n)
;file = file_search(fsupath+'PACMAN_GEN_RECORD*.fits',count=n)


if (calquest eq '2') then pathlabcal = fsupath+'LabCalibration/'
if (calquest eq '1') then pathlabcal = fsupath+'FringeScan/'

readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', cn, id, scical, nf, mf, photA, photB, pf, maskf, phot, format = 'a,a,a,d,a,a,a,a,a,a', /silent
idx = where(cn eq obsnum)
if (idx[0] ne -1) then begin
  id = id[idx]
  nf = nf[idx]
  scical = scical[idx]
  mf = mf[idx]
  pf = pf[idx]
  maskf = maskf[idx]
  phot = phot[idx]
  photA = photA[idx]
  photB = photB[idx]
endif else begin
  print, ''
  print, 'Commissioning number does not match with current data set in observation.txt.'
  print, ''
  return
endelse

nobs = n_elements(idx)
miditime = strmid(mf,16,24)

readcol, '/home/amueller/work/MIDI_FSUA/observation_FSUA.txt', comm, fid, kmag, fsu, bg, ff1, ff2, labcal, format='a,a,d,a,a,a,a,a', skipline=1, /silent

idx = where(comm eq obsnum)
if (idx[0] ne -1) then begin
  comm = comm[idx]
  fid = fid[idx]
  kmag = kmag[idx]
  fsu = fsu[idx]
  bg = bg[idx]
  ff1 = ff1[idx]
  ff2 = ff2[idx]
  labcal = labcal[idx]
endif

if (array_equal(id, fid) ne 1) then begin
  print, 'Check obsveration lists.'
  stop
endif

;**********************************************************************************************

openw, lun, path+'FSUA_FringeTrackPerformance_wrappedPD_'+obsnum+'.txt', width=1400,/get_lun;, /append
printf, lun, '## The RMS value of the phasedelay is the STDDEV of written PD values, i.e. computed on a wrapped phase!'
printf, lun, '                             PACMAN file                Star   Freq. [Hz]   LR St.7   St.5-7 [%]     rmsGD_7   rmsGD_5-7     rmsPD_7   rmsPD_5-7 [nm] IRIS DIT [ms]'
;printf, lun, '                             PACMAN file                Star   Freq. [Hz]   LR St.7   St.5-7 [%]       rmsGD       rmsPD [nanometer]'



for i=0,nobs-1 do begin
; for i=7,nobs-1 do begin
;for i=42,nobs-1 do begin

  counter, i+1, nobs, 'Processing observation '
;   print, 'Observation '+strcompress(i+1,/rem)+' / '+strcompress(nobs,/rem)

  file = file_search(fsupath+pf[i],count=nfiles)
  file = file[0]

  dummy = readfits(file,hdr,/silent)
  ;available DITs are 1,2,4
  irisdit = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS IAS IRIS DIT'))*1000.d0	;in ms
  freq = get_eso_keyword(hdr, 'HIERARCH ESO ISS PRI FSU1 FREQ')
  freq = round(double(freq))
;   dit = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PRI FSU1 DIT'))
;   dit = uint(round(dit*1000.d0))

;   if (obsnum eq '6' or obsnum eq '8' or obsnum eq '9' or obsnum eq '14') then $
;     st11 = mrdfits(file,11,/silent) else st11 = mrdfits(file,11,/silent)
  st11 = mrdfits(file, 'OPDC', /silent)

;   if (obsnum eq '6' or obsnum eq '8' or obsnum eq '9' or obsnum eq '14') then $
;     st7 = mrdfits(file,9,/silent) else st7 = mrdfits(file,7,/silent)
  st7 = mrdfits(file, 'IMAGING_DATA_FSUA', /silent)

  gd = st7.gd
  pd = st7.pd
  time = st7.time

  state_orig = st11.state
  time_double = st11.time

  index_time = closest(time_double,time)
  state = state_orig[index_time]


;unwrap phase from FSU
  caldir = pathlabcal+'LabCal_'+labcal[i]+'/'
  if (calquest eq '2') then caldir = pathlabcal+'LabCal_'+labcal[i]+'/'
  if (calquest eq '1') then caldir = pathlabcal

  if (calquest eq '2') then begin

    lc_file1 = 'FSUA_Calibration_ABCD.dat'
    lc_file2 = 'FSUA_Calibration_Waves.dat'
    lc_file3 = 'FSUA_Calibration_Lambda.dat'
    lc_file4 = 'FSUA_Calibration_Angles.dat'

  endif

  if (calquest eq '1') then begin

    lc_file1 = 'FSUA_SkyCalibration_ABCD.dat'
;     lc_file2 = 'FSUA_SkyCalibration_Waves.dat'
    lc_file3 = 'FSUA_SkyCalibration_Lambda.dat'
;     lc_file4 = 'FSUA_SkyCalibration_Angles.dat'

  endif

  readcol, caldir+lc_file3, wltmp, format='d', /silent
  wl_center = median(wltmp[0:3])

  nstate = n_elements(time)
  idx7 = where(state eq '7')
  if (idx7[0] eq -1) then begin
    lr7 = 0
  endif else begin
    lr7 = double(n_elements(idx7))/double(nstate)
  endelse

  idx57 = where(state eq '7' or state eq '5')
  if (idx57[0] eq -1) then begin
    lr57 = 0
  endif else begin
    lr57 = double(n_elements(idx57))/double(nstate)
  endelse


  rmsGD7 = stddev(gd[idx7])*1.d9
  rmsGD57 = stddev(gd[idx57])*1.d9
  rmsPD7 = stddev(pd[idx7])*wl_center*1.d9/(2.*!DPI)
  rmsPD57 = stddev(pd[idx57])*wl_center*1.d9/(2.*!DPI)



  printf, lun, pf[i], id[i], freq, lr7*100., lr57*100., rmsGD7, rmsGD57, rmsPD7, rmsPD57, irisdit, format='(a40,a20,i10,2f13.2,4f12.2,i13)'


endfor

close, lun
free_lun, lun

;stop
end
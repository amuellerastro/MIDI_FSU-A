@counter.pro
@fifteenb.pro
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
@ts_diff.pro
@uniq.pro
@reverse.pro
@stddev.pro
@fsc_color.pro
@tvread.pro
@setdifference.pro
@setintersection.pro
@unwrap_disp_juwe.pro
@unwrap_disp.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@mean.pro
@moment.pro
@mima.pro
@mima.pro
@uniq.pro
@ts_diff.pro
@reverse.pro
@numlines.pro
@repchr.pro
@is_ieee_big.pro
@ieee_to_host.pro




;we know that phase has to be constant and has to be around zero or at a constant value between -pi and +pi
;can we say that if fsu_state 1 or 5 its a random phase but not the fringe phase

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


pro calcMIDI_Koresko;, star, midi_file_gd, midi_file, fsu_file

;VAARIABLES
;**********************************************************************************************
obsnum = ''
read, 'Enter Observing Number: ', obsnum
;obsnum = 'P92CHOQUET'

calquest = ''
read, 'Skycal [1] or Labcal [2]: ', calquest

upquest = ''
read, 'Qi [1] or Juwe [2] phase unwrapping: ', upquest

readcol, '/home/amueller/work/MIDI_FSUA/observation_FSUA.txt', comm, id, kmag, fsu, bg, ff1, ff2, labcal, format='a,a,d,a,a,a,a,a', skipline=1, /silent

idx = where(comm eq obsnum)
if (idx[0] ne -1) then begin

  comm = comm[idx]
  id = id[idx]
  kmag = kmag[idx]
  fsu = fsu[idx]
  bg = bg[idx]
  ff1 = ff1[idx]
  ff2 = ff2[idx]
  labcal = labcal[idx]

endif else begin

  print, ''
  print, 'Observing number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

;paths
pathdata = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'
if (calquest eq '2') then pathlabcal = pathdata+'LabCalibration/'
if (calquest eq '1') then pathlabcal = pathdata+'FringeScan/'
resultpath1 = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_Koresko_StandardMask/'
file_mkdir, resultpath1
resultpath2 = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_Koresko_CustomMask/'
file_mkdir, resultpath2
intpath1 = resultpath1+'IntermediateResults/'
file_mkdir, intpath1
intpath2 = resultpath2+'IntermediateResults/'
file_mkdir, intpath2
workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir

;**********************************************************************************************

;read in observation.txt to get time stamp from MIDI files
readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', mcomm, mid, scical, nfl, midi, photA, photB, mfsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', skipline=1, /silent

idx = where(mcomm eq obsnum)
if (idx[0] ne -1) then begin
  mcomm = mcomm[idx]
  mid = mid[idx]
  scical = scical[idx]
  nfl = nfl[idx]
  midi = midi[idx]
  mfsu = mfsu[idx]
  maskname = maskname[idx]
  photA = photA[idx]
  photB = photB[idx]
endif


for i=0,n_elements(fsu)-1 do begin
;for i=2,2 do begin
;for i=3,3 do begin
;for i=0,0 do begin

  counter, i+1, n_elements(fsu), 'Processing observation '


  if (upquest eq '1') then boxvalue = 2001	;default
  if (upquest eq '2') then boxvalue = 11	;default

  ;get MIDI time stamp
  idx = where(fsu[i] eq mfsu)
  miditime = strmid(midi[idx], 16, 8)

;read in of LabCal data
;   caldir = pathlabcal+'LabCal_'+labcal[i]+'/'
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


;FSU
  print, ' '
  print, 'read in FSU file'
;   print, ' '
  fsu_file = pathdata+fsu[i]
  dum = readfits(fsu_file,fsu_hdr,/silent)


  d7 = mrdfits(fsu_file,'IMAGING_DATA_FSUA',/silent)
  ;post-processed FSU GD and PD
  fsu_time_orig = d7.time	;usec, this should be the reference
  fsu_gd_orig = d7.gd	;meter, FSU GD
  fsu_phase_orig = d7.pd

  fsu_time = fsu_time_orig
  fsu_phase = fsu_phase_orig
  fsu_gd = fsu_gd_orig



; http://www.davidpace.com/physics/graduate-school/idl-butterworth.htm
; create butterworth filter for phase and GD / adjust cutoff frequency?
  print, 'Butterworth filtering of FSU GD'
;   print, ' '
  ;compute cutoff freq for GD, maximum frequency to survive (Hz)
  maxfreqGD = 1.0d0
  acq = double(n_elements(fsu_time))/((double(fsu_time[n_elements(fsu_time)-1])-double(fsu_time[0]))/1.d6)	  ;determine frequency resolution of FFT 
  fftsize = n_elements(fsu_gd)	  ;number of data points to be passed into FFT function
  deltafreq = acq / fftsize
  cutfreqGD = fix(maxfreqGD/deltafreq)	;determine index value of cutoff frequency
  filterGD = butterworth(n_elements(fsu_gd),cutoff=cutfreqGD)	;generate Butterworth function with the desired cutoff
  fsu_gd_filtered = real_part(FFT(FFT(fsu_gd, -1, /double)*filterGD, 1, /double))	;filter data
 
  fsu_gd = fsu_gd_filtered
  fsu_gd_orig = fsu_gd_filtered


;   ;compute cutoff freq for PD, maximum frequency to survive (Hz)
;   maxfreqPD = 1.d3
;   acq = double(n_elements(fsu_time))/((double(fsu_time[n_elements(fsu_time)-1])-double(fsu_time[0]))/1.d6)	  ;determine frequency resolution of FFT 
;   fftsize = n_elements(fsu_phase)	  ;number of data points to be passed into FFT function
;   deltafreq = acq / fftsize	
;   cutfreqPD = fix(maxfreqPD/deltafreq)	;determine index value of cutoff frequency
;   filterPD = butterworth(n_elements(fsu_phase),cutoff=cutfreqPD)	;generate Butterworth function with the desired cutoff
;   fsu_pd_filtered = real_part(FFT(FFT(fsu_phase, -1)*filterPD, 1))	;filter data
; 
;   fsu_pd = fsu_pd_filtered
;   fsu_pd_orig = fsu_pd_filtered


;stop
;11th table
  d11 = mrdfits(fsu_file,'OPDC',/silent)
  fsu_time_double = d11.time	;usec, double frequency of fsu_time2
  fsu_state_orig = d11.state	;OPDC state
  fsu_rtoffset_orig = d11.rtoffset
  ;fsu_uwphase_orig = d13.uwphase

  idxdouble = closest(fsu_time_double,fsu_time)
  fsu_state = fsu_state_orig[idxdouble]
  ;fsu_rtoffset = fsu_rtoffset_orig[idxdouble]

  ;ndata_orig = n_elements(fsu_time_orig)
  ;ndata = ndata_orig

  quest = 'y'
  if (quest eq 'y') then begin

    repeat begin 

      print, 'Phase and Dispersion unwrapping'

      if (upquest eq '1') then $
	  unwrap_disp, id[i], miditime, workingdir, fsu_time, fsu_phase_orig, fsu_gd_orig, fsu_state, wl_center, /sav, /hak, box=boxvalue

      if (upquest eq '2') then $
	  unwrap_disp_juwe, id[i], miditime, workingdir, fsu_time, fsu_phase_orig, fsu_gd_orig, fsu_state, wl_center, /sav, /hak, box=boxvalue

      read, 'Happy with result (y/n)?: ', quest
      if (quest ne 'y') then read, 'Enter new box size: ', boxvalue
      if (boxvalue mod 2) ne 1 then boxvalue = boxvalue + 1

    endrep until (quest eq 'y')

  endif

  openw, lun, intpath1+'BoxSize_for_Unwrapping.txt', width=1400, /get_lun, /append
    printf, lun, id[i], miditime, boxvalue, format='(a20, a10, i8)'
  close, lun
  free_lun, lun
  openw, lun, intpath2+'BoxSize_for_Unwrapping.txt', width=1400, /get_lun, /append
    printf, lun, id[i], miditime, boxvalue, format='(a20, a10, i8)'
  close, lun
  free_lun, lun
  if (upquest eq '2') then spawn, 'tiff2ps '+workingdir+id[i]+'_'+miditime+'_unwrap.tif > '+workingdir+id[i]+'_'+miditime+'_unwrap.ps'
  ;spawn, 'gv '+workingdir+id[i]+'_'+miditime+'_unwrap.ps'
  spawn, 'cp '+workingdir+id[i]+'_'+miditime+'_unwrap.ps '+intpath1
  spawn, 'mv '+workingdir+id[i]+'_'+miditime+'_unwrap.ps '+intpath2
  if (upquest eq '2') then spawn, 'rm '+workingdir+id[i]+'_'+miditime+'_unwrap.tif '

  restore, workingdir+id[i]+'_'+miditime+'_unwrap.sav';, /verbose

  pdgd = disp



;*******************
;*     Korekso     *
;*******************


; now execute calculations from Koresko
  ;eqn. 9
  ;dispersion_K = (pdgd-median(pdgd))*100.d0
  dispersion_K = pdgd*100.d0

  ; eqn. 11
  gamma_K = -1.095290252D24	;cm^-3

  sigma = gamma_K*dispersion_K	;cm^-2


  ; eqn 12
  ;refract_N = 6.816D-24	;cm^3
  ;from Koresko2006, p4, table 1
  refract_N = 6.491D-24	;cm^3
  refract_K = 9.029D-24	;cm^3

  fsu_phase_uw = pdel	;in meter
  midi_phasedelay = fsu_phase_uw*100.d0 + (refract_N-refract_K)*sigma
  midi_phasedelay = midi_phasedelay/100.d0	;in micron
;  midi_phasedelay = midi_phasedelay*2.d0*!DPI/10.34d-6	;<- central wavelength MIDI
  ;see http://www.eso.org/sci/facilities/paranal/instruments/midi/inst/filters.html
;  midi_phasedelay = midi_phasedelay/!dtor	;in degree for MIDI



  ; eqn. 13
  gamma_N = -1.36600825D23	;computed using wavelength range 8-13micron and central wavelength of 10.5micron. Practically its 10.34 micron.
  psi = refract_N - refract_K - (1.d0/gamma_N) + (1.d0/gamma_K)

;   midi_gd_calc = (((fsu_gd-median(fsu_gd))*100.d0)  + (psi*sigma))/100.d0	;meter
  midi_gd_calc = ((fsu_gd*100.d0)  + (psi*sigma))/100.d0	;meter

;   window, 5, xs=1500, ys=400
;   plot, fsu_time, midi_gd_calc, xst=1


  openw, lun, workingdir+id[i]+'_'+miditime+'.KoreskoGD_orig', width=1400, /get_lun
    for j=0L,n_elements(midi_gd_calc)-1 do begin
      printf, lun, midi_gd_calc[j]
    endfor
  close, lun
  free_lun, lun


  openw, lun, workingdir+id[i]+'_'+miditime+'.KoreskoPD_orig', width=1400, /get_lun
    for j=0L,n_elements(midi_phasedelay)-1 do begin
      printf, lun, midi_phasedelay[j]
    endfor
  close, lun
  free_lun, lun

  ;compute dispersion for MIDI in degree and wrap it
  disp = midi_phasedelay-midi_gd_calc	;in micrometer
  disp_orig = disp
  disp = disp*2.d0*!DPI/10.34d-6	;<- central wavelength MIDI, in radian
  ;see http://www.eso.org/sci/facilities/paranal/instruments/midi/inst/filters.html
  ;convert in degree and wrap it for MIDI phase file
  disp = disp/!dtor
; stop
;   for j=0L,n_elements(disp)-1 do begin
;     if (disp[j] gt !DPI) then disp[j] = disp[j]-2.*!DPI
;     if (disp[j] lt -!DPI) then disp[j] = disp[j]+2.*!DPI
;   endfor

  openw, lun, workingdir+id[i]+'_'+miditime+'.KoreskoDISP_orig', width=1400, /get_lun
    for j=0L,n_elements(disp)-1 do begin
      printf, lun, disp[j]
    endfor
  close, lun
  free_lun, lun

  ;move all files to both result directories
  spawn, 'cp '+workingdir+id[i]+'_'+miditime+'* '+resultpath1
  spawn, 'mv '+workingdir+id[i]+'_'+miditime+'* '+resultpath2

;stop
endfor

spawn, 'rm -r '+workingdir

;stop
end
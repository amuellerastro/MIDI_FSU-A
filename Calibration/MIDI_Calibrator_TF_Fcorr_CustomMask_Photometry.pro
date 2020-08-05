@readcol.pro
@numlines.pro
@remchar.pro
@gettok.pro
@repchr.pro
@strsplit.pro
@strnumber.pro
@valid_num.pro
@fxpar.pro
@uniq.pro
@fxbtform.pro
@fxbfind.pro
@ieee_to_host.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@mrd_skip.pro
@mrd_struct.pro
@is_ieee_big.pro
@interpol.pro
@readfits.pro
@sxpar.pro
@mean.pro
@moment.pro
@counter.pro
@fifteenb.pro
@ts_diff.pro
@caldat.pro
@stddev.pro
@oploterror.pro
@match.pro
@setdefaultvalue.pro   
@tag_exist.pro
@cgquery.pro
@cgplot.pro
@setdecomposedstate.pro
@decomposedcolor.pro
@cgdefaultcolor.pro
@getdecomposedstate.pro
@colorsareidentical.pro
@cgcolor.pro
@cgdefcharsize.pro
@fsc_color.pro
@mima.pro
@get_eso_keyword.pro
@get_baseline.pro
@get_visibility.pro


pro MIDI_Calibrator_TF_Fcorr_CustomMask_Photometry
;FOR CALIBRATORS
;**********************************************************************************************
;==========================================

if (obsnum eq 'P92CHOQUET') then begin

	  ;HD76111 HD80934
  diameter = [0.411d0, 0.467d0]	;mas
  base = ['U1U3','U1U4']	;used baselines for this program

endif

;==========================================


;ewsversion = 'MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013Feb02'
;ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'

obsnum = ''
read, 'Enter Observation ID: ', obsnum
;obsnum = 'P92CHOQUET'

basequest = ''
for i=0,n_elements(base)-1 do print, strcompress(i,/rem), '  ', base[i]
read, 'Which baseline?: ', basequest


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
;   photA = photA[idx]
;   photB = photB[idx]

endif else begin

  print, ''
  print, 'Commissioning number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

time = strmid(midi,16,24)

idxsc = where(scical eq 'Cal')
 if (idxsc[0] ne -1) then begin
  comm = comm[idxsc]
  id = id[idxsc]
  scical = scical[idxsc]
  nfl = nfl[idxsc]
  midi = midi[idxsc]
  fsu = fsu[idxsc]
  maskname = maskname[idxsc]
  time = time[idxsc]
  n = n_elements(idxsc)
endif else begin
  print, 'No calibrator targets observed in that run.'
  stop
endelse


;do not include the 4th calibrator (HD1014 03:30:57) bad one...
if (obsnum eq 'P88') then begin

  idx = where(time ne '03:30:57')
  comm = comm[idx]
  id = id[idx]
  scical = scical[idx]
  nfl = nfl[idx]
  midi = midi[idx]
  fsu = fsu[idx]
  maskname = maskname[idx]
  time = time[idx]
  n = n_elements(idx)

endif



;select reduction method
redmethod = ''
print, ''
print, '     1 EWS'
print, '     2 CohInt'
print, '     3 Koresko'
read, 'For which reduction method?: ', redmethod
if (redmethod eq '1') then $
  resultsdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_EWS_CustomMask/'
if (redmethod eq '2') then $
  resultsdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_CohInt_CustomMask/'
if (redmethod eq '3') then $
  resultsdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_Koresko_CustomMask/'
print, ''

photdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/Photometry_CustomMask/'

midipath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/'
pcrdir = resultsdir+'Calibrator_PhotonCountRate/'
file_mkdir, pcrdir
workingdir = '/media/disk_MIDIFSU/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir



openw, lun, pcrdir+'PhotonCountRate'+'_'+base[basequest]+'.txt', width=1400, /get_lun;, /append
  printf, lun, 'Median photon count rate between 10.2 and 12 micrometer'
close, lun
free_lun, lun

;**********************************************************************************************

for i=0,n-1 do begin
;for i=0,0 do begin

  counter, i+1, n, 'Processing observation '

  if (id[i] eq 'HD76111') then diam = diameter[0]/1.d3
  if (id[i] eq 'HD80934') then diam = diameter[1]/1.d3
  diamerr = 0.

  flux12 = nfl[i]



;   model_Ftot_orig = transpose((*specData[0]).flux)	;Jy
;   model_wl_orig = transpose((*specData[0]).wavelength)	;meter
  ;plot, (*specData[0]).wavelength ,(*specData[0]).flux

  st = mrdfits(photdir+id[i]+'_'+time[i]+'.photometry.fits', 2, /silent)
  raw_Ftot = st[5].data1
  raw_Ftoterr = st[11].data1


;compute projected baseline
  bl_info = get_baseline(midipath, midi[i])
  bl = bl_info.bl	;meter
  pa = bl_info.padeg	;deg


;get rawFcorr and wavelength from corresponding .corr.fits file
  st1 = mrdfits(resultsdir+id[i]+'_'+time[i]+'.corr.fits', 1, /silent)	;wavelength information
  st3 = mrdfits(resultsdir+id[i]+'_'+time[i]+'.corr.fits', 3, /silent)	;raw corr. flux information

  wl = st1.eff_wave	;meter
  raw_Fcorr = st3.visamp
  raw_Fcorrerr = st3.visamperr

;intrpolate model flux with respect to measured wavelength
;   model_Ftot = interpol(model_Ftot_orig, model_wl_orig, wl, /spline)

;compute model Visibility for each wavelength
;   model_Vis = get_visibility(bl, wl, diam, diamerr)

;from the definition V = Fcorr / Ftot we can compute theoretical Fcorr
;   model_Fcorr = model_Vis*raw_Ftot

;photon count rate, ADU/s/Jy
  pcr = raw_Ftot/flux12
  pcrerr = raw_Ftoterr/flux12
;   pcr = raw_Fcorr/model_Fcorr

  openw, lun, workingdir+id[i]+'_'+time[i]+'.pcr', width=1400, /get_lun
    wlout = wl*1.d6
    printf, lun, 'photon count rate in ADU/s/Jy as a function of wavelength'
    for j=0,n_elements(wl)-1 do begin
      printf, lun, wlout[j], pcr[j], pcrerr[j], format="(f10.5,2f20.10)"
    endfor

  close, lun
  free_lun, lun


  ;get JD
  dum = readfits(resultsdir+id[i]+'_'+time[i]+'.corr.fits', hdr, /silent)
  jd = 0.5d0+double(get_eso_keyword(hdr, 'MJD-OBS'))

  idx = where(wl gt 10.2d-6 and wl lt 12.d-6)
  median_pcr = mean(pcr[idx])
  openw, lun, pcrdir+'PhotonCountRate'+'_'+base[basequest]+'.txt', width=1400, /get_lun, /append
    printf, lun, id[i], midi[i], jd, median_pcr, format="(a20, a28, f18.8, f12.6)"
  close, lun
  free_lun, lun

;stop

  spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'* '+resultsdir


endfor

spawn, 'rm -r '+workingdir


;plot average photon count rate as function of time or the transfer function of Fcorr as a function of wavelength and of night

readcol, pcrdir+'PhotonCountRate'+'_'+base[basequest]+'.txt', calid, midifile, jd, pcr, format='a,a,d,d', /silent


wlall = dblarr(n,n_elements(wl)) & pcrall = wlall & pcrallerr = pcrall
;window, 0, xs=1000, ys=800
for i=0,n-1 do begin
  readcol, resultsdir+id[i]+'_'+time[i]+'.pcr', wltmp, pcrtmp, pcrtmperr, format='d,d,d', /silent
  wlall[i,*] = wltmp
  pcrall[i,*] = pcrtmp
  pcrallerr[i,*] = pcrtmperr
;  if (i eq 0) then $
;    plot, wlall[i,*], pcrall[i,*], xr=[7.5,13], xst=1, yr=[0,140], yst=1 else $
;    oplot, wlall[i,*], pcrall[i,*]
endfor

;micron = string(181B) + 'm'

;Calculate the aspect ratio of display window.
 aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE 
 xsize = 8.0
  ysize = xsize * aspectRatio
  IF ysize GT 10.5 THEN BEGIN
    ysize = 10.5
    xsize = ysize / aspectRatio
  ENDIF
  ; Calculate the offsets, so the output window is not off the page.
  xoffset = (8.5 - xsize) / 2.0
  yoffset = (11.0 - ysize) / 2.0

!p.font=0
!p.thick=4
!x.thick=3
!y.thick=3

micron = '!Mm!X'+'m'


;identify single nights
diff = abs(ts_diff(jd,1))
idxnight = where(diff gt 0.5d0)

;if (idxnight[0] ne -1) then begin

  if (idxnight[0] ne -1) then nnights = n_elements(idxnight)+1
  if (idxnight[0] eq -1) then nnights = 1
  for i=0,nnights-1 do begin

    if (nnights gt 1) then begin	;if there are more than one night

      if (i eq 0) then begin
	jdnight = jd[0:idxnight[0]]
	wlnight = wlall[0:idxnight[0],*]
	pcrnight = pcrall[0:idxnight[0],*]
	pcrnighterr = pcrallerr[0:idxnight[0],*]
      endif

      if (i eq nnights-1) then begin
	jdnight = jd[idxnight[i-1]+1:*]
	wlnight = wlall[idxnight[i-1]+1:*,*]
	pcrnight = pcrall[idxnight[i-1]+1:*,*]
	pcrnighterr = pcrallerr[idxnight[i-1]+1:*,*]
      endif

      if (i ne 0 and i ne nnights-1) then begin
	jdnight = jd[idxnight[i-1]+1:idxnight[i]]
	wlnight = wlall[idxnight[i-1]+1:idxnight[i],*]
	pcrnight = pcrall[idxnight[i-1]+1:idxnight[i],*]
	pcrnighterr = pcrallerr[idxnight[i-1]+1:idxnight[i],*]
      endif

    endif else begin	;if there is only on night

      jdnight = jd
      wlnight = wlall
      pcrnight = pcrall
      pcrnighterr = pcrallerr

    endelse


    nobs = n_elements(jdnight)

    ;decide which "date" is the night, e.g. observations started on 25.11.2011 at 01:00 a.m. the date of the night is actually 24.11.2011
    jdtmp = jdnight	;just to be sure
    jdtmp = jdtmp[0]-0.5d0
    if ((jdtmp-floor(jdtmp)) gt 0.5d0) then begin
      caldat, 2400000.d0+jdtmp-1.d0, month, day, year, hour, minute, second
      night = strcompress([day, month, year], /rem)
    endif else begin
      caldat, 2400000.d0+jdtmp, month, day, year, hour, minute, second
      night = strcompress([day, month, year], /rem)
    endelse

    if (uint(night[1] lt 10)) then night[1] = '0'+strcompress(night[1],/rem)
    if (uint(night[0] lt 10)) then night[0] = '0'+strcompress(night[0],/rem)

    ;compute mean TF with error
    pcrave = dblarr(n_elements(wlnight[0,*]))
    pcrave_err = pcrave
    if (n_elements(wlnight[*,0]) gt 1) then begin
      for j=0,n_elements(wlnight[0,*])-1 do begin
	pcrave[j] = mean(pcrnight[*,j])
	pcrave_err[j] = stddev(pcrnight[*,j])
      endfor
    endif else begin
	stop
	pcrave = pcrnight[0,*]
	;USING 10% ERROR, DON'T KNOW HOW TO MAKE IT BETTER FOR THE MOMENT
	pcrave_err = 0.1*pcrnight[0,*]
    endelse

    openw, lun, pcrdir+night[2]+night[1]+night[0]+'_'+base[basequest]+'_meanTF_Fcorr_Calibrator.txt', width=1400, /get_lun
      for j=0,n_elements(pcrave)-1 do begin
	printf, lun, wlnight[0,j], pcrave[j], pcrave_err[j], format="(f10.5, 2f20.10)"
      endfor
    close, lun
    free_lun, lun


    set_plot, 'ps'
    device, isolatin=1
    device, filename=pcrdir+night[2]+night[1]+night[0]+'_'+base[basequest]+'_TF_Fcorr_Calibrator.ps', /color,XSIZE=23, YSIZE=17, XOffset=xoffset, YOffset=yoffset

;       yr = mima(pcrnight)
;       yr[0] = 0.
	yr = [0,2000]	;UTs
; 	yr = [0,130]	;ATs

      plot, /nodata, wlnight[0,*], pcrnight[0,*], color=fsc_color('black'), background=fsc_color('white'), $
	xtitle='Wavelength ['+micron+']', ytitle='TF [ADU/s/Jy]', $
	xr=[8,13], xst=1, yr=yr, yst=1, xminor=2, yminor=2, charsize=1.5, $
	title='Baseline: '+base[basequest]+' / Night: '+night[0]+'.'+night[1]+'.'+night[2]
; 	title='Transfer Function for F!Dcorr!N / Night: '+night[0]+'.'+night[1]+'.'+night[2] 
      for j=0,nobs-1 do oploterror, wlnight[j,*], pcrnight[j,*], pcrnighterr[j,*], color=fsc_color('powder blue')

      oploterror, wlnight[0,*], pcrave, pcrave_err, color=fsc_color('blue'), errcolor=fsc_color('blue'), /nohat


    device,/close
    set_plot,'x'

    spawn, 'gv '+pcrdir+night[2]+night[1]+night[0]+'_'+base[basequest]+'_TF_Fcorr_Calibrator.ps'

  endfor

;endif
!P.Font=0
!p.thick=1
!x.thick=1
!y.thick=1


;rename PhotonCountRate because date of night was determined later

spawn, 'mv '+pcrdir+'PhotonCountRate'+'_'+base[basequest]+'.txt '+pcrdir+night[2]+night[1]+night[0]+'_'+base[basequest]+'_PhotonCountRate.txt'


stop
end

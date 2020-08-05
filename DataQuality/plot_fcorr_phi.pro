pro plot_fcorr_phi

obsnum = ''
read, 'Obsnum?: ', obsnum

mainpath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/'
datadir = mainpath+'MIDIreduced_EWS_CustomMask/'
resultsdir = datadir+'Plots/'
file_mkdir, resultsdir
resultsdirraw = datadir+'Plots/uncalibrated/'
file_mkdir, resultsdirraw
resultsdircal = datadir+'Plots/calibrated/'
file_mkdir, resultsdircal

fileraw = file_search(datadir+'*.corr.fits', count=nraw)
filecal = file_search(datadir+'*.calcorr.fits', count=ncal)

;=======================================================================================
;raw correlated flux and dphase

for i=0,nraw-1 do begin

  st1 = mrdfits(fileraw[i], 1)	;wavelength
  st2 = mrdfits(fileraw[i], 3)	;corr flux and dphase

  tmp = strmid(fileraw[i], 27, strpos(fileraw[i], '.', /reverse_search)-32)
  pos = strpos(tmp, '_')
  id = strmid(tmp, 0, pos)
  time = strmid(tmp, pos+1, 8)


  wl = st1.eff_wave*1.d6
  fc = st2.visamp
  fce = st2.visamperr
  ph = st2.visphi
  phe = st2.visphierr

  idx = where(wl ge 8. and wl le 13.)
  yr = [0,ceil(max(fc[idx]+fce[idx]))]

  set_plot, 'ps'
  device, isolatin=1
  device, filename=resultsdirraw+id+'_'+time+'_raw_Fcorr.ps', /color,XSIZE=25, YSIZE=15, XOffset=xoffset, YOffset=yoffset

  micron = '!Mm!X'+'m'
  !p.font=0

  plot, /nodata, [0,0], [0,0], xr=[8,13], xst=1, yr=yr, yst=1, $
    title=id+' '+time, xtitle='Wavelength ['+micron+']', ytitle='Raw Correlated Flux [counts]'

  ;ozone region
  polyfill, [9.3,10,10,9.3],[0,0,yr[1],yr[1]] ,/data, color=fsc_color('light gray')
  plot, /noerase, /nodata, [0,0], [0,0], xr=[8,13], xst=1, yr=yr, yst=1

  oploterror, wl, smooth(fc,3), fce, /nohat, color=fsc_color('black'), errcolor=fsc_color('black')

  device,/close
  set_plot,'x'

  idx = where(wl ge 8. and wl le 13.)
  yr = [floor(min(ph[idx]-phe[idx])),ceil(max(ph[idx]+phe[idx]))]

  set_plot, 'ps'
  device, isolatin=1
  device, filename=resultsdirraw+id+'_'+time+'_raw_dPhase.ps', /color,XSIZE=25, YSIZE=15, XOffset=xoffset, YOffset=yoffset

  micron = '!Mm!X'+'m'
  !p.font=0

  plot, /nodata, [0,0], [0,0], xr=[8,13], xst=1, yr=yr, yst=1, $
    title=id+' '+time, xtitle='Wavelength ['+micron+']', ytitle='Raw Differential Phase [deg]'

  oploterror, wl, smooth(ph,3), phe, /nohat, color=fsc_color('black'), errcolor=fsc_color('black')

  device,/close
  set_plot,'x'

endfor

;=======================================================================================
;calibrated correlated flux and dphase

for i=0,ncal-1 do begin

  st1 = mrdfits(filecal[i], 1)	;wavelength
  st2 = mrdfits(filecal[i], 3)	;corr flux and dphase

  tmp = strmid(filecal[i], 27, strpos(filecal[i], '.', /reverse_search)-32)
  pos = strpos(tmp, '_')
  id = strmid(tmp, 0, pos)
  time = strmid(tmp, pos+1, 8)


  wl = st1.eff_wave*1.d6
  fc = st2.visamp
  fce = st2.visamperr
  ph = st2.visphi
  phe = st2.visphierr

  idx = where(wl ge 8. and wl le 13.)
  yr = [0,ceil(max(fc[idx]+fce[idx])*10.)/10.]

  set_plot, 'ps'
  device, isolatin=1
  device, filename=resultsdircal+id+'_'+time+'_cal_Fcorr.ps', /color,XSIZE=25, YSIZE=15, XOffset=xoffset, YOffset=yoffset

  micron = '!Mm!X'+'m'
  !p.font=0

  plot, /nodata, [0,0], [0,0], xr=[8,13], xst=1, yr=yr, yst=1, $
    title=id+' '+time, xtitle='Wavelength ['+micron+']', ytitle='Correlated Flux [Jy]'

  ;ozone region
  polyfill, [9.3,10,10,9.3],[0,0,yr[1],yr[1]] ,/data, color=fsc_color('light gray')
  plot, /noerase, /nodata, [0,0], [0,0], xr=[8,13], xst=1, yr=yr, yst=1

  oploterror, wl, smooth(fc,3), fce, /nohat, color=fsc_color('black'), errcolor=fsc_color('black')

  device,/close
  set_plot,'x'

  idx = where(wl ge 8. and wl le 13.)
  yr = [floor(min(ph[idx]-phe[idx])),ceil(max(ph[idx]+phe[idx]))]

  set_plot, 'ps'
  device, isolatin=1
  device, filename=resultsdircal+id+'_'+time+'_cal_dPhase.ps', /color,XSIZE=25, YSIZE=15, XOffset=xoffset, YOffset=yoffset

  micron = '!Mm!X'+'m'
  !p.font=0

  plot, /nodata, [0,0], [0,0], xr=[8,13], xst=1, yr=yr, yst=1, $
    title=id+' '+time, xtitle='Wavelength ['+micron+']', ytitle='Differential Phase [deg]'

  ;ozone region
  polyfill, [9.3,10,10,9.3],[!y.crange[0],!y.crange[0],!y.crange[1],!y.crange[1]] ,/data, color=fsc_color('light gray')
  plot, /noerase, /nodata, [0,0], [0,0], xr=[8,13], xst=1, yr=yr, yst=1

  oploterror, wl, smooth(ph,3), phe, /nohat, color=fsc_color('black'), errcolor=fsc_color('black')

  device,/close
  set_plot,'x'

endfor

stop
end
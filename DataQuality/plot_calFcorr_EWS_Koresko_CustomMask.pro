pro plot_calFcorr_EWS_Koresko_CustomMask

micron = string(181B) + 'm'
lambda = '!7'+string(107B)+'!X'
phi = '!Mf!X'
alpha = '!Ma!X'
beta = '!Mb!X'
gamma = '!Mg!X'
delta = '!Md!X'

commnum = ''
read, 'Enter Observation ID: ', commnum

mainpath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/'
resultpath = mainpath+'Plots_CalFcorr_EWS_Koresko_CustomMask/'
file_mkdir, resultpath

;data path
path1 = mainpath+'MIDIreduced_EWS_CustomMask/'
path2 = mainpath+'MIDIreduced_Koresko_CustomMask/'
; path3 = mainpath+'MIDIreduced_CohInt_CustomMask/'


obsfile = '/home/amueller/work/MIDI_FSUA/observation.txt'

readcol, obsfile, comm, id, scical, flux, midi, phota, photb, fsu, mask, phot, format='a,a,a,d,a,a,a,a,a,a', /silent


idx = where(comm eq commnum[0])

id = id[idx]
scical = scical[idx]
flux = flux[idx]
midi = midi[idx]
fsu = fsu[idx]

time = strmid(midi, 16, 8)

nobs = n_elements(id) ;should be 44 observations

corr1 = file_search(path1+'*.calcorr', count=nfiles)

if (nfiles ne nobs) then begin
  print, 'number of listet observations do not match the number of reduced files'
  stop
endif


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

for i=0,nobs-1 do begin
;for i=3,3 do begin

  set_plot, 'ps'
  device, isolatin=1
  device, filename=resultpath+id[i]+'_'+time[i]+'_CalFcorr_EWS_Koresko_CustomMask.ps', /color,XSIZE=17, YSIZE=10, XOffset=xoffset, YOffset=yoffset

;read in data in chronological order - extract time stamps
  ;wavelength
  readcol, path1+id[i]+'_'+time[i]+'.calcorr', wl, fcorr1, fcorrerr1, format='d,d,d'
  readcol, path2+id[i]+'_'+time[i]+'.calcorr', wl, fcorr2, fcorrerr2, format='d,d,d'
;   readcol, path3+id[i]+'_'+time+'.calcorr', wl, fcorr3, fcorrerr3, format='d,d,d'

  idx = where(wl ge 8. and wl le 13.)
  wl = wl[idx]
  fcorr1 = fcorr1[idx]
  fcorrerr1 = fcorrerr1[idx]
  fcorr2 = fcorr2[idx]
  fcorrerr2 = fcorrerr2[idx]
;   fcorr3 = fcorr3[idx]
;   fcorrerr3 = fcorrerr3[idx]

  ;yr = mima([fcorr1+fcorrerr1,fcorr2+fcorrerr2,fcorr3+fcorrerr3])
  ;yr = mima([fcorr1+fcorrerr1,fcorr2+fcorrerr2])
  yr = mima([(fcorr2+fcorrerr2)+(fcorr2+fcorrerr2)*0.1])
  yr[0] = 0.


  plot, /nodata, wl, fcorr1, xr=[8.0,13], xst=1, xticklen=0.03, xminor=2, yminor=2, yr=yr, yst=1, charsize=1.5, color=fsc_color('black'), background=fsc_color('white'), $
  xtitle='Wavelength ['+micron+']', ytitle='F!Dcorr!N [Jy]', $
  pos=[0.12,0.165,0.97,0.97]

  polyfill, [9.3,10,10,9.3],[0,0,!y.crange[1],!y.crange[1]] ,/data, color=fsc_color('light gray')

  plot, /noerase, wl, fcorr1, xr=[8.0,13], xst=1, xticklen=0.03, xminor=2, yminor=2, yr=yr, yst=1, charsize=1.5, color=fsc_color('black'), background=fsc_color('white'), $
  xtitle='Wavelength ['+micron+']', ytitle='F!Dcorr!N [Jy]', $
  pos=[0.12,0.165,0.97,0.97]

  oploterror, wl, fcorr1, fcorrerr1, color=fsc_color('dark green'), errcolor=fsc_color('pale green'), /nohat
  oploterror, wl, fcorr2, fcorrerr2, color=fsc_color('red'), errcolor=fsc_color('pink'), /nohat

  oplot, wl, fcorr1, color=fsc_color('dark green')
  oplot, wl, fcorr2, color=fsc_color('red')

;   oploterror, wl, fcorr3, fcorrerr3, color=fsc_color('red'), errcolor=fsc_color('red'), /nohat

  if (id[i] eq 'WDSJ05320-0018A') then id[i] = delta+'OriA'
  if (id[i] eq 'betaPic') then id[i] = beta+'Pic'

  time[i] = strmid(time[i], 0, 5)
  legend, ['#'+strcompress(i+1,/rem)+' '+id[i]+' '+time[i]], /right, textcolors=fsc_color('black'), box=0, margin=0, charsize=1.5

  device,/close
  set_plot,'x'

;  spawn, 'gv '+resultpath+id[i]+'_'+time[i]+'_CalFcorr_EWS_Koresko_CustomMask.ps'

endfor

!p.multi = [0,1,0]


stop
end
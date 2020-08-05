;is v_theo^2 white light true? difference if no bandwidth smearing?
;is there a better way to fit/remove the bias power

;surface, waveac, opdAC, wavenumberAC, charsize=2, /zlog

;FILTER OUT FRINGE SIGNAL TO GET PHOTOMETRY by computing A+C and B+D
;using Merand2006 approach to normalize fringe signal by photometry

@interpol.pro
@multiplot.pro
@cgplot.pro
@colorsareidentical.pro
@cgtext.pro
@stddev.pro
@boot_mean.pro
@linspace.pro
@readcol.pro
@remchar.pro
@strnumber.pro
@mpfit2dpeak.pro
@mpfit2dfun.pro
@setdecomposedstate.pro
@cgcolorbar.pro
@setdefaultvalue.pro
@cgdefcharsize.pro
@congrid.pro
@cgdefaultcolor.pro
@cgloadct.pro
@exponent.pro
@histoplot.pro
@convert_to_type.pro
@getdecomposedstate.pro
@decomposedcolor.pro
@cgGreek.pro
@binsizeFSU.pro
@sigfig.pro
@cgcolor.pro
@tvread.pro
@legend.pro
@get_visibility.pro
@get_baseline.pro
@wiener.pro
@readfits.pro
@sxpar.pro
@gettok.pro
@get_eso_keyword.pro
@strsplit.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@closest.pro
@ts_diff.pro
@mpfitfun.pro
@mpfit.pro
@proceeding_text.pro
@wavelet.pro
@reverse.pro
@moment.pro
@chisqr_cvf.pro
@chisqr_pdf.pro
@igamma.pro
@bisect_pdf.pro
@resistant_mean.pro
@mean.pro
@poly.pro
@polyfitw.pro
@real_part.pro
@loadct.pro
@filepath.pro
@path_sep.pro
@cgcolor.pro
@cgsnapshot.pro
@rgb.pro
@hak.pro
@integral.pro
@showsym.pro
@legend.pro
@poly_fit.pro
@fsc_color.pro
@int_tabulated.pro
@uniq.pro

;=================================================================================================

;Here, A0 is the height of the Gaussian, A1 is the center of the Gaussian, A2 is the width (the standard deviation) of the Gaussian, A3 is the constant term, A4 is the linear term, and A5 is the quadratic term. 
function fit_gauss_3terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)
  return, fit
end
function fit_gauss_4terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]
  return, fit
end
function fit_gauss_5terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)
  return, fit
end
function fit_gauss_6terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)+(p[5]*x^2.)
  return, fit
end
function fit_gauss_7terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)+(p[5]*x^2.)+(p[6]*x^3.)
  return, fit
end
function fit_gauss_8terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)+(p[5]*x^2.)+(p[6]*x^3.)+(p[7]*x^4.)
  return, fit
end
;=================================================================================================

function normalize, x

  norm = x
  idx = where(finite(x) eq 1)
  norm[idx] = (x[idx]-min(x[idx]))/(max(x[idx])-min(x[idx]))

  return, norm

end

;=================================================================================================

function func_biaspower, x, p

  bp = p[0]+p[1]*(sin(2.d0*!DPI*x/max(x)))^2.	;bias power based on LeBouquin+2011
  ga = p[2]*exp(-0.5d0*((x-p[3])/p[4])^2.)	;Gauss
  fit = bp+ga

  return, fit

end

;=================================================================================================

function func_offset, x, p

  fit = x*p[0]
  return, fit

end

;=================================================================================================

function snrfun, x, p
  fit = p[7]*( exp(-0.5d0*((x-p[0])/p[1])^2.) + p[4]*exp(-0.5d0*((x-p[5])/p[6])^2.)*sin(p[2]*x-p[3]) )+p[8]
  return, fit
end

;=================================================================================================

function loadFSU, path, file, pxmode

;read in FSU data
;   if (fsuid eq 'A') then obs = mrdfits(path+file,6, /silent)
;   if (fsuid eq 'B') then obs = mrdfits(path+file,8, /silent)

  if (strmatch(file, 'PACMAN*') eq 1) then obs = mrdfits(path+file,7, /silent) else $
       obs = mrdfits(path+file,6, /silent)

  obs2 = mrdfits(path+file, 11, /silent)

  idxdouble = closest(obs2.time, obs.time)

  time = obs.time
  snr = obs.pdsnr
  gd = obs.gd
  A = obs.data1	;white light is A[0], spectral channel 1,...,5 is A[1...5]
  B = obs.data2
  C = obs.data3
  D = obs.data4
  darkwhite = obs.data5
  fuoffset = obs2[idxdouble].fuoffset
  rtoffset = obs2[idxdouble].rtoffset
  state = obs2[idxdouble].state
  dlpos = obs2[idxdouble].dl4

  for i=0,5 do begin

    A[i,*] = A[i,*]-darkwhite[0,*]
    B[i,*] = B[i,*]-darkwhite[1,*]
    C[i,*] = C[i,*]-darkwhite[2,*]
    D[i,*] = D[i,*]-darkwhite[3,*]

  endfor

  ;remove non interesting states, where OPDC was not in 20 or 21
  idx = where(state eq 20 or state eq 21)
  time = time[idx]
  snr = snr[idx]
  gd = gd[idx]
  A = A[*,idx]
  B = B[*,idx]
  C = C[*,idx]
  D = D[*,idx]
  fuoffset = fuoffset[idx]
  rtoffset = rtoffset[idx]
  state = state[idx]
  dlpos = dlpos[idx]

  if (pxmode eq '3') then begin

    ;do not consider pixel 1 and 5 
    A = A[[0,2,3,4],*]
    B = B[[0,2,3,4],*]
    C = C[[0,2,3,4],*]
    D = D[[0,2,3,4],*]

    ;sum up pixel 2,3,4 to create artificial white light
    tmpA = A[0,*]
    tmpB = B[0,*]
    tmpC = C[0,*]
    tmpD = D[0,*]

    for i=0L,long(n_elements(tmpA))-1L do begin

      tmpA[i] = total(A[1:3,i])
      tmpB[i] = total(B[1:3,i])
      tmpC[i] = total(C[1:3,i])
      tmpD[i] = total(D[1:3,i])

    endfor

    A[0,*] = tmpA
    B[0,*] = tmpB
    C[0,*] = tmpC
    D[0,*] = tmpD

  endif



  data = {time:time, snr:snr, gd:gd, A:A, B:B, C:C, D:D, fuoffset:fuoffset, rtoffset:rtoffset, state:state, dlpos:dlpos}

  return, data

end

;=================================================================================================

function getFSUA_DarkFlat, pathdf, darkfile, flatfile1, flatfile2, pxmode

  if (strmatch(darkfile, 'PACMAN*') eq 1) then begin

    st_d1 = mrdfits(pathdf+darkfile,7, /silent)
    st_f1 = mrdfits(pathdf+flatfile1,7, /silent)
    st_f2 = mrdfits(pathdf+flatfile2,7, /silent)

  endif else begin

    st_d1 = mrdfits(pathdf+darkfile,6, /silent)
    st_f1 = mrdfits(pathdf+flatfile1,6, /silent)
    st_f2 = mrdfits(pathdf+flatfile2,6, /silent)

  endelse

  dark = dblarr(4,6)	;4 quadrants containing 6 pixels where 1st pixel is white light pixel
  flat1 = dblarr(4,6)	;flat beam1
  flat2 = dblarr(4,6)	;flat beam2

;   idxd1 = where(st_d1.time ne 0.d0)	;don't remember why implemented


;   fa1 = st_f1.data1[0]
;   fB1 = st_f1.data2[0]
;   fc1 = st_f1.data3[0]
;   fd1 = st_f1.data4[0]
;   fa2 = st_f2.data1[0]
;   fB2 = st_f2.data2[0]
;   fc2 = st_f2.data3[0]
;   fd2 = st_f2.data4[0]

;compute dark and flat

  for i=0,3 do begin	;loop over quadrants

    for j=0,5 do begin	;loop over channels

      if (i eq 0) then begin
	dark[i,j] = median(st_d1.data1[j]-st_d1.data5[0], /even)
	flat1[i,j] = median(st_f1.data1[j]-st_f1.data5[0], /even)-dark[i,j]
	flat2[i,j] = median(st_f2.data1[j]-st_f2.data5[0], /even)-dark[i,j]
      endif

      if (i eq 1) then begin 
	dark[i,j] = median(st_d1.data2[j]-st_d1.data5[1], /even)
	flat1[i,j] = median(st_f1.data2[j]-st_f1.data5[1], /even)-dark[i,j]
	flat2[i,j] = median(st_f2.data2[j]-st_f2.data5[1], /even)-dark[i,j]
      endif

      if (i eq 2) then begin 
	dark[i,j] = median(st_d1.data3[j]-st_d1.data5[2], /even)
	flat1[i,j] = median(st_f1.data3[j]-st_f1.data5[2], /even)-dark[i,j]
	flat2[i,j] = median(st_f2.data3[j]-st_f2.data5[2], /even)-dark[i,j]
      endif

      if (i eq 3) then begin
	dark[i,j] = median(st_d1.data4[j]-st_d1.data5[3], /even)
	flat1[i,j] = median(st_f1.data4[j]-st_f1.data5[3], /even)-dark[i,j]
	flat2[i,j] = median(st_f2.data4[j]-st_f2.data5[3], /even)-dark[i,j]
      endif

    endfor

  endfor

  if (pxmode eq '3') then begin

    dark = dark[*,[0,2,3,4]]
    flat1 = flat1[*,[0,2,3,4]]
    flat2 = flat2[*,[0,2,3,4]]

    ;sum up pixel 2,3,4 to create artificial white light
    tmpd = dark[*,0]
    tmpf1 = flat1[*,0]
    tmpf2 = flat2[*,0]

    for i=0,3 do begin

      tmpd[i] = total(dark[i,1:3])
      tmpf1[i] = total(flat1[i,1:3])
      tmpf2[i] = total(flat2[i,1:3])

    endfor

    dark[*,0] = tmpd
    flat1[*,0] = tmpf1
    flat2[*,0] = tmpf2

  endif


  data = {dark:dark, flat1:flat1, flat2:flat2}

  return, data

end

;=================================================================================================

; function func_kappa, x, p, _extra=struc
; 
;   fit = p[0]*struc.PA + p[1]*struc.PB
; 
;   return, fit
; 
; end

;=================================================================================================

function filter50, data, filter

  d = reform(data)
  tmpA = fft(d,-1)
  tmpA2 = tmpA*filter
  filtd = fft(tmpA2,1)

  return, filtd

end

;=================================================================================================


function check2filter, data, time, step

  ;before I find a smart way lets do it stupidly
  d = reform(data)
  ndata = n_elements(d)

  dfreq = 1.d0/(step*ndata)
  auxgrid = dindgen(ndata)
  wavenum = dfreq*auxgrid

  ;create artificial sinus with 50Hz frequency (=0.02s)
  t = reform(time)/1.d6
  s = sin(t*2.d0*!DPI/0.02d0)

  fs = abs(fft(s,-1))
  fs = fs[0:round(n_elements(fs)/2.)]

  ;find the position of the 50Hz peak in the FFT as iti is a function of ndata
  dum = max(fs, idxfreq)
  refwn = double(sigfig(wavenum[idxfreq],3))	;refeence wavenumber
;  refwn = wavenum[idxfreq]	;refeence wavenumber

  ;create filter
  filter = dblarr(ndata)
  filter[*] = 1.d0
  nwin = fix(max(wavenum)/refwn)	;total number of expected peaks
  nwin2 = ceil(double(nwin)/2.d0)
  tmpgrid = indgen(nwin)+1

  for i=0,nwin2-1 do begin

    idx1 = where(wavenum ge tmpgrid[i]*refwn-7.d3 and wavenum le tmpgrid[i]*refwn+7.d3)
    idx2 = ndata-idx1
    idx = [[idx1],[idx2]]
    filter[idx] = 0.d0

  endfor

  fd = abs(fft(d,-1))
  md = median(fd, /even)
  idxf = where(filter eq 0.d0)
  idx = where(fd[idxf] gt 3.*md)
  if (idx[0] ne -1) then flag = 1

;   window, 0, xs=1600, ys=600
;   plot, fd, yr=[0,50]
;   oplot, filter*50.d0, color=rgb(255,0,0)
;   plots, !x.crange, [md, md], color=rgb(0,255,0)
;   plots, !x.crange, [3.*md, 3.*md], color=rgb(255,255,0)

  return, [flag, filter]

end


;=================================================================================================

function comp_wavelet, x, y, plotflag

  ndata = n_elements(x)
  t = reform(x)
  dt = median(abs(ts_diff(reform(t),1)), /even)
  y = reform(y)
  wave = WAVELET(y, dt, PERIOD=period, COI=coi, /PAD, SIGNIF=signif, siglvl=0.95d0, s0=dt, dj=0.1, MOTHER='Morlet', scale=scale)
  wave = abs(wave)^2.
  nscale = n_elements(period)
  wavenumber = 1.d0/(100.d0*period) ;wavenumber cm^-1

  idxsort = sort(wavenumber)
  wavenumber = wavenumber[idxsort]
  wave = wave[*,idxsort]
  signifo = signif[idxsort]

  if (plotflag eq '1') then begin

    window, 0, xs=1000, ys=800
    !x.thick=2
    !y.thick=2

    SetDecomposedState, 0, CurrentState=state
;     !p.multi=[0,1,2]
    !p.multi=[0,1,0]
    ;levels = 12
    levels = 12
    step = (Max(alog10(wave)) - Min(alog10(wave))) / levels
    userLevels = IndGen(levels) * step + Min(alog10(wave))

    cgLoadCT, 33, NColors=12, Bottom=3
    contour, alog10(wave), t, wavenumber, $
      xtitle='OPD / [m]', ytitle='Wavenumber / [cm!U-1!N]', title='Noise Wavelet', $
      Background=cgColor('white'), Color=cgColor('black'), /fill, C_Colors=Indgen(levels)+3, level=userLevels, $ ;, NLevels=levels
      ystyle=8, ytickformat='exponent', xtickformat='exponent', xstyle=1, $;, xr=[-2.2d-4, 2.2d-4]
      yrange=[min(wavenumber),max(wavenumber)], /ytype, /yl, $
      charsize=2, position=[0.12, 0.1, 0.89, 0.95] ;,position=[0.12, 0.55, 0.88, 0.95]

    signif = rebin(transpose(signifo), ndata, nscale)
    contour, (wave/signif), t, wavenumber, /Overplot, levels=1, c_annot='95%', Color=cgColor('white')

    axis, yaxis=1, ylog=1, YRANGE=[max(period), min(period)], /save, ytitle='Wavelength / [m]', charsize=2, ytickformat='exponent', Color=cgColor('black'), yst=1

;     plots,t,coi,NOCLIP=0, color=cgcolor('black'), thick=2   ;*** anything "below" this line is dubious

    SetDecomposedState, state

    cgColorBar, NColors=12, Bottom=3, Divisions=8, $	;color bar
       Range=[round(min(wave)), round(max(wave))], Format='(f5.2)', $
       Position = [0.17, 0.92, 0.83, 0.93], AnnotateColor='black'

  ;----------------------------------------------------------------

  ;if ps file is required

;     set_plot, 'ps'
;     device, filename='wavelet.ps',/color,XSIZE=30, YSIZE=12, XOffset=xoffset, YOffset=yoffset
;     !p.thick=4
;     !x.thick=3
;     !y.thick=3
; 
;     wavenumberps = wavenumber
; 
;     idx = where(wavenumberps ge 1.9d3 and wavenumberps le 1.02d4)
;     wavenumberps = wavenumberps[idx]
;     periodps = 1./(100.*wavenumberps)
;     wavenumberps = wavenumberps/1.d3
;     wave = wave[*,idx]
; 
;       SetDecomposedState, 0, CurrentState=state
;       !p.multi=[0,1,0]
;       levels = 12
;       step = (Max(alog10(wave)) - Min(alog10(wave))) / levels
;       userLevels = IndGen(levels) * step + Min(alog10(wave))
; 
;       ticks = [2, 3, 4, 5, 6, 7, 8, 9]
;       nticks = n_elements(ticks)
;       wlticks = [1.d-6, 2.d-6, 3.d-6, 4.d-6, 5.d-6, 6.d-6, 7.d-6, 8.d-6, 9.d-6, 10.d-6]*1.d6
;       wlnticks = n_elements(wlticks)
; 
;       micron = '!4' + String("154B) + '!Xm'
; 
;       cgLoadCT, 33, NColors=12, Bottom=3
;       contour, alog10(wave), t, wavenumberps, $
; 	xtitle='OPD / [m]', ytitle='Wavenumber x10!U3!N / [cm!U-1!N]', $
; 	Background=cgColor('white'), Color=cgColor('black'), /fill, C_Colors=Indgen(levels)+3, level=userLevels, $
; 	ystyle=9, xtickformat='exponent', xstyle=1, /yl, $;, xr=[-2.2d-4, 2.2d-4]
; 	charsize=1.7, position=[0.08, 0.19, 0.92, 0.96], yticks=nticks-1, ytickv=ticks;, ytickformat='exponent'
; 
;       axis, yaxis=1, YRANGE=[max(periodps*1.d6), min(periodps*1.d6)], /save, ytitle='Wavelength / ['+micron+']', charsize=1.7, Color=cgColor('black'), yst=1, yticks=wlnticks-1, ytickv=wlticks
; 
;       SetDecomposedState, state
; 
; 
;     device,/close
;     set_plot,'x'
;     !p.thick=1
;     !x.thick=1
;     !y.thick=1
; 
;     spawn, 'gv wavelet.ps'

  ;----------------------------------------------------------------

    hak

  endif

  return, {wave:wave, opd:t, freq:wavenumber}

end



function get_uvw_bl, ha, latv, dec, t1x, t1y, t1z, t2x, t2y, t2z

  x = t1x - t2x
  y = t1y - t2y
  z = t1z - t2z

  latv  = latv/!radeg

  cha = cos(ha)
  sha = sin(ha)
  cdec = cos(dec)
  sdec = sin(dec)
  clat = cos(latv)
  slat = sin(latv)

  ;1st baseline T1 - T2

  ; compute u,v,w
  u = x*cha - y*slat*sha + z*clat*sha
  v = x*sha*sdec + y*(slat*cha*sdec + clat*cdec) - $
	  z*(clat*cha*sdec - slat*cdec)
  w = -x*sha*cdec - y*(slat*cha*cdec - clat*sdec) + $
		z*(clat*cha*cdec + slat*sdec)

  SinAlt =   (clat*cha*cdec + slat*sdec)
  alt = asin(SinAlt) ; <<<  in Radians!
  pa = atan(y, x) ; also in radians
  para = atan(sha, (slat/clat)*cdec-sdec*cha)

  alt = alt*!radeg ; to degrees

  bl = sqrt((u)^2.+(v)^2.)
  padeg = (atan(u, v))*!radeg

  if (padeg lt 0. and padeg ge -180.d0) then padeg = padeg + 180.d0; else stop

  return, [u,v,w,bl,padeg,alt]

end

;=================================================================================================


pro FSUA_CoherenceFactor

print, ''
path = '/home/amueller/work/MIDI_FSUA/Pipeline/V_FSUA/'
pathdf = '/home/amueller/work/MIDI_FSUA/Pipeline/V_FSUA/'	;but will be in SkyCalibration

;=========================================================================

; obsnum = ''
; read, 'Enter ObsID: ', obsnum

dsselect = ''
dataset = ['1: HD155826 - 20130706', $
           '2:    42Cet - 20131025', $
           '3:    24Psc - 20131028', $
           '4: VelaX1 U13 - 20140113', $
           '5: VelaX1 U14 - 20140113', $
           '6: P92Rivi - 20140306_1', $
           '7: P92Rivi - 20140306_2', $
           '8: P92Rivi - 20140306_3', $
           '9: P92Rivi - 20140307_1', $
          '10: P92Rivi - 20140307_2', $
          '11: P92Rivi - 20140307_3', $
          '12: P93KRAUSU14 - 20140412_1', $
          '13: P93KRAUSU14 - 20140412_2', $
          '14: P93KRAUSU14 - 20140412_3', $
          '15: P93KRAUSU14 - 20140412_4', $
          '16: P93MENU - 20140415', $
          '17: P93KRAUSU23 - 20140413_1', $
          '18: P93KRAUSU23 - 20140413_2', $
          '19: P93KRAUSU23 - 20140413', $
          '20: P93MONNIERA1C1 - 20140706', $
          '21: P93MONNIERB2C1 - 20140706', $
          '22: P94KOEHLERH0I1 - 20141023', $
          '23: P94KOEHLERD0I1 - 20141024']

print, 'Available Data Sets:'
for i=0,n_elements(dataset)-1 do print, dataset[i]
read, 'Select Data Set: ', dsselect

if (dsselect eq '1') then begin
  obsnum = 'FSUAscans_20130706'	;HD155826
  caldata = path+'20130706_HD155826/allFringeScans_HD155826.sav'
endif
if (dsselect eq '2') then begin
  obsnum = 'FSUAscans_20131025'  ;42Cet
  caldata = path+'20131025_42Cet/allFringeScans_42Cet.sav'
endif
if (dsselect eq '3') then begin
  obsnum = 'FSUAscans_20131028'	;24Psc
  caldata = path+'20131028_24Psc/allFringeScans_24Psc.sav'
endif
if (dsselect eq '4') then begin
  obsnum = 'FSUAscans_20140113_U13'	;24Psc
  caldata = path+'20140113_VelaX1/allFringeScans_VelaX1.sav'
endif
if (dsselect eq '5') then begin
  obsnum = 'FSUAscans_20140113_U14'	;VelaX1
  caldata = path+'20140113_VelaX1/allFringeScans_VelaX1.sav'
endif
if (dsselect eq '6') then begin
  obsnum = 'FSUAscans_20140306_1'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
endif
if (dsselect eq '7') then begin
  obsnum = 'FSUAscans_20140306_2'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
endif
if (dsselect eq '8') then begin
  obsnum = 'FSUAscans_20140306_3'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
endif
if (dsselect eq '9') then begin
  obsnum = 'FSUAscans_20140307_1'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
endif
if (dsselect eq '10') then begin
  obsnum = 'FSUAscans_20140307_2'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
endif
if (dsselect eq '11') then begin
  obsnum = 'FSUAscans_20140307_3'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
endif
if (dsselect eq '12') then begin
  obsnum = 'FSUAscans_20140412_1'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
endif
if (dsselect eq '13') then begin
  obsnum = 'FSUAscans_20140412_2'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
endif
if (dsselect eq '14') then begin
  obsnum = 'FSUAscans_20140412_3'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
endif
if (dsselect eq '15') then begin
  obsnum = 'FSUAscans_20140412_4'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
endif
if (dsselect eq '16') then begin
  obsnum = 'FSUAscans_20140415'
  caldata = path+'20140415_P93MENU/allFringeScans_P93MENU.sav'
endif
if (dsselect eq '17') then begin
  obsnum = 'FSUAscans_20140413_1'
  caldata = path+'20140413_P93KRAUSU23/allFringeScans_P93KRAUSU23.sav'
endif
if (dsselect eq '18') then begin
  obsnum = 'FSUAscans_20140413_2'
  caldata = path+'20140413_P93KRAUSU23/allFringeScans_P93KRAUSU23.sav'
endif
if (dsselect eq '20') then begin
  obsnum = 'FSUAscans_20140706_A1C1'
  caldata = path+'20140706_P93MONNIER/allFringeScans_P93MONNIER.sav'
endif
if (dsselect eq '21') then begin
  obsnum = 'FSUAscans_20140706_B2C1'
  caldata = path+'20140706_P93MONNIER/allFringeScans_P93MONNIER.sav'
endif
if (dsselect eq '22') then begin
  obsnum = 'FSUAscans_P94KOEHLERH0I1'
  caldata = path+'20141023_P94KOEHLERH0I1/allFringeScans_P94KOEHLERH0I1.sav'
endif
if (dsselect eq '23') then begin
  obsnum = 'FSUAscans_P94KOEHLERD0I1'
  caldata = path+'20141024_P94KOEHLERD0I1/allFringeScans_P94KOEHLERD0I1.sav'
endif
;=========================================================================

pxmode = ''
print, ''
print, '3: 3 pixel'
print, '5: 5 pixel'
read, 'Select Pixel mode: ', pxmode
print, ''

if (pxmode ne '3' and pxmode ne '5') then begin

  print, ''
  print, 'Choice does not exist. Stop.'
  print, ''
  return

endif

if (pxmode eq '3') then npix = 4
if (pxmode eq '5') then npix = 6

readcol, path+'FSUA_FringeScans.txt', obsid, scan, bg, ff1, ff2, starflag, diam, format='a,a,a,a,a,a,d', /silent

idx = where(obsid eq obsnum)

if (idx[0] ne -1) then begin

  obsid = obsid[idx]
  scan = scan[idx]
  bg = bg[idx]
  ff1 = ff1[idx]
  ff2 = ff2[idx]
  starflag = starflag[idx]
  diam = diam[idx]

endif else begin

  print, 'No such observation ID found. Stop.'
  stop

endelse


;   wl_cal = 2.2037d-6	;wavelength from Sky or Lab Calibration
  ;ave_wl_cal = [2.2037946d-06, 2.0088231d-06, 2.1150856d-06, 2.2310914d-06, 2.3544076d-06, 2.4471105d-06]	;median combined wavlengths for the individual spectral pixel, calib file should be read in later from Sky or LabCal

  restore, caldata
  if (pxmode eq '5') then begin

    ave_wl_cal = [median(lam2[*,0],/even), median(lam2[*,1],/even), median(lam2[*,2],/even), median(lam2[*,3],/even), median(lam2[*,4],/even), median(lam2[*,5],/even)]
    bw = median(bandwidth[*,*,0],/even)	;WHITE LIGHT
    wl_cal = ave_wl_cal[0]	;needed for reference fringe

  endif

  if (pxmode eq '3') then begin

    ave_wl_cal = [median(lam2[*,0],/even), median(lam2[*,2],/even), median(lam2[*,3],/even), median(lam2[*,4],/even)]
    bw = median(bandwidth[*,*,0],/even)	;WHITE LIGHT
    wl_cal = ave_wl_cal[0]	;needed for reference fringe

  endif


nfiles = n_elements(obsid)

;=================================================================================================

for xx=0,nfiles-1 do begin
;for xx=4,nfiles-1 do begin
;for xx=4,4 do begin
; for xx=10,10 do begin

  file = path+obsid[xx]+'/'+scan[xx]+'.fits'
  darkfile =  'FSUA_SKY_DARK_'+bg[xx]+'.fits'
  flatfile1 = 'FSUA_SKY_FLAT_'+ff1[xx]+'.fits'
  flatfile2 = 'FSUA_SKY_FLAT_'+ff2[xx]+'.fits'


  ;remove full path in front of file
    pos1 = strpos(file, 'FSUA_SECOND_FRINGE')
    pos2 = strpos(file, '.fits')
    file = strmid(file, pos1, pos2-pos1+5)


  ;extract file ID
    pos1 = strpos(file, 'SCAN_')
    pos2 = strpos(file, '.fits')
    fileid = strmid(file, pos1+5, pos2-pos1-5)

  ;extract keywords present in OBJ_SCAN files
  dum = readfits(path+obsid[xx]+'/'+file, hdr, /silent)
    target = get_eso_keyword(hdr, 'HIERARCH ESO OBS TARG NAME')
    base = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN BASE'))	;center position
    channel = (get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN CHANNEL'))
    nslew = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN NSLEW'))	;number of slews
    samprate = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN SAMPRATE'));number of points/micron
  ;   stroke = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN STROKE'))	;meter
    mjd = double(get_eso_keyword(hdr, 'MJD-OBS'))
;     altitude = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS ALT'))
    azimut = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS AZ'))
    airmass_start = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AIRM START'))
    airmass_end = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AIRM END'))
    seeing_start = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI FWHM START'))
    seeing_end = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI FWHM END'))
    tau0_start = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI TAU0 START'))*1.d3	;in [msec]
    tau0_end = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI TAU0 END'))*1.d3	;in [msec]
    opl = abs(double(get_eso_keyword(hdr, 'HIERARCH ESO DEL DLT1 OPL START'))-double(get_eso_keyword(hdr, 'HIERARCH ESO DEL DLT2 OPL START')))
    ;bl = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PBL12 START'))
    object = get_eso_keyword(hdr, 'HIERARCH ESO OBS TARG NAME')
    station = get_eso_keyword(hdr, 'HIERARCH ESO ISS CONF STATION1')+'-'+get_eso_keyword(hdr, 'HIERARCH ESO ISS CONF STATION2')
    fsudit = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PRI FSU1 FREQ'))	;ms

    if (seeing_start ge 0. and seeing_end ge 0.) then seeing = (seeing_start+seeing_end)/2.d0
    if (airmass_start ge 0. and airmass_end ge 0.) then airmass = (airmass_start+airmass_end)/2.
    if (tau0_start ge 0. and tau0_end ge 0.) then tau0 = (tau0_start+tau0_end)/2.

    if (seeing_start le 0. and seeing_end ge 0.) then seeing = seeing_end 
    if (seeing_start ge 0. and seeing_end le 0.) then seeing = seeing_start
    if (seeing_start ge 0. and seeing_end ge 0.) then seeing = mean([seeing_start,seeing_end])
    if (seeing_start le 0. and seeing_end le 0.) then seeing = 999.

    if (airmass_start le 0. and airmass_end ge 0.) then airmass = airmass_end 
    if (airmass_start ge 0. and airmass_end le 0.) then airmass = airmass_start
    if (airmass_start ge 0. and airmass_end ge 0.) then airmass = mean([airmass_start,airmass_end])
    if (airmass_start le 0. and airmass_end le 0.) then airmass = 999.

    if (tau0_start le 0. and tau0_end ge 0.) then tau0 = tau0_end
    if (tau0_start ge 0. and tau0_end le 0.) then tau0 = tau0_start
    if (tau0_start ge 0. and tau0_end ge 0.) then tau0 = mean([tau0_start,tau0_end])
    if (tau0_start le 0. and tau0_end le 0.) then tau0 = 999.


    radeg = double(strtrim(sxpar(hdr,'RA'),2))	;degrees
    dec = double(strtrim(sxpar(hdr,'DEC'),2))	;degrees
    lstsec = double(strtrim(sxpar(hdr,'LST'),2))	;seconds

    ;get the middle of the observation
    tmid = double(get_eso_keyword(hdr, 'EXPTIME'))/2.d0
    lstsec = lstsec+tmid
  ;   if (lstsec gt 86400.d0) then lstsec = lstsec-86400.d0	;needed?
    utc = double((get_eso_keyword(hdr, 'UTC')+tmid)/3600.d0)

    t1x = double(get_eso_keyword(hdr,'HIERARCH ESO ISS CONF T1X'))
    t1y = double(get_eso_keyword(hdr,'HIERARCH ESO ISS CONF T1Y'))
    t1z = double(get_eso_keyword(hdr,'HIERARCH ESO ISS CONF T1Z'))
    t2x = double(get_eso_keyword(hdr,'HIERARCH ESO ISS CONF T2X'))
    t2y = double(get_eso_keyword(hdr,'HIERARCH ESO ISS CONF T2Y'))
    t2z = double(get_eso_keyword(hdr,'HIERARCH ESO ISS CONF T2Z'))
    lat_vlti = double(get_eso_keyword(hdr,'HIERARCH ESO ISS GEOLAT'))	;latitude VLTI
    pa_start = double(get_eso_keyword(hdr,'HIERARCH ESO ISS PBLA12 START'))
    pa_end = double(get_eso_keyword(hdr,'HIERARCH ESO ISS PBLA12 END'))
    pa_mid = (pa_start+pa_end)/2.d0

    dec = dec*!dtor
    lst = (lstsec)*15.d0/(3600.d0*!radeg)
    ra = radeg*!dtor
    ha = lst - ra

    ; if (ha lt 0.) then ha = ha + 360.0d0
  ;   if (ha gt 180.d0) then ha = ha - 360.0d0
  ;   if (ha[i] gt 2.*!DPI or ha[i] lt -2.*!DPI) then stop	;something to do here?
    if (ha lt -!DPI) then ha = ha+2.*!DPI
    if (ha gt !DPI) then ha = ha-2.*!DPI

    result = get_uvw_bl(ha, lat_vlti, dec, t1x, t1y, t1z, t2x, t2y, t2z)
    u = result[0]
    v = result[1]
    w = result[2]
    bl = result[3]
    padeg = result[4]
    altitude = result[5]
    if (padeg lt pa_mid-180.) then padeg = padeg+180.


  ;define path for results
;     resultpath = path+'Vis_'+fileid+'_'+target+'/'
    resultpath0 = path+'Vis_'+obsid[xx]+'_'+pxmode+'pixel/'
    resultpath = resultpath0+'/CoherenceFactor/'
    spawn, 'mkdir -p '+resultpath

  ; print, '     sampling       stroke    nslew         seeing             tau'
  ; print, round(samprate), round(stroke*1.d6), round(nslew), seeing, tau0

  ;=================================================================================================

  print, '> Processing file '+strcompress(xx+1,/rem)+' of '+strcompress(nfiles,/rem), ': ', target, ' ', starflag[xx]
  print, '> Read in Data'

  ;LoadFSUA
  rawdata = loadFSU(path+obsid[xx]+'/', file, pxmode)
  ;(real) dark already subtracted from A,B,C,D

  ;=================================================================================================

  ;LoadDarkFlat

  ;flat fields are dark corrected and corrected for 'column 5'
  datadf = getFSUA_DarkFlat(pathdf+obsid[xx]+'/', darkfile, flatfile1, flatfile2, pxmode)

  ;=================================================================================================

  ;compute kappa matrix

  kappa = dblarr(4,npix)

  for i=0,3 do begin

    for j=0,npix-1 do begin

      kappa[i,j] = datadf.flat1[i,j]/datadf.flat2[i,j]

    endfor

  endfor

  print, '> Kappa Matrix'

  ;=================================================================================================

  ;extract scans

  diff = ts_diff(rawdata.state,1)
  idx1 = where(diff eq -1.d0)
  idx2 = where(diff eq 1.d0)

  tmp1 = abs(ts_diff(idx1,1))
  tmp1 = tmp1[0:n_elements(tmp1)-2]
  tmp2 = abs(ts_diff(idx2,1))
  tmp2 = tmp2[0:n_elements(tmp2)-2]

  idx = where(tmp1 lt median(tmp1)-0.1*median(tmp1))
  if (idx[0] ne -1) then begin

    for i=0,n_elements(idx)-1 do begin
      idx1[idx[i]] = -99
      idx2[idx[i]] = -99
    endfor

    idx99 = where(idx1 ne -99)
    idx1 = idx1[idx99]
    idx99 = where(idx2 ne -99)
    idx2 = idx2[idx99]

  endif

  idxscan = [idx1,idx2]
  idxscan = idxscan[sort(idxscan)]

;   idxscan = idxscan[0:10]
;   print, ''
;   print, 'REMOVE LINE 790'
;   print, ''
  

  nscans = n_elements(idxscan)-1

  if (total(idxscan) eq -2.) then begin
    print, 'Error. No Fringe Scans Found.'
    return
  endif

  ;number of elements per scan
  nval = dblarr(nscans)	;number of values per scan
  for i=0,nscans-1 do begin
; 
    snr = rawdata.snr[idxscan[i]+1:idxscan[i+1]]
    rtoffset = rawdata.rtoffset[idxscan[i]+1:idxscan[i+1]]
    nval[i] = n_elements(snr)

  endfor

  time = dblarr(nscans, max(nval))
  snr = dblarr(nscans, max(nval))
  gd = dblarr(nscans, max(nval))
  rawA = dblarr(nscans, npix, max(nval))
  rawB = dblarr(nscans, npix, max(nval))
  rawC = dblarr(nscans, npix, max(nval))
  rawD = dblarr(nscans, npix, max(nval))
  fuoffset = dblarr(nscans, max(nval))
  rtoffset = dblarr(nscans, max(nval))
  state = dblarr(nscans, max(nval))

  for i=0,nscans-1 do begin

    idxval = idxscan[i+1]-idxscan[i]-1	;number of elements per scan is not neccessarily the same
    time[i, 0:idxval] = rawdata.time[idxscan[i]+1:idxscan[i+1]]
    snr[i, 0:idxval] = rawdata.snr[idxscan[i]+1:idxscan[i+1]]
    gd[i, 0:idxval] = rawdata.gd[idxscan[i]+1:idxscan[i+1]]
    fuoffset[i, 0:idxval] = rawdata.fuoffset[idxscan[i]+1:idxscan[i+1]]
    rtoffset[i, 0:idxval] = rawdata.rtoffset[idxscan[i]+1:idxscan[i+1]]
    state[i, 0:idxval] = rawdata.state[idxscan[i]+1:idxscan[i+1]]

    rawA[i, *, 0:idxval] = rawdata.A[*,idxscan[i]+1:idxscan[i+1]]
    rawB[i, *, 0:idxval] = rawdata.B[*,idxscan[i]+1:idxscan[i+1]]
    rawC[i, *, 0:idxval] = rawdata.C[*,idxscan[i]+1:idxscan[i+1]]
    rawD[i, *, 0:idxval] = rawdata.D[*,idxscan[i]+1:idxscan[i+1]]

  endfor

  stroke = round(abs(rtoffset[0,n_elements(rtoffset[0,*])-1]-rtoffset[0,0])*1.d6)
  if (stroke lt 200. or stroke gt 500.) then stroke = round(abs(rtoffset[1,n_elements(rtoffset[1,*])-1]-rtoffset[1,0])*1.d6)

  ;=================================================================================================

  ;shorten arrays because there can be zeros in the arrays because of unequal number of values of the individual scans

  idx = where(time eq 0.d0)
  if (idx[0] ne -1) then begin

    for i=0,nscans-1 do begin

      idx = where(time[i,*] eq 0.d0)
;       if (idx[0] ne -1 and n_elements(idx) lt 2) then begin
      if (idx[0] ne -1) then begin

	time = time[*,0:idx[0]-1]
	snr = snr[*,0:idx[0]-1]
	gd = gd[*,0:idx[0]-1]
	fuoffset = fuoffset[*,0:idx[0]-1]
	rtoffset = rtoffset[*,0:idx[0]-1]
	state = state[*,0:idx[0]-1]
	rawA = rawA[*,*,0:idx[0]-1]
	rawB = rawB[*,*,0:idx[0]-1]
	rawC = rawC[*,*,0:idx[0]-1]
	rawD = rawD[*,*,0:idx[0]-1]

      endif

    endfor

  endif

  ;=================================================================================================

  ;check if array has even number of elements - needed for WIENER filtering
  
  ndata = double(n_elements(rawA[0,0,*]))
  
  if ((ndata mod 2.d0) ne 0.d0) then begin
  
    time = time[*,0:ndata-2]
    state = state[*,0:ndata-2]
    snr = snr[*,0:ndata-2]
    gd = gd[*,0:ndata-2]
    fuoffset = fuoffset[*,0:ndata-2]
    rtoffset = rtoffset[*,0:ndata-2]
    state = state[*,0:ndata-2]
    rawA = rawA[*,*,0:ndata-2]
    rawB = rawB[*,*,0:ndata-2]
    rawC = rawC[*,*,0:ndata-2]
    rawD = rawD[*,*,0:ndata-2]
  
  endif


  ;=================================================================================================

  ;Normalization

  for i=0,nscans-1 do begin

    for j=0,npix-1 do begin

      ;subtract background from interferometric data and divide by flat

      rawA[i,j,*] = (rawA[i,j,*]-datadf.dark[0,j])/(datadf.flat1[0,j]+datadf.flat2[0,j])
      rawB[i,j,*] = (rawB[i,j,*]-datadf.dark[1,j])/(datadf.flat1[1,j]+datadf.flat2[1,j])
      rawC[i,j,*] = (rawC[i,j,*]-datadf.dark[2,j])/(datadf.flat1[2,j]+datadf.flat2[2,j])
      rawD[i,j,*] = (rawD[i,j,*]-datadf.dark[3,j])/(datadf.flat1[3,j]+datadf.flat2[3,j])

    endfor

    proceeding_text,loop=(nscans), i=i, prompt='> Sky and Dark correction        '+string(i+1,form='(I4)')

  endfor

  ;=================================================================================================

  ;recenter fringe using positions found in the WSP
;   for i=0,nscans-1 do rtoffset[i,*] = rtoffset[i,*]+sign*1.d0*base	;sign=OPDC sign

  ;=================================================================================================

  ;-------------------------------------------------------------------------------------------------

  ;Definition of Arrays and Parameters

  ndata = n_elements(rtoffset[0,*])
  photA1 = dblarr(nscans, npix, ndata) & photB1 = photA1 & photC1 = photA1 & photD1 = photA1
  photA2 = dblarr(nscans, npix, ndata) & photB2 = photA2 & photC2 = photA2 & photD2 = photA2

  A = dblarr(nscans, npix, ndata) & B=A & C=A & D=A
  filtphotA1 = A & filtphotB1 = A & filtphotC1 = A & filtphotD1 = A
  filtphotA2 = A & filtphotB2 = A & filtphotC2 = A & filtphotD2 = A

  diffAC = A & diffBD = A

  wavelet_pos = dblarr(nscans) & wavelet_pos_err = wavelet_pos

  wavelet_height = dblarr(nscans)

  wavelet_lambda = dblarr(nscans) & wavelet_lambda_err = wavelet_lambda
  wavelet_lambda_width = dblarr(nscans) & wavelet_pos_width = dblarr(nscans)

  qualflag = intarr(nscans)
  qualflag[*] = 1
  qualflagscan = qualflag	;flag if too many scans were rejected because zpd offset was not set proper

  cohAC = dblarr(nscans,npix)
  cohBD = dblarr(nscans,npix)
  bw_cohAC = dblarr(npix) & bw_cohBD = bw_cohAC	;binwidth for histogram

  intrange_wavenum = [2.d3, 10.d3]	;integration limits over wavenumber

  if (pxmode eq '5') then fringeposlimit = 77.d0	;accepted minimum distance of fringe packet from the edge of the scan, if too close scan will be rejected
  if (pxmode eq '3') then fringeposlimit = 50.d0
  ;77 limit comes from the coherence length of the 5th pixel, wl~2.43d-6m and bw~1.55d-7m
  ;50 limit comes from the coherence length of the 4th pixel, wl~2.37d-6m and bw~1.47d-7m
  ;-------------------------------------------------------------------------------------------------

  ;create theoretical whitle light fringe of Visibility 1

  ;intensity = 522400.d0 	;using 1-0.4776
  ;x = reform(rtoffset[0,*])
  dx = median(ts_diff(reform(rtoffset[0,*]),1))
  x = linspace(-1.d0*stroke/1.d6/2.d0, stroke/1.d6/2.d0, abs(round(stroke/1.d6/dx)))
  lcoh = (wl_cal^2.-(bw^2.)/4.d0)/bw

  k0 = 2.d0*!DPI/wl_cal
  sinc = (sin(!DPI*(x)/lcoh))/(!DPI*(x)/lcoh)
  idxnan = where(finite(sinc) ne 1)	;in case there is a NAN when fringe going through 0 OPD
  if (idxnan[0] ne -1) then sinc[idxnan] = 1.d0
  fringe_theo = sinc*sin(k0*x)*2.d0

  wpsd = comp_wavelet(x, fringe_theo, '0')
  wave = wpsd.wave	;already abs(wave)^2
  opd = wpsd.opd
  wavenumber = wpsd.freq

  idx1 = where(wavenumber gt 2.d3 and wavenumber lt 8.d3)
  wavenumber = wavenumber[idx1]
  wave = wave[*,idx1]
  idx2 = where(opd lt 40.d-6 and opd gt -40.d-6)
  opd = opd[idx2]
  wave = wave[idx2,*]

  sum1 = dblarr(n_elements(opd))
  sum2 = dblarr(n_elements(wavenumber))
  for j=0,n_elements(opd)-1 do sum1[j] = total(wave[j,*])
  for j=0,n_elements(wavenumber)-1 do sum2[j] = total(wave[*,j])

  startval = [30.d0, 0.d0, 1.d-5, 0.d0]
  gauss_nterms = '4'
  dummy_err = dblarr(n_elements(opd))
  dummy_err[*] = 1.d0
  fitparams_gauss1theo = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', opd, sum1, dummy_err, startval, weights=sum1, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit1theo, perror=perror1, dof=dof, /quiet)

  wavelet_pos_theo = abs(fitparams_gauss1theo[1])
  wavelet_pos_width_theo = abs(2.d0*sqrt(2.d0*alog(2.d0))*fitparams_gauss1theo[2])

  startval = [50.d0, 1./(wl_cal*1.d2), 5.5d2, 0.d0]
  gauss_nterms = '4'
  dummy_err = dblarr(n_elements(wavenumber))
  dummy_err[*] = 1.d0
  fitparams_gauss2theo = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', wavenumber, sum2, dummy_err, startval, weights=sum2, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit2theo, perror=perror2, dof=dof, /quiet)

  wavelet_lambda_theo = abs(fitparams_gauss2theo[1])
  wavelet_lambda_width_theo = abs(2.d0*sqrt(2.d0*alog(2.d0))*fitparams_gauss2theo[2])

  plotflag = '0'
  if (plotflag eq '1') then begin
    window, 3, xs=900, ys=600
    !p.multi = [0,1,2]
    plot, opd, sum1, xtitle='OPD / [m]', xst=1, charsize=2
    oplot, opd, yfit1theo, color=cgcolor('red')
    plots, [0,0], !y.crange
    plot, wavenumber, sum2, xtitle='Wavenumber / [cm!U-1!N]', xst=1, charsize=2
    oplot, wavenumber, yfit2theo, color=cgcolor('red')
    !p.multi = [0,1,0]
    stop
  endif

  ;-------------------------------------------------------------------------------------------------

  for i=0,nscans-1 do begin

  ;-------------------------------------------------------------------------------------------------

  ;reverse scan direction to have monotonic increasing vectors

    ;the scans going in opposite direction have to be reversed
    ;state 20 gives increasing rtoffset, reverse scans with state eq 21
    if (state[i,0] eq 21.d0) then begin

      rtoffset[i,*] = reverse(transpose(rtoffset[i,*]))
  ;     grid[i,*] = reverse(transpose(grid[i,*]))

      for j=0,npix-1 do begin

	rawA[i,j,*] = reverse(transpose(rawA[i,j,*]))
	rawB[i,j,*] = reverse(transpose(rawB[i,j,*]))
	rawC[i,j,*] = reverse(transpose(rawC[i,j,*]))
	rawD[i,j,*] = reverse(transpose(rawD[i,j,*]))

      endfor

    endif

  ;-------------------------------------------------------------------------------------------------

  ;check for 50Hz noise, not needed anymore because A-C and B-D removes the peaks

;   step = abs(median(ts_diff(transpose(rtoffset[0,*]),1)))	;needed for wavenumbers
;   ;check if there is 50Hz noise in data
;   dum = max(max_snr, idxmax)
;   flagA = check2filter(rawA[idxmax,0,*], time[idxmax,*], step)
;   flagB = check2filter(rawB[idxmax,0,*], time[idxmax,*], step)
;   flagC = check2filter(rawC[idxmax,0,*], time[idxmax,*], step)
;   flagD = check2filter(rawD[idxmax,0,*], time[idxmax,*], step)
; 
;   if (flagA[0]+flagB[0]+flagC[0]+flagD[0] gt 0) then begin
; 
;     print, '> 50Hz noise present in data, filtering'
;     filter = flagA[1:ndata]
; 
;   endif

  ;-------------------------------------------------------------------------------------------------

  ;Photometry

    for j=0,npix-1 do begin

      ;because of flat and dark subtraction no offset fit needed
      rawphotAC = (rawA[i,j,*]+rawC[i,j,*])/4.d0
      rawphotBD = (rawB[i,j,*]+rawD[i,j,*])/4.d0

      photA1[i,j,*] = rawphotAC
      photA2[i,j,*] = rawphotAC/(kappa[0,j])
      photB1[i,j,*] = rawphotBD
      photB2[i,j,*] = rawphotBD/(kappa[1,j])
      photC1[i,j,*] = rawphotAC
      photC2[i,j,*] = rawphotAC/(kappa[2,j])
      photD1[i,j,*] = rawphotBD
      photD2[i,j,*] = rawphotBD/(kappa[3,j])


;       window, 0, xs=1300,ys=800
;       !p.multi=[0,2,2]
;       plot, rawA[i,j,*], title='A'
;         oplot, photA1[i,j,*], color=rgb(255,0,0)
;         oplot, photA2[i,j,*], color=rgb(255,255,0)
;         oplot, photA1[i,j,*]+kappa[0,j]*photA2[i,j,*], color=rgb(0,255,0)
;       plot, rawB[i,j,*], title='B'
;         oplot, photB1[i,j,*], color=rgb(255,0,0)
;         oplot, photB2[i,j,*], color=rgb(255,255,0)
;         oplot, photB1[i,j,*]+kappa[1,j]*photB2[i,j,*], color=rgb(0,255,0)
;       plot, rawC[i,j,*], title='C'
;         oplot, photC1[i,j,*], color=rgb(255,0,0)
;         oplot, photC2[i,j,*], color=rgb(255,255,0)
;         oplot, photC1[i,j,*]+kappa[2,j]*photC2[i,j,*], color=rgb(0,255,0)
;       plot, rawD[i,j,*], title='D'
;         oplot, photD1[i,j,*], color=rgb(255,0,0)
;         oplot, photD2[i,j,*], color=rgb(255,255,0)
;         oplot, photD1[i,j,*]+kappa[3,j]*photD2[i,j,*], color=rgb(0,255,0)
;       !p.multi=[0,1,0]
;       hak

  ;-------------------------------------------------------------------------------------------------

  ;Photometric Calibration

      ;Wiener filtering

;       tmpA1 = reform(photA1[i,j,*])
;       tmpB1 = reform(photB1[i,j,*])
;       tmpC1 = reform(photC1[i,j,*])
;       tmpD1 = reform(photD1[i,j,*])
; 
;       tmpA2 = reform(photA2[i,j,*])
;       tmpB2 = reform(photB2[i,j,*])
;       tmpC2 = reform(photC2[i,j,*])
;       tmpD2 = reform(photD2[i,j,*])
; 
; 
;       filterA1 = wiener(tmpA1, /quiet)
;       filterB1 = wiener(tmpB1, /quiet)
;       filterC1 = wiener(tmpC1, /quiet)
;       filterD1 = wiener(tmpD1, /quiet)
;       filterA2 = wiener(tmpA2, /quiet)
;       filterB2 = wiener(tmpB2, /quiet)
;       filterC2 = wiener(tmpC2, /quiet)
;       filterD2 = wiener(tmpD2, /quiet)
; 
;       fA1 = fft(tmpA1,-1)
;       fB1 = fft(tmpB1,-1)
;       fC1 = fft(tmpC1,-1)
;       fD1 = fft(tmpD1,-1)
; 
;       fA2 = fft(tmpA2,-1)
;       fB2 = fft(tmpB2,-1)
;       fC2 = fft(tmpC2,-1)
;       fD2 = fft(tmpD2,-1)
; 
;       filtphotA1[i,j,*] = real_part(fft(fA1*filterA1,+1))
;       filtphotB1[i,j,*] = real_part(fft(fB1*filterB1,+1))
;       filtphotC1[i,j,*] = real_part(fft(fC1*filterC1,+1))
;       filtphotD1[i,j,*] = real_part(fft(fD1*filterD1,+1))
; 
;       filtphotA2[i,j,*] = real_part(fft(fA2*filterA2,+1))
;       filtphotB2[i,j,*] = real_part(fft(fB2*filterB2,+1))
;       filtphotC2[i,j,*] = real_part(fft(fC2*filterC2,+1))
;       filtphotD2[i,j,*] = real_part(fft(fD2*filterD2,+1))


      ;filtphotX1/2  IS NOT NEEDED IF WE ARE USING THE DEFINITION FOR PHOTOMETRIC CALIBRATION OF MERAND2006


      ;Merand2006, SPIE, 6268
      A[i,j,*] = (1.d0/(sqrt(kappa[0,j]))) * $
		  (rawA[i,j,*]-photA1[i,j,*]-kappa[0,j]*photA2[i,j,*])/ $
		  sqrt(mean(photA1[i,j,*])*mean(photA2[i,j,*]))

      B[i,j,*] = (1.d0/(sqrt(kappa[1,j]))) * $
		  (rawB[i,j,*]-photB1[i,j,*]-kappa[1,j]*photB2[i,j,*])/ $
		  sqrt(mean(photB1[i,j,*])*mean(photB2[i,j,*]))

      C[i,j,*] = (1.d0/(sqrt(kappa[2,j]))) * $
		  (rawC[i,j,*]-photC1[i,j,*]-kappa[2,j]*photC2[i,j,*])/ $
		  sqrt(mean(photC1[i,j,*])*mean(photC2[i,j,*]))

      D[i,j,*] = (1.d0/(sqrt(kappa[3,j]))) * $
		  (rawD[i,j,*]-photD1[i,j,*]-kappa[3,j]*photD2[i,j,*])/ $
		  sqrt(mean(photD1[i,j,*])*mean(photD2[i,j,*]))


      ; 5.3. Alternative normalization methods
      ; If the SNR of the photometric channels PA and PB reaches too
      ; low values over the scan length, we choose to normalize the
      ; interferograms simply by averaging P over the fringe length,
      ; instead of using the Wiener filtered signal. This allows us to
      ; significantly reduce the amplification of the noise due to the
      ; normalization division. The limit between the two regimes is
      ; usually set to 5 times the readout noise. For interferograms that
      ; present very low photometric signal over the fringe packet itself,
      ; we discard the scan as a significant bias can be expected on
      ; the modulated power. Both averaging and Wiener filtering are
      ; almost equivalent on the final calibrated interferograms, with
      ; a slight advantage to the Wiener filtering when

      ;account for very low fluxes which lead to inifinte values by division of photometry
      idxA = where(finite(A[i,j,*]) ne 1)
      if (idxA[0] ne -1) then A[i,j,*] = (1.d0/(sqrt(kappa[0,j]))) * $
		  (rawA[i,j,*]-photA1[i,j,*]-kappa[0,j]*photA2[i,j,*])/ $
		  (total(photA1[i,j,*]+photA2[i,j,*])/double(ndata))

      idxB = where(finite(B[i,j,*]) ne 1)
      if (idxB[0] ne -1) then B[i,j,*] = (1.d0/(sqrt(kappa[1,j]))) * $
		  (rawB[i,j,*]-photB1[i,j,*]-kappa[1,j]*photB2[i,j,*])/ $
		  (total(photB1[i,j,*]+photB2[i,j,*])/double(ndata))

      idxC = where(finite(C[i,j,*]) ne 1)
      if (idxC[0] ne -1) then C[i,j,*] = (1.d0/(sqrt(kappa[2,j]))) * $
		  (rawC[i,j,*]-photC1[i,j,*]-kappa[2,j]*photC2[i,j,*])/ $
		  (total(photC1[i,j,*]+photC2[i,j,*])/double(ndata))

      idxD = where(finite(D[i,j,*]) ne 1)
      if (idxD[0] ne -1) then D[i,j,*] = (1.d0/(sqrt(kappa[3,j]))) * $
		  (rawD[i,j,*]-photD1[i,j,*]-kappa[3,j]*photD2[i,j,*])/ $
		  (total(photD1[i,j,*]+photD2[i,j,*])/double(ndata))


;       diffAC[i,j,*] = (A[i,j,*]-C[i,j,*])/2.d0
;       diffBD[i,j,*] = (B[i,j,*]-D[i,j,*])/2.d0

      diffac[i,j,*] = a[i,j,*]
      diffbd[i,j,*] = b[i,j,*]


      ; diffac und A the same
      ; stop
    ;   window, 0, xs=1000, ys=700
    ;   !p.multi=[0,2,2]
    ;   plot, a[i,j,*], xst=1, title='AC'
    ;     oplot, c[i,j,*], color=rgb(0,255,0)
    ;     plots, !x.crange, [0,0], color=rgb(0,0,255)
    ;   plot, diffAC[i,j,*], title='AC', xst=1
    ;     plots, !x.crange, [0,0], color=rgb(0,0,255)
    ;   plot, b[i,j,*], xst=1, title='BD'
    ;     oplot, d[i,j,*], color=rgb(0,255,0)
    ;     plots, !x.crange, [0,0], color=rgb(0,0,255)
    ;   plot, diffBD[i,j,*], title='BD', xst=1
    ;     plots, !x.crange, [0,0], color=rgb(0,0,255)
    ;   !p.multi=[0,1,0]
    ;   hak

      ;------------------------------------------------------------------------------------------------

      ;WAVELET to get fringe position and quality check

      if (j eq 0) then begin	;consider white light pixel only

	wpsdAC = comp_wavelet(rtoffset[i,*], (diffAC[i,j,*]), '0')
	waveAC = wpsdAC.wave	;already abs(wave)^2
	opdAC = wpsdAC.opd
	wavenumberAC = wpsdAC.freq

  ;       wpsdBD = comp_wavelet(rtoffset[i,*], (diffBD[i,j,*]), '0')
  ;       waveBD = wpsdBD.wave	;already abs(wave)^2
  ;       opdBD = wpsdBD.opd
  ;       wavenumberBD = wpsdBD.freq

      ;------------------------------------------------

	idxAC = where(wavenumberAC ge intrange_wavenum[0] and wavenumberAC le intrange_wavenum[1])
  ;       idxBD = where(wavenumberBD ge intrange_wavenum[0] and wavenumberBD le intrange_wavenum[1])

      ;------------------------------------------------

	;fringe position

	waveAC_tmp = waveAC[*,idxAC]
	sum1 = dblarr(n_elements(opdAC))
	sum2 = dblarr(n_elements(wavenumberAC[idxAC]))
	for k=0,n_elements(opdAC)-1 do sum1[k] = total(waveAC_tmp[k,*])
	for k=0,n_elements(wavenumberAC[idxAC])-1 do sum2[k] = total(waveAC_tmp[*,k])

	;fit gaussian to get ZPD offset for each scan
	maxsum = max(sum1, idxpos)
	startval = [maxsum, opdAC[idxpos], 1.d-5, 0.d0]
	gauss_nterms = '4'
	dummy_err = dblarr(n_elements(opdAC))
	dummy_err[*] = 1.d0
	pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)	;force amplitude to be positive
	pi[0].limited[0] = 1
	pi[0].limits[0] = 0.1d0
	fitparams_gauss = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', opdAC, sum1, dummy_err, startval, weights=sum1, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, parinfo=pi, /quiet)

	wavelet_pos[i] = fitparams_gauss[1]
	wavelet_pos_err[i] = perror[1]

	plotflag = '0'
	if (plotflag eq '1') then begin

	  window, 2, xs=600, ys=600
	  !p.multi = [0,1,2]
	  plot, wavenumberAC[idxAC], sum2, xtitle='Wavenumber / [cm!U-1!N]', xst=1
	  plot, opdAC, sum1, xtitle='OPD / [m]', xst=1
	  oplot, opdAC, yfit, color=cgcolor('red')
	  hak
	  !p.multi = [0,1,0]

	endif

	;-------------------------------------------------

	;if fringes are too close at the edge of a scan, reject them

	dum = max(diffac[i,0,*], idxfm)	;position of the maximum white light flux in the scan, i.e. fringe amplitude assuming photometric variations are corrected properly
	diff1 = (rtoffset[i,n_elements(rtoffset[i,*])-1]-rtoffset[i,idxfm])*1.d6
	diff2 = (abs(rtoffset[i,0]-rtoffset[i,idxfm]))*1.d6

	if (diff1 lt fringeposlimit or diff2 lt fringeposlimit) then begin

	  qualflag[i] = 0
	  qualflagscan[i] = 0

	endif

	;------------------------------------------------

	if (qualflag[i] eq 1) then begin

	  ;scan selection based on quality, i.e. fringe shape etc.

	  waveAC_tmp = waveAC[*,idxAC]
	  sum1 = dblarr(n_elements(opdAC))
	  sum2 = dblarr(n_elements(wavenumberAC[idxAC]))
	  for k=0,n_elements(opdAC)-1 do sum1[k] = total(waveAC_tmp[k,*])
	  for k=0,n_elements(wavenumberAC[idxAC])-1 do sum2[k] = total(waveAC_tmp[*,k])

	  maxsum = max(sum1, idxpos)
	  startval = [maxsum, opdAC[idxpos], 1.d-5, 0.d0]
	  ;startval = [2.d0, 0.d0, 1.d-5, 0.d0]
	  gauss_nterms = '4'
	  dummy_err = dblarr(n_elements(opdAC))
	  dummy_err[*] = 1.d0
	  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)	;force amplitude to be positive
	  pi[0].limited[0] = 1
	  pi[0].limits[0] = 0.01d0
	  fitparams_gauss1 = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', opdAC, sum1, dummy_err, startval, weights=sum1, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit1, perror=perror1, dof=dof, parinfo=pi, /quiet)

	  wavelet_height[i] = fitparams_gauss1[0]

	  wavelet_pos[i] = fitparams_gauss1[1]
	  wavelet_pos_err[i] = perror1[1]
	  wavelet_pos_width[i] = abs(2.d0*sqrt(2.d0*alog(2.d0))*fitparams_gauss1[2])

	  maxsum = max(sum2, idxpos)
	  startval = [maxsum, wavenumberAC[idxAC[idxpos]], 5.5d2, 0.d0]
	  gauss_nterms = '4'
	  dummy_err = dblarr(n_elements(wavenumberAC[idxAC]))
	  dummy_err[*] = 1.d0
	  fitparams_gauss2 = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', wavenumberAC[idxAC], sum2, dummy_err, startval, weights=sum2, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit2, perror=perror2, dof=dof, /quiet)

	  wavelet_lambda[i] = abs(fitparams_gauss2[1])
	  wavelet_lambda_err[i] = perror2[1]
	  wavelet_lambda_width[i] = abs(2.d0*sqrt(2.d0*alog(2.d0))*fitparams_gauss2[2])

	  plotflag = '0'
	  if (plotflag eq '1') then begin
	    window, 2, xs=900, ys=600
	    !p.multi = [0,1,2]
	    plot, opdAC, sum1, xtitle='OPD / [m]', xst=1, charsize=2
	    oplot, opdAC, yfit1, color=cgcolor('red')
	    plots, [0,0], !y.crange
	    plot, wavenumberAC[idxAC], sum2, xtitle='Wavenumber / [cm!U-1!N]', xst=1, charsize=2
	    oplot, wavenumberAC[idxAC], yfit2, color=cgcolor('red')
	    !p.multi = [0,1,0]
	    hak
	  endif

	  ; Kervella+2004
	  ; As the fringe packet has been recentered before the calibration, its position in the time domain is zero. Three parameters are then checked for quality:

	  if (wavelet_lambda[i] lt (1.d0/(wl_cal*1.d2)-0.3*1.d0/(wl_cal*1.d2)) or wavelet_lambda[i] gt (1.d0/(wl_cal*1.d2)+0.3*1.d0/(wl_cal*1.d2))) then tmp1 = 1 else tmp1 = 0	; -peak position in the frequency (i.e. wavenumber) domain (±30%);

	  if (wavelet_lambda_width[i] lt (wavelet_lambda_width_theo-0.4*wavelet_lambda_width_theo) or wavelet_lambda_width[i] gt (wavelet_lambda_width_theo+0.4*wavelet_lambda_width_theo)) then tmp2 = 1 else tmp2 = 0	; -peak width in the frequency (i.e. wavenumber) domain (±40%).

	  if (wavelet_pos_width[i] lt (wavelet_pos_width_theo-0.5*wavelet_pos_width_theo) or wavelet_pos_width[i] gt (wavelet_pos_width_theo+0.5*wavelet_pos_width_theo)) then tmp3 = 1 else tmp3 = 0	; -peak width in the time domain (OPD) (typically ±50% around the theoretical value is acceptable);

; 	  print, (1.d0/(wl_cal*1.d2)-0.3*1.d0/(wl_cal*1.d2)), wavelet_lambda[i], (1.d0/(wl_cal*1.d2)+0.3*1.d0/(wl_cal*1.d2))
; 	  print, (wavelet_lambda_width_theo-0.4*wavelet_lambda_width_theo), wavelet_lambda_width[i], (wavelet_lambda_width_theo+0.4*wavelet_lambda_width_theo)
; 	  print, (wavelet_pos_width_theo-0.5*wavelet_pos_width_theo), wavelet_pos_width[i], (wavelet_pos_width_theo+0.5*wavelet_pos_width_theo)
; 	  hak

	  if ((tmp1 + tmp2 + tmp3) gt 0) then qualflag[i] = 0

	endif

      endif

    endfor

  ;-------------------------------------------------------------------------------------------------

    ;recenter fringe using positions found in the WSP
    rtoffset[i,*] = rtoffset[i,*]-wavelet_pos[i]

  ;-------------------------------------------------------------------------------------------------

    proceeding_text,loop=(nscans), i=i, prompt='> Wavelet and Quality Check      '+string(i+1,form='(I4)')

  endfor


  ;=================================================================================================
  ;=================================================================================================
  ;=================================================================================================

  ;additional quality check based on height of the signal and wavelet position
  mpos = median(wavelet_pos)
  mheight = median(wavelet_height)
  for i=0,nscans-1 do begin

    if (wavelet_pos[i] lt (mpos-abs(0.5*mpos)) or wavelet_pos[i] gt (mpos+abs(0.5*mpos))) then qualflag[i] = 0
    if (wavelet_height[i] lt (mheight-0.3*mheight)) then qualflag[i] = 0

  endfor


  ;=================================================================================================
  ;=================================================================================================
  ;=================================================================================================


  ;=================================================================================================

  ;rescale white light fringe to get proper amplitude
  ;sum over pixel 1 to 5, rescale oroginal white light with the sum
  ;sum is not used becaused it has lower SNR

  factorAC = dblarr(nscans) & factorBD = dblarr(nscans)

  tmpAC = diffac[*,0,*] & tmpBD = diffBD[*,0,*]

  for i=0,nscans-1 do begin

    if (qualflag[i] eq 1) then begin

      for j=0,ndata-1 do begin

	if (pxmode eq '3') then begin
	  tmpAC[i,0,j] = total(diffAC[i,1:3,j])
	  tmpBD[i,0,j] = total(diffBD[i,1:3,j])
	endif

	if (pxmode eq '5') then begin
	  tmpAC[i,0,j] = total(diffAC[i,1:5,j])
	  tmpBD[i,0,j] = total(diffBD[i,1:5,j])
	endif

      endfor

      ;cut out central part of the fringe assuming ZPD offset is properly determined
      idx = where(rtoffset[i,*] ge -5.d-6 and rtoffset[i,*] le 5.d-6)

      startval = [5.d0]
      dummy_err = dblarr(ndata)
      dummy_err[*] = 1.d0

      factor = mpfitfun('func_offset', diffAC[i,0,idx], tmpAC[i,0,idx], dummy_err, startval, weights=abs(tmpAC[i,0,idx]), maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
      factorAC[i] = factor[0]

      factor = mpfitfun('func_offset', diffBD[i,0,idx], tmpBD[i,0,idx], dummy_err, startval, weights=abs(tmpBD[i,0,idx]), maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
      factorBD[i] = factor[0]

      diffAC[i,0,*] = diffAC[i,0,*]*factorAC[i]
      diffBD[i,0,*] = diffBD[i,0,*]*factorBD[i]

  ;     window, 0
  ;     plot, rtoffset[i,*], tmpAC[i,0,*], xr=[-2d-5,2d-5]
  ;     oplot, rtoffset[i,*], diffAC[i,0,*]*factor, color=rgb(255,0,0)
  ;     hak

    endif

    proceeding_text,loop=(nscans), i=i, prompt='> Scaling of White Light Pixel   '+string(i+1,form='(I4)')

  endfor

  openw, lun, resultpath+fileid+'_'+starflag[xx]+'_'+target+'_WhiteLightScalingFactor.txt', width=1400, /get_lun

    idx = where(qualflag eq 1)

    printf, lun, '           AC           BD'
    for i=0,n_elements(idx)-1 do printf, lun, factorAC[idx[i]], factorBD[idx[i]], format='(2f13.9)'

  close, lun
  free_lun, lun

  ;=================================================================================================

  ;in case fringe is close to the edge of the scanrange because of bad OPD model, piston etc. adjust intrange_opd

  intrange_opd = 100.d-6	;default value

  idx = where(qualflag eq 1)
  mintmp = dblarr(n_elements(idx)) & maxtmp = mintmp
  for i=0,n_elements(idx)-1 do begin

    mintmp[i] = min(rtoffset[idx[i],*])
    maxtmp[i] = max(rtoffset[idx[i],*])

  endfor

  minopd = max(mintmp)
  maxopd = min(maxtmp)

  if (minopd ge -1.*intrange_opd) then intrange_opd = abs(minopd)
  if (maxopd le intrange_opd) then intrange_opd = abs(maxopd)
  ;------------------------------------------------

  for i=0,nscans-1 do begin

    if (qualflag[i] eq 1) then begin

      idx = where(rtoffset[i,*] lt intrange_opd and rtoffset[i,*] gt -1.*intrange_opd)

      ;------------------------------------------------

      for j=0,npix-1 do begin

	;computation of wavelets AC / BD

	wpsdAC = comp_wavelet(rtoffset[i,idx], (diffAC[i,j,idx]), '0')
	waveAC = wpsdAC.wave	;already abs(wave)^2
	opdAC = wpsdAC.opd
	wavenumberAC = wpsdAC.freq

	wpsdBD = comp_wavelet(rtoffset[i,idx], (diffBD[i,j,idx]), '0')
	waveBD = wpsdBD.wave	;already abs(wave)^2
	opdBD = wpsdBD.opd
	wavenumberBD = wpsdBD.freq

	idxAC = where(wavenumberAC ge intrange_wavenum[0] and wavenumberAC le intrange_wavenum[1])
	idxBD = where(wavenumberBD ge intrange_wavenum[0] and wavenumberBD le intrange_wavenum[1])

	;------------------------------------------------


	;inegration A-C, B-D over time (OPD)

	np = n_elements(wavenumberAC)
	intACfreq = dblarr(np) & intBDfreq = dblarr(np)

	for k=0,np-1 do begin

	  intACfreq[k] = integral(reform(rtoffset[i,idx]), abs(waveAC[*,k]))
	  intBDfreq[k] = integral(reform(rtoffset[i,idx]), abs(waveBD[*,k]))
	  ;intACfreq[k] = tsum(reform(rtoffset[i,idx]), abs(waveAC[*,k]))
	  ;intBDfreq[k] = tsum(reform(rtoffset[i,idx]), abs(waveBD[*,k]))

	endfor

	;------------------------------------------------

; 	;remove bias power from AC
; 
	startval = [2., 4500.d0, 1200.d0, -9.5d0, 8.d-4, -7.d-8]
	pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
	pi[0].limited[0] = 1
	pi[0].limited[1] = 1
	pi[0].limits[0] = 0.d0
	pi[0].limits[1] = 1.d3
	pi[1].limited[0] = 1
	pi[1].limited[1] = 1
	pi[1].limits[0] = intrange_wavenum[0]
	pi[1].limits[1] = intrange_wavenum[1]
	pi[2].limited[0] = 1
	pi[2].limited[1] = 1
	pi[2].limits[0] = 100.d0
	pi[2].limits[1] = 3.d3

	dummy_err = dblarr(n_elements(idxAC))
	dummy_err[*] = 1.d0

	fitparams = mpfitfun('fit_gauss_6terms', wavenumberAC[idxAC], alog10(intACfreq[idxAC]), dummy_err, startval, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet, parinfo=pi)

	;window, 0
	;plot, wavenumberAC[idxAC], intACfreq[idxAC], charsize=2, /yl
	;oplot, wavenumberAC[idxAC], 10.^yfit, color=fsc_color('blue')
	;oplot, wavenumberAC[idxAC], 10.^(fitparams[3]+(fitparams[4]*wavenumberAC[idxAC])+(fitparams[5]*wavenumberAC[idxAC]^2.)), color=fsc_color('green')
	;oplot, wavenumberAC[idxAC], 10.^yfit - 10.^(fitparams[3]+(fitparams[4]*wavenumberAC[idxAC])+(fitparams[5]*wavenumberAC[idxAC]^2.)), color=fsc_color('red')
	;legend, ['noise', 'signal', 's+n'], box=0, margin=0, /right, /top, color=[fsc_color('green'), fsc_color('red'), fsc_color('blue')], psym=[sym(1), sym(1), sym(1)], charsize=2
 	;;window, 2
	;;surface, waveac, opdac, wavenumberac, charsize=5, az=30, ax=40, min_value=0.005
	;hak

	;intACfreq[idxAC] = 10.^yfit - 10.^(fitparams[3]+(fitparams[4]*wavenumberAC[idxAC])+(fitparams[5]*wavenumberAC[idxAC]^2.))	;gauss
 	intACfreq[idxAC] = intACfreq[idxAC] - 10.^(fitparams[3]+(fitparams[4]*wavenumberAC[idxAC])+(fitparams[5]*wavenumberAC[idxAC]^2.))	;gauss

	;------------------------------------------------

	;remove bias power from BD

	startval = [2., 4500.d0, 1200.d0, -9.5d0, 8.d-4, -7.d-8]
	dummy_err = dblarr(n_elements(idxBD))
	dummy_err[*] = 1.d0

	fitparams = mpfitfun('fit_gauss_6terms', wavenumberBD[idxBD], alog10(intBDfreq[idxBD]), dummy_err, startval, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet, parinfo=pi)

	;window, 1
	;plot, wavenumberBD[idxBD], intBDfreq[idxBD], charsize=2, /yl
	;oplot, wavenumberBD[idxBD], 10.^yfit, color=fsc_color('blue')
	;oplot, wavenumberBD[idxBD], 10.^(fitparams[3]+(fitparams[4]*wavenumberBD[idxBD])+(fitparams[5]*wavenumberBD[idxBD]^2.)), color=fsc_color('green')
	;;oplot, wavenumberBD[idxBD], intBDfreq[idxBD] - 10.^(fitparams[3]+(fitparams[4]*wavenumberBD[idxBD])+(fitparams[5]*wavenumberBD[idxBD]^2.)), color=fsc_color('red')
	;legend, ['noise', 'signal', 's+n'], box=0, margin=0, /right, /top, color=[fsc_color('green'), fsc_color('red'), fsc_color('blue')], psym=[sym(1), sym(1), sym(1)], charsize=2

	;intBDfreq[idxBD] = 10.^yfit - 10.^(fitparams[3]+(fitparams[4]*wavenumberBD[idxBD])+(fitparams[5]*wavenumberBD[idxBD]^2.))	;gauss
	intBDfreq[idxBD] = intBDfreq[idxBD] - 10.^(fitparams[3]+(fitparams[4]*wavenumberBD[idxBD])+(fitparams[5]*wavenumberBD[idxBD]^2.))

	;------------------------------------------------

	;integrate over frequency axis of WPSD

	cohAC[i,j] = integral(wavenumberAC[idxAC], intACfreq[idxAC])
	cohBD[i,j] = integral(wavenumberBD[idxBD], intBDfreq[idxBD])
	;cohAC[i,j] = tsum(wavenumberAC[idxAC], intACfreq[idxAC])
	;cohBD[i,j] = tsum(wavenumberBD[idxBD], intBDfreq[idxBD])

	;check if we get NAN values. Reason is that the polynomial fit did not converge and 10^x will result in infinity
	idxNAN_AC = where(finite(cohAC[i,j]) ne 1)
	if idxNAN_AC[0] ne -1 then stop
	idxNAN_BD = where(finite(cohBD[i,j]) ne 1)
	if idxNAN_BD[0] ne -1 then stop

      endfor

    endif

  proceeding_text,loop=(nscans), i=i, prompt='> Integration of WPSD            '+string(i+1,form='(I4)')

  endfor

  ;=================================================================================================

  ;rescale coherence factor
  for i=0,nscans-1 do begin

    cohAC[i,0] = cohAC[i,0]/factorAC[i]
    cohBD[i,0] = cohBD[i,0]/factorBD[i]

  endfor

  ;=================================================================================================

  ;discard scans which were found to be too close at the edge of a scan

  idxqual = where(qualflag eq 1)
  if (idxqual[0] eq -1) then begin

    print, ''
    print, 'All scans are too bad quality.'
    return

  endif else begin

    rtoffset = rtoffset[idxqual,*]
    A = A[idxqual,*,*]
    B = B[idxqual,*,*]
    C = C[idxqual,*,*]
    D = D[idxqual,*,*]
    diffAC = diffAC[idxqual,*,*]
    diffBD = diffBD[idxqual,*,*]
    cohAC = cohAC[idxqual,*]
    cohBD = cohBD[idxqual,*]
    wavelet_pos = wavelet_pos[idxqual,*]
    factorAC = factorAC[idxqual]
    factorBD = factorBD[idxqual]

    idxr0 = where(qualflagscan eq 0)
    idxr1 = where(qualflagscan eq 1)
    ratio = 100.d0*double(n_elements(idxr0))/double(n_elements(idxr1))

    print, '> Discarded '+strcompress(nscans-n_elements(idxqual),/rem)+' scans based on quality check.'
    if (ratio gt 50.) then print, '> !!! More than 50% of fringe packages are too close to the edge !!!'

    nscans = n_elements(idxqual)

  endelse

  ;=================================================================================================

  ;some histogram statistics
  skewness_cohAC = dblarr(npix) & kurtosis_cohAC = dblarr(npix)
  sd_cohAC = dblarr(npix) & mean_cohAC = dblarr(npix)
  bootsd_cohAC = dblarr(npix) & bootmean_cohAC = dblarr(npix)
  skewness_cohBD = dblarr(npix) & kurtosis_cohBD = dblarr(npix)
  sd_cohBD = dblarr(npix) & mean_cohBD = dblarr(npix)
  bootsd_cohBD = dblarr(npix) & bootmean_cohBD = dblarr(npix)

  for j=0,npix-1 do begin

    tmpAC = moment(cohAC[*,j])
    ;mean_cohAC[j] = tmpAC[0]
    ;mean_cohAC[j] = median(cohAC[*,j], /even)
    ;sd_cohAC[j] = sqrt(tmpAC[1])
    skewness_cohAC[j] = tmpAC[2]
    kurtosis_cohAC[j] = tmpAC[3]

    tmpBD = moment(cohBD[*,j])
    ;mean_cohBD[j] = tmpBD[0]
    ;mean_cohBD[j] = median(cohBD[*,j], /even)
    ;sd_cohBD[j] = sqrt(tmpBD[1])
    skewness_cohBD[j] = tmpBD[2]
    kurtosis_cohBD[j] = tmpBD[3]

    ;bootstrap mean
    nboot = 10000.d0
    boot_mean, cohAC[*,j], nboot, tmp
    bootmean_cohAC[j] = mean(tmp)
    bootsd_cohAC[j] = stddev(tmp)
    boot_mean, cohBD[*,j], nboot, tmp
    bootmean_cohBD[j] = mean(tmp)
    bootsd_cohBD[j] = stddev(tmp)

  endfor

  ;=================================================================================================

;   ;AVERAGING DOES NOT WORK, WASHES OUT THE FRINGE

;   ;Quality check of the remaining fringe scans w.r.t. the mean AC white light fringe package
; 
;   ;-------------------------------------------------------------------------------------------------
; 
;   ;before averaging fringe package select the scan with the most centered fringe at use this rtoffset as referecne
; 
;   rangeval = dblarr(nscans)
;   for i=0,nscans-1 do rangeval[i] = abs(rtoffset[i,0]+rtoffset[i,n_elements(rtoffset[i,*])-1])*1.d6
;   dum = min(rangeval, idxbest)	;fringes are centered well here
; 
;   ;-------------------------------------------------------------------------------------------------
; 
;   ;interpolate to compute mean fringe package
; 
;   refgrid = rtoffset[idxbest,*]
;   diffAC_interpol = dblarr(nscans, 6, ndata) & diffBD_interpol = diffAC_interpol
; 
;   for i=0,nscans-1 do begin
; 
;     for j=0,5 do begin
; 
;       diffAC_interpol[i,j,*] = interpol(reform(diffAC[i,j,*]), reform(rtoffset[i,*]), reform(refgrid), /spline)
;       diffBD_interpol[i,j,*] = interpol(reform(diffBD[i,j,*]), reform(rtoffset[i,*]), reform(refgrid), /spline)
; 
;     endfor
; 
;     proceeding_text,loop=(nscans), i=i, prompt='> Interpolation of Fringe Package'+string(i+1,form='(I4)')
; 
;   endfor
; 
;   ;-------------------------------------------------------------------------------------------------
; 
;   ;compute mean fringe package
; 
;   diffACave = dblarr(6,ndata) & diffBDave = diffACave
;   for i=0,5 do begin
; 
;     for j=0,ndata-1 do begin
; 
;     diffACave[i,j] = mean(diffAC_interpol[*,i,j], /nan)
;     diffBDave[i,j] = mean(diffBD_interpol[*,i,j], /nan)
; 
;     endfor
; 
;   endfor
; 
;   ;AVERAGING DOES NOT WORK, WASHES OUT THE FRINGE
; 


  ;=================================================================================================

  ;visual check if fringes are superimposed and how they look like

  if (pxmode eq '3') then begin

    print, '> Plotting Normalized Fringes'

    multiplot,/default

    set_plot, 'ps'
    device, filename=resultpath+fileid+'_'+starflag[xx]+'_'+target+'_NormalizedFringeScan.ps',/color,XSIZE=30, YSIZE=20, XOffset=xoffset, YOffset=yoffset
    !p.thick=4
    !x.thick=3
    !y.thick=3

    multiplot, [4,2], mxtitle='OPD*10!U-6!N / [m]', mytitle='Normalized Flux', mTitOffset=-0.5, mxTitSize=1.5, mxTitOffset=1, myTitSize=1.5, myTitOffset=2, gap=0.0

    pixel = ['0', '2', '3', '4']

      for j=0,3 do begin

	if (j eq 0) then begin

	  plot, rtoffset[0,*]*1.d6, diffAC[0,j,*]/factorAC[0], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffAC[i,j,*]/factorAC[i]
	  if (j eq 0) then legend, ['AC'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5
; 	  xyouts, 0.13, 0.58, 'Flux / Scale Factor', /normal, charsize=1.5

	endif else begin

	  plot, rtoffset[0,*]*1.d6, diffAC[0,j,*], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffAC[i,j,*]
	  if (j eq 0) then legend, ['AC'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5

	endelse

      multiplot

      endfor

      for j=0,3 do begin

	if (j eq 0) then begin

	  plot, rtoffset[0,*]*1.d6, diffBD[0,j,*]/factorBD[0], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffBD[i,j,*]/factorBD[i]
	  if (j eq 0) then legend, ['BD'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5
; 	  xyouts, 0.13, 0.16, 'Flux / Scale Factor', /normal, charsize=1.5

	endif else begin

	  plot, rtoffset[0,*]*1.d6, diffBD[0,j,*], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffBD[i,j,*]
	  if (j eq 0) then legend, ['BD'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5

	endelse

      multiplot

      endfor

      multiplot, /default

    device,/close
    set_plot,'x'
    !p.multi=[0,1,0]
    !p.thick=1
    !x.thick=1
    !y.thick=1

  endif

  ;-----------------------------------------

  if (pxmode eq '5') then begin

    print, '> Plotting Normalized Fringes'

    multiplot,/default

    set_plot, 'ps'
    device, filename=resultpath+fileid+'_'+starflag[xx]+'_'+target+'_NormalizedFringeScan.ps',/color,XSIZE=40, YSIZE=20, XOffset=xoffset, YOffset=yoffset
    !p.thick=4
    !x.thick=3
    !y.thick=3

    multiplot, [6,2], mxtitle='OPD*10!U-6!N / [m]', mytitle='Normalized Flux', mTitOffset=-0.5, mxTitSize=1.5, mxTitOffset=1, myTitSize=1.5, myTitOffset=2, gap=0.0

    pixel = ['0', '1', '2', '3', '4', '5']

      for j=0,5 do begin

	if (j eq 0) then begin

	  plot, rtoffset[0,*]*1.d6, diffAC[0,j,*]/factorAC[0], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata;, xr=[-100,100]
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffAC[i,j,*]/factorAC[i]
	  if (j eq 0) then legend, ['AC'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5
; 	  xyouts, 0.09, 0.58, 'Flux / Scale Factor', /normal, charsize=1.5

	endif else begin

	  plot, rtoffset[0,*]*1.d6, diffAC[0,j,*], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata;, xr=[-100,100]
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffAC[i,j,*]
	  if (j eq 0) then legend, ['AC'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5

	endelse

      multiplot

      endfor

      for j=0,5 do begin

	if (j eq 0) then begin

	  plot, rtoffset[0,*]*1.d6, diffBD[0,j,*]/factorBD[0], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata;, xr=[-100,100]
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffBD[i,j,*]/factorBD[i]
	  if (j eq 0) then legend, ['BD'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5
; 	  xyouts, 0.09, 0.16, 'Flux / Scale Factor', /normal, charsize=1.5

	endif else begin

	  plot, rtoffset[0,*]*1.d6, diffBD[0,j,*], yr=[-2.5,2.5], yst=1, xst=1, charsize=1.2, /nodata;, xr=[-100,100]
	  for i=0,nscans-1 do oplot, rtoffset[i,*]*1.d6, diffBD[i,j,*]
	  if (j eq 0) then legend, ['BD'], box=0, margin=-0.5, /left, /top, charsize=1.5
	  legend, ['Pixel: '+pixel[j]], box=0, margin=-0.5, /left, /bottom, charsize=1.5

	endelse

      multiplot

      endfor

      multiplot, /default

    device,/close
    set_plot,'x'
    !p.multi=[0,1,0]
    !p.thick=1
    !x.thick=1
    !y.thick=1

    ; spawn, 'gv '+resultpath+fileid+'_'+starflag[xx]+'_'+target+'_NormalizedFringeScan.ps'

  endif

  ;=================================================================================================

  ;Histogram of coherence factor per pixel

  plotflag = '1'

  if (pxmode eq '3') then begin

    if (plotflag eq '1') then begin

      print, '> Writing Histogram'

      set_plot, 'ps'
      device, filename=resultpath+fileid+'_'+starflag[xx]+'_'+target+'_HistogramCoherenceFactor.ps',/color,XSIZE=35, YSIZE=25, XOffset=xoffset, YOffset=yoffset
      ;   !P.Font=0
	!p.thick=4
	!x.thick=3
	!y.thick=3
	!p.multi=[0,2,2]

	for i=0,3 do begin

	  bw_cohAC[i] = binsizeFSU(cohAC[*,i])
  ; 	bw_cohBD[i] = binsizeFSU(cohBD[*,i])
	  bw_cohBD[i] = bw_cohAC[i]

	  if (i eq 0) then title = 'White Light Pixel'
	  if (i eq 1) then title = '2nd Pixel'
	  if (i eq 2) then title = '3rd Pixel'
	  if (i eq 3) then title = '4th Pixel'

	  histoplot, cohAC[*,i], binsize=bw_cohAC[i], charsize=1.5, backcolorname=cgcolor('white'), $
	  axiscolorname=cgcolor('black'), datacolorname=cgcolor('black'), polycolor=cgcolor('blue'), /fillpolygon, $
	  xtitle='Coherence Factor / '+cgGreek('mu')+'!U2!N [%]', title=title, /nan, /frequency, ytickformat='(f4.2)';ytitle='Number of values', yr=[0,1], yst=1
	  plots, bootmean_cohAC[i]*[1.,1.], !y.crange, color=fsc_color('powder blue')

	  histoplot, cohBD[*,i], binsize=bw_cohBD[i], /line_fill, orientation=45, /oplot, polycolor=cgcolor('gray'), datacolorname=cgcolor('gray'), spacin=0.1, /nan, /frequency
	  plots, bootmean_cohBD[i]*[1.,1.], !y.crange, color=fsc_color('dark gray')

	  if (i eq 0) then legend, ['AC', 'BD'], /right, /top, box=0, margin=0, psym=[sym(5),sym(5)], color=[fsc_color('blue'), fsc_color('gray')], symsize=[1.0,1.0]
	endfor

      device,/close
      set_plot,'x'

      !p.multi=[0,1,0]
      ; !P.Font=0
      !p.thick=1
      !x.thick=1
      !y.thick=1

    ;   spawn, 'gv '+resultpath+target+'_HistrogramCoherenceFactor.ps'

    endif

  endif

  ;---------------------------------------------------------

  if (pxmode eq '5') then begin

    if (plotflag eq '1') then begin

      print, '> Writing Histogram'

      ; bw_cohAC = [0.002, 0.004, 0.006, 0.008, 0.01, 0.012]
      ; bw_cohAC = [0.005, 0.005, 0.008, 0.008, 0.01, 0.01]
      ; bw_cohAC = [0.0125, 0.012, 0.011, 0.012, 0.0388, 0.021]	;computed for 333_01
  ;     bw_coh = [0.005,0.0125,0.0125,0.0125,0.0125,0.05]

      set_plot, 'ps'
      device, filename=resultpath+fileid+'_'+starflag[xx]+'_'+target+'_HistrogramCoherenceFactor.ps',/color,XSIZE=40, YSIZE=20, XOffset=xoffset, YOffset=yoffset
      ;   !P.Font=0
	!p.thick=4
	!x.thick=3
	!y.thick=3
	!p.multi=[0,3,2]

	for i=0,5 do begin

	  bw_cohAC[i] = binsizeFSU(cohAC[*,i])
  ; 	bw_cohBD[i] = binsizeFSU(cohBD[*,i])
	  bw_cohBD[i] = bw_cohAC[i]

	  if (i eq 0) then title = 'White Light Pixel'
	  if (i eq 1) then title = '1st Pixel'
	  if (i eq 2) then title = '2nd Pixel'
	  if (i eq 3) then title = '3rd Pixel'
	  if (i eq 4) then title = '4th Pixel'
	  if (i eq 5) then title = '5th Pixel'

	  histoplot, cohAC[*,i], binsize=bw_cohAC[i], charsize=2, backcolorname=cgcolor('white'), $
	  axiscolorname=cgcolor('black'), datacolorname=cgcolor('black'), polycolor=cgcolor('blue'), /fillpolygon, $
	  xtitle='Coherence Factor / '+cgGreek('mu')+'!U2!N [%]', title=title, /nan, /frequency, ytickformat='(f4.2)';ytitle='Number of values', yr=[0,1], yst=1
	  plots, bootmean_cohAC[i]*[1.,1.], !y.crange, color=fsc_color('powder blue')

	  histoplot, cohBD[*,i], binsize=bw_cohBD[i], /line_fill, orientation=45, /oplot, polycolor=cgcolor('gray'), datacolorname=cgcolor('gray'), spacin=0.1, /nan, /frequency
	  plots, bootmean_cohBD[i]*[1.,1.], !y.crange, color=fsc_color('dark gray')

	  if (i eq 0) then legend, ['AC', 'BD'], /right, /top, box=0, margin=0, psym=[sym(5),sym(5)], color=[fsc_color('blue'), fsc_color('gray')], symsize=[1.0,1.0]
	endfor

      device,/close
      set_plot,'x'

      !p.multi=[0,1,0]
      ; !P.Font=0
      !p.thick=1
      !x.thick=1
      !y.thick=1

    ;   spawn, 'gv '+resultpath+target+'_HistogramCoherenceFactor.ps'

    endif

  endif


  ;=================================================================================================

  ;theoretical visibility of calibrator

  ;projected baseline
  config = get_baseline(path+obsid[xx]+'/', file)

  ;-------------------------------------------------------------------------------------------------

  if (starflag[xx] eq 'Cal') then begin

    theovis2 = dblarr(npix)	;squared visibility for AC

    for i=0,npix-1 do theovis2[i] = (get_visibility(config.bl, ave_wl_cal[i], diam[xx]/1000.d0, 0.))^2.

    tmp = int_tabulated(ave_wl_cal[1:*],theovis2[1:*]*ave_wl_cal[1:*]^2.,/double)/int_tabulated(ave_wl_cal[1:*], ave_wl_cal[1:*]^2.)	;theoretical Visibility due to bandpass smearing, Aufdenberg 2006
    theovis2[0] = tmp

    T2AC = bootmean_cohAC/theovis2	;squared transfer function
    T2BD = bootmean_cohBD/theovis2

    T2AC_err = T2AC*bootsd_cohAC/bootmean_cohAC
    T2BD_err = T2BD*bootsd_cohBD/bootmean_cohBD

;     print, 'Theo. Vis^2 ', theovis2

  endif

  ;V^2 = V_cal^2 * coh^2 / coh_cal^2
  ;V^2 = coh^2 / T^2	->  T^2 = coh_cal^2/V_cal^2


  ;=================================================================================================

  print, '> Writing Output File'

  openw, lun, resultpath+fileid+'_'+starflag[xx]+'_'+target+'_results.txt', width=1400, /get_lun

    printf, lun, target
    if (starflag[xx] eq 'Sci') then printf, lun, 'Science' else printf, lun, 'Calibrator'
    printf, lun, ''
    printf, lun, 'JD: ', mjd+0.5d0+2400000.d0, format='(a4, f16.8)'
    printf, lun, ''
    printf, lun, 'Files'
    printf, lun, '-----'
    printf, lun, 'Fringe scan file     : ', file
    printf, lun, 'Background file used : ', darkfile
    printf, lun, 'Flat filed files used: ', flatfile1, ' / ', flatfile2
    printf, lun, ''
    printf, lun, 'Scan Setup'
    printf, lun, '----------'
    printf, lun, 'Total number of scans : ', strcompress(round(nslew), /rem)

    if (ratio gt 50.) then printf, lun, 'Used scans            : ', strcompress(round(nscans), /rem), ' ','!!! More than 50% of fringe packages are too close to the edge !!!' else $
    printf, lun, 'Used scans            : ', strcompress(round(nscans), /rem)

    printf, lun, 'Sampling rate         : ', strcompress(round(samprate), /rem)
    printf, lun, 'Stroke / [micrometer] : ', strcompress(stroke, /rem)
    printf, lun, 'FSU-A frequency / [Hz]: ', strcompress(round(fsudit), /rem)
    printf, lun, ''
    printf, lun, 'Observational Conditions'
    printf, lun, '------------------------'
    printf, lun, 'Seeing / [arcsec]    : ', sigfig(seeing, 4)
    printf, lun, 'Coherence time / [ms]: ', sigfig(tau0, 4)
    printf, lun, 'Airmass              : ', sigfig(airmass, 4)
    printf, lun, 'Altitude / [deg]     : ', sigfig(altitude, 5)
    printf, lun, ''
    printf, lun, 'Telescope Setup'
    printf, lun, '---------------'
    printf, lun, 'Station: ', station
    printf, lun, 'uvw / [m]               : ', u, v, w, format='(a26,3f12.5)'
    printf, lun, 'Projected baseline / [m]: ', sigfig(bl, 6)
    printf, lun, 'Position angle / [deg]  : ', sigfig(padeg, 6)
    printf, lun, ''
    printf, lun, ''
    if (pxmode eq '3') then printf, lun, 'Bootstrap Mean coherence factor W-2-3-4 for AC'
    if (pxmode eq '3') then printf, lun, '----------------------------------------------'
    if (pxmode eq '5') then printf, lun, 'Bootstrap Mean coherence factor W-1-2-3-4-5 for AC'
    if (pxmode eq '5') then printf, lun, '--------------------------------------------------'
    printf, lun, bootmean_cohAC
    if (pxmode eq '3') then printf, lun, 'Bootstrap Mean coherence factor W-2-3-4 for BD'
    if (pxmode eq '3') then printf, lun, '----------------------------------------------'
    if (pxmode eq '5') then printf, lun, 'Bootstrap Mean coherence factor W-1-2-3-4-5 for BD'
    if (pxmode eq '5') then printf, lun, '---------------------------------------------------'
    printf, lun, bootmean_cohBD
    printf, lun, 'Bootstrap Standard Deviation AC, BD'
    printf, lun, '-----------------------------------'
    printf, lun, bootsd_cohAC, bootsd_cohBD
    printf, lun, 'Skewness AC, BD'
    printf, lun, '---------------'
    printf, lun, skewness_cohAC, skewness_cohBD
    printf, lun, 'Kurtosis AC, BD'
    printf, lun, '---------------'
    printf, lun, kurtosis_cohAC, kurtosis_cohBD
    printf, lun, ''
    printf, lun, ''

    if (starflag[xx] eq 'Cal') then begin

      printf, lun, 'V^2_theo    : '
      printf, lun, theovis2
      printf, lun, 'T^2 AC     : '
      printf, lun, T2AC
      printf, lun, 'T^2 AC ERR : '
      printf, lun, T2AC_err
      printf, lun, 'T^2 BD     : '
      printf, lun, T2BD
      printf, lun, 'T^2 BD ERR  : '
      printf, lun, T2BD_err
      printf, lun, ''

    endif

    n = indgen(nscans)+1

    if (pxmode eq '3') then begin

      printf, lun, '#Scan                                                  Coh.Fac.^2 AC                                                                                              Coh.Fac.^2 BD'
      printf, lun, '              White AC            2nd AC            3rd AC            4th AC          White BD            2nd BD            3rd BD            4th BD'
      for i=0,nscans-1 do begin

	  printf, lun, n[i], cohAC[i,*], cohBD[i,*], format='(i4, 8f18.10)'

      endfor

    endif

    if (pxmode eq '5') then begin

      printf, lun, '#Scan                                                  Coh.Fac.^2 AC                                                                                              Coh.Fac.^2 BD'
      printf, lun, '              White AC            1st AC            2nd AC            3rd AC            4th AC            5th AC          White BD            1st BD            2nd BD            3rd BD            4th BD            5th BD'
      for i=0,nscans-1 do begin

	  printf, lun, n[i], cohAC[i,*], cohBD[i,*], format='(i4, 12f18.10)'

      endfor

    endif


  close, lun
  free_lun, lun


  print, '> Done.'
  print, ''

endfor


stop
end
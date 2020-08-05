@closest.pro
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
@ts_diff.pro
@mpfitfun.pro
@mpfit.pro
@reverse.pro
@interpol.pro
@mean.pro
@moment.pro
@rgb.pro
@linspace.pro
@c_correlate.pro


;===============================================================================================

;no improvement if used instead auf Gauss function
; function fit_snr, x, p
;   fit = p[7]*( exp(-0.5d0*((x-p[0])/p[1])^2.) + p[4]*exp(-0.5d0*((x-p[5])/p[6])^2.)*sin(p[2]*x-p[3]) )+p[8]
;   return, fit
; end

function fit_fringe, x, p

  ;p0: I0, p1: l0, p2: dl0, p3: doffset, p4: Ioffset, p5: Phase

  lcoh = (p[1]^2.-(p[2]^2.)/4.d0)/p[2]
  eta0 = 1.d0

  k0 = 2.d0*!DPI/p[1]
  sinc = (sin(!DPI*(x-p[3])/lcoh))/(!DPI*(x-p[3])/lcoh)
  idxnan = where(finite(sinc) ne 1)	;in case there is a NAN when fringe going through 0 OPD
  if (idxnan[0] ne -1) then sinc[idxnan] = 1.d0
;   Iint = (2.d0*p[0]*p[2]*eta0*sinc*sin(k0*(x-p[3])))+p[4]
  Iint = (2.d0*p[0]*p[2]*eta0*sinc*sin(k0*(x-p[3])+p[5]))+p[4]

  return, Iint

end

function create_fringe, x, p

  ;p0: I0, p1: l0, p2: dl0, p3: doffset, p4: Ioffset, p5: Phase

  lcoh = (p[1]^2.-(p[2]^2.)/4.d0)/p[2]
  eta0 = 1.d0

  k0 = 2.d0*!DPI/p[1]
  sinc = (sin(!DPI*(x-p[3])/lcoh))/(!DPI*(x-p[3])/lcoh)
  idxnan = where(finite(sinc) ne 1)	;in case there is a NAN when fringe going through 0 OPD
  if (idxnan[0] ne -1) then sinc[idxnan] = 1.d0
  ;Iint = (2.d0*p[0]*p[2]*eta0*sinc*sin(k0*(x-p[3])))+p[4]
  Iint = (2.d0*p[0]*p[2]*eta0*sinc*sin(k0*(x-p[3])+p[5]))+p[4]

  return, Iint

end

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

function normalize, x

  return, (x-min(x))/(max(x)-min(x))

end

function snrfun, x, p
  fit = p[7]*( exp(-0.5d0*((x-p[0])/p[1])^2.) + p[4]*exp(-0.5d0*((x-p[5])/p[6])^2.)*sin(p[2]*x-p[3]) )+p[8]
  return, fit
end

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


pro calibrate_FSUA_SKY_SIMULATION


; ;define path for results
;   resultpath = path+'SkyCal_'+fileid+'/'
;   spawn, 'mkdir '+resultpath

;parameters
;sampling = [44, 110, 122, 134, 146, 158]	;sampling similar to whats used for real data
l0 = [2.25d-6, 2.0d-6, 2.1d-6, 2.2d-6, 2.3d-6, 2.4d-6]	;lambda
dl0 = [4.5d-7, 1.4d-7, 1.4d-7, 1.4d-7, 1.4d-7, 1.4d-7 ]	;delta lambda, bandwidth
p0 = ([0., 90., 180., 270.])*!dtor	;phase shifts between ABCD
base = -3.d-3	;in meter, ZOPD Offset
step = 5.d-7	;according to 2 samples per micrometer
;cut = [1.11d-05, 2.755d-05, 3.038d-05, 3.334d-05, 3.645d-05, 3.969d-05]	;+/- micrometer
cut =  [11.d-6, 27.d-6, 20.d-6, 33.d-6, 35.d-6, 35.d-6]


;white light pixe should be not so narrow cutted, has influence on wavelength fit

Ioffset = 0.5d0
I0 = 1000000.d0

; lcoh = (l0^2.-(dl0^2.)/4.d0)/dl0
; k0 = 2.d0*!DPI/l0
; sinc = (sin(!DPI*(delta-base)/lcoh))/(!DPI*(delta-base)/lcoh)
; Iint = (2.d0*I0*dl0*eta0*sinc*sin(k0*(delta-base)+p0))+Ioffset
; idx = where(finite(Iint) ne 1)
; Iint[idx] = Ioffset


;reference grid
refgridW = dindgen(2.*cut[0]/step)*step
refgrid1 = dindgen(2.*cut[1]/step)*step
refgrid2 = dindgen(2.*cut[2]/step)*step
refgrid3 = dindgen(2.*cut[3]/step)*step
refgrid4 = dindgen(2.*cut[4]/step)*step
refgrid5 = dindgen(2.*cut[5]/step)*step

refgridW = refgridW-max(refgridW)/2.d0+base
refgrid1 = refgrid1-max(refgrid1)/2.d0+base
refgrid2 = refgrid2-max(refgrid2)/2.d0+base
refgrid3 = refgrid3-max(refgrid3)/2.d0+base
refgrid4 = refgrid4-max(refgrid4)/2.d0+base
refgrid5 = refgrid5-max(refgrid5)/2.d0+base


;create fringe

avenormW = dblarr(4, n_elements(refgridW))
avenorm1 = dblarr(4, n_elements(refgrid1))
avenorm2 = dblarr(4, n_elements(refgrid2))
avenorm3 = dblarr(4, n_elements(refgrid3))
avenorm4 = dblarr(4, n_elements(refgrid4))
avenorm5 = dblarr(4, n_elements(refgrid5))

for i=0,3 do begin

  startval = [I0, l0[0], dl0[0], base, Ioffset, p0[i]]
  avenormW[i,*] = normalize(create_fringe(refgridW, startval))
  startval = [I0, l0[1], dl0[1], base, Ioffset, p0[i]]
  avenorm1[i,*] = normalize(create_fringe(refgrid1, startval))
  startval = [I0, l0[2], dl0[2], base, Ioffset, p0[i]]
  avenorm2[i,*] = normalize(create_fringe(refgrid2, startval))
  startval = [I0, l0[3], dl0[3], base, Ioffset, p0[i]]
  avenorm3[i,*] = normalize(create_fringe(refgrid3, startval))
  startval = [I0, l0[4], dl0[4], base, Ioffset, p0[i]]
  avenorm4[i,*] = normalize(create_fringe(refgrid4, startval))
  startval = [I0, l0[5], dl0[5], base, Ioffset, p0[i]]
  avenorm5[i,*] = normalize(create_fringe(refgrid5, startval))

endfor

;save original values without noise
foW = avenormW
fo1 = avenorm1
fo2 = avenorm2
fo3 = avenorm3
fo4 = avenorm4
fo5 = avenorm5
loW = refgridW
lo1 = refgrid1
lo2 = refgrid2
lo3 = refgrid3
lo4 = refgrid4
lo5 = refgrid5

;add noise to the fringe package
nf = 0.2
;nf = 0.0
for i=0,3 do begin

  ;(s/abs(s)) is only to get the sign

  z = randomu(seed, n_elements(avenormW[0,*]), /double)
  s = randomn(seed, n_elements(avenormW[0,*]), /double)
  avenormW[i,*] = avenormW[i,*]+avenormW[i,*]*(s/abs(s))*z*nf

  z = randomu(seed, n_elements(avenorm1[0,*]), /double)
  s = randomn(seed, n_elements(avenorm1[0,*]), /double)
  avenorm1[i,*] = avenorm1[i,*]+avenorm1[i,*]*(s/abs(s))*z*nf

  z = randomu(seed, n_elements(avenorm2[0,*]), /double)
  s = randomn(seed, n_elements(avenorm2[0,*]), /double)
  avenorm2[i,*] = avenorm2[i,*]+avenorm2[i,*]*(s/abs(s))*z*nf

  z = randomu(seed, n_elements(avenorm3[0,*]), /double)
  s = randomn(seed, n_elements(avenorm3[0,*]), /double)
  avenorm3[i,*] = avenorm3[i,*]+avenorm3[i,*]*(s/abs(s))*z*nf

  z = randomu(seed, n_elements(avenorm4[0,*]), /double)
  s = randomn(seed, n_elements(avenorm4[0,*]), /double)
  avenorm4[i,*] = avenorm4[i,*]+avenorm4[i,*]*(s/abs(s))*z*nf

  z = randomu(seed, n_elements(avenorm5[0,*]), /double)
  s = randomn(seed, n_elements(avenorm5[0,*]), /double)
  avenorm5[i,*] = avenorm5[i,*]+avenorm5[i,*]*(s/abs(s))*z*nf

endfor


;add noise to the wavelength grid
nf = 0.0001
;nf = 0.0
z = randomu(seed, n_elements(refgridW), /double)
s = randomn(seed, n_elements(refgridW), /double)
refgridW = refgridW+refgridW*(s/abs(s))*z*nf
z = randomu(seed, n_elements(refgrid1), /double)
s = randomn(seed, n_elements(refgrid1), /double)
refgrid1 = refgrid1+refgrid1*(s/abs(s))*z*nf
z = randomu(seed, n_elements(refgrid2), /double)
s = randomn(seed, n_elements(refgrid2), /double)
refgrid2 = refgrid2+refgrid2*(s/abs(s))*z*nf
z = randomu(seed, n_elements(refgrid3), /double)
s = randomn(seed, n_elements(refgrid3), /double)
refgrid3 = refgrid3+refgrid3*(s/abs(s))*z*nf
z = randomu(seed, n_elements(refgrid4), /double)
s = randomn(seed, n_elements(refgrid4), /double)
refgrid4 = refgrid4+refgrid4*(s/abs(s))*z*nf
z = randomu(seed, n_elements(refgrid5), /double)
s = randomn(seed, n_elements(refgrid5), /double)
refgrid5 = refgrid5+refgrid5*(s/abs(s))*z*nf

;normalize the fringe package
for i=0,3 do begin

  avenormW[i,*] = normalize(avenormW[i,*])
  avenorm1[i,*] = normalize(avenorm1[i,*])
  avenorm2[i,*] = normalize(avenorm2[i,*])
  avenorm3[i,*] = normalize(avenorm3[i,*])
  avenorm4[i,*] = normalize(avenorm4[i,*])
  avenorm5[i,*] = normalize(avenorm5[i,*])

endfor


;plot fringes with and w/o noise for all 24 pixels
  window, 0, xs=1800, ys=1000
  !p.multi = [0,4,6]
  for i=0,3 do begin
     plot, (refgridW-base)*1.d6, avenormW[i,*], xst=1, charsize=2, xr=[-30,30]
    oplot, (loW-base)*1.d6, foW[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid1-base)*1.d6, avenorm1[i,*], xst=1, charsize=2, xr=[-30,30]
    oplot, (lo1-base)*1.d6, fo1[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid2-base)*1.d6, avenorm2[i,*], xst=1, charsize=2, xr=[-30,30]
    oplot, (lo2-base)*1.d6, fo2[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid3-base)*1.d6, avenorm3[i,*], xst=1, charsize=2, xr=[-30,30]
    oplot, (lo3-base)*1.d6, fo3[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid4-base)*1.d6, avenorm4[i,*], xst=1, charsize=2, xr=[-30,30]
    oplot, (lo4-base)*1.d6, fo4[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid5-base)*1.d6, avenorm5[i,*], xst=1, charsize=2, xr=[-30,30]
    oplot, (lo5-base)*1.d6, fo5[i,*], color=rgb(255,0,0)
  endfor
  !p.multi = [0,1,0]

;===========================================================================================

;fit fringe packacge

fitparams = dblarr(4,6,6)	;channel-pixel-params
fringefitW = dblarr(4,n_elements(refgridW))
fringefit1 = dblarr(4,n_elements(refgrid1))
fringefit2 = dblarr(4,n_elements(refgrid2))
fringefit3 = dblarr(4,n_elements(refgrid3))
fringefit4 = dblarr(4,n_elements(refgrid4))
fringefit5 = dblarr(4,n_elements(refgrid5))

  ;p0: Intensity, p1: wavelength, p2: bandwidth, p3: ZOPD offset, p4: Ioffset, p5: Phase
startvalW = [500000.d0, 2.25d-6,5.d-7, base, 0.5d0 , 0.d0]
startval1 = [500000.d0, 2.0d-6, 1.4d-7, base, 0.5d0, 0.d0]
startval2 = [500000.d0, 2.1d-6, 1.4d-7, base, 0.5d0, 0.d0]
startval3 = [500000.d0, 2.2d-6, 1.4d-7, base, 0.5d0, 0.d0]
startval4 = [500000.d0, 2.3d-6, 1.4d-7, base, 0.5d0, 0.d0]
startval5 = [500000.d0, 2.4d-6, 1.4d-7, base, 0.5d0, 0.d0]
;phasestart = ([0.d0, 90.d0, 180.d0, 270.d0])*!dtor
phasestart = ([0.d0, 90.d0, 180.d0, 270.d0])*!dtor
;phasestart = ([0.d0, 0.d0, 0.d0, 0.d0])*!dtor

;weights, computed from sinc function
lcoh = (startvalW[1]^2.-(startvalW[2]^2.)/4.d0)/startvalW[2]
weightW = (sin(!DPI*(refgridW-startvalW[3])/lcoh))/(!DPI*(refgridW-startvalW[3])/lcoh)
idx = where(weightW lt 0.)
if (idx[0] ne -1) then weightW[idx] = 0.d0

lcoh = (startval1[1]^2.-(startval1[2]^2.)/4.d0)/startval1[2]
weight1 = (sin(!DPI*(refgrid1-startval1[3])/lcoh))/(!DPI*(refgrid1-startval1[3])/lcoh)
idx = where(weight1 lt 0.)
if (idx[0] ne -1) then weight1[idx] = 0.d0

lcoh = (startval2[1]^2.-(startval2[2]^2.)/4.d0)/startval2[2]
weight2 = (sin(!DPI*(refgrid2-startval2[3])/lcoh))/(!DPI*(refgrid2-startval2[3])/lcoh)
idx = where(weight2 lt 0.)
if (idx[0] ne -1) then weight2[idx] = 0.d0

lcoh = (startval3[1]^2.-(startval3[2]^2.)/4.d0)/startval3[2]
weight3 = (sin(!DPI*(refgrid3-startval3[3])/lcoh))/(!DPI*(refgrid3-startval3[3])/lcoh)
idx = where(weight3 lt 0.)
if (idx[0] ne -1) then weight3[idx] = 0.d0

lcoh = (startval4[1]^2.-(startval4[2]^2.)/4.d0)/startval4[2]
weight4 = (sin(!DPI*(refgrid4-startval4[3])/lcoh))/(!DPI*(refgrid4-startval4[3])/lcoh)
idx = where(weight4 lt 0.)
if (idx[0] ne -1) then weight4[idx] = 0.d0

lcoh = (startval5[1]^2.-(startval5[2]^2.)/4.d0)/startval5[2]
weight5 = (sin(!DPI*(refgrid5-startval5[3])/lcoh))/(!DPI*(refgrid5-startval5[3])/lcoh)
idx = where(weight5 lt 0.)
if (idx[0] ne -1) then weight5[idx] = 0.d0

;if fringe fit should be performed without wheight, uncomment this block
; weightW[*] = 1.d0
; weight1[*] = 1.d0
; weight2[*] = 1.d0
; weight3[*] = 1.d0
; weight4[*] = 1.d0
; weight5[*] = 1.d0


;pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},5)
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
pi[1].limits[0] = 1.8d-6	;wavelength limited between 1.8 and 2.7 micrometer
pi[1].limits[1] = 2.7d-6
pi[1].limited[0] = 1
pi[1].limited[1] = 1
pi[2].limits[0] = 0.5d-7	;bandbass limited between 0.5 and 1 micrometer
pi[2].limits[1] = 10.d-7
pi[2].limited[0] = 1
pi[2].limited[1] = 1
pi[4].limits[0] = 0.d0		;Ioffset limited between 0 and 1
pi[4].limits[1] = 1.d0
pi[4].limited[0] = 1
pi[4].limited[1] = 1
pi[5].limited[0] = 1
pi[5].limited[1] = 1

for i=0,3 do begin

  startvalW[5] = phasestart[i]
  startval1[5] = phasestart[i]
  startval2[5] = phasestart[i]
  startval3[5] = phasestart[i]
  startval4[5] = phasestart[i]
  startval5[5] = phasestart[i]

  if (i eq 0) then begin
    pi[5].limits[0] = -100.*!dtor
    pi[5].limits[1] = 100.*!dtor
  endif
  if (i eq 1) then begin
    pi[5].limits[0] = 20.*!dtor
    pi[5].limits[1] = 120.*!dtor
  endif
  if (i eq 2) then begin
    pi[5].limits[0] = 120.*!dtor
    pi[5].limits[1] = 250.*!dtor
  endif
  if (i eq 3) then begin
    pi[5].limits[0] = 220.*!dtor
    pi[5].limits[1] = 340.*!dtor
  endif

  dummy_err = dblarr(n_elements(refgridW))
  dummy_err[*] = 1.d0
  paramW = mpfitfun('fit_fringe', refgridW, avenormW[i,*], dummy_err, startvalW, weights=weightW, maxiter=2000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, yfit=yfitW, perror=perror, dof=dof, /quiet, /nan)
  fitparams[i,0,*] = paramW
  fringefitW[i,*] = yfitW

  dummy_err = dblarr(n_elements(refgrid1))
  dummy_err[*] = 1.d0
  param1 = mpfitfun('fit_fringe', refgrid1, avenorm1[i,*], dummy_err, startval1, weights=weight1, maxiter=2000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit1, perror=perror, dof=dof, /quiet, /nan)
  fitparams[i,1,*] = param1
  fringefit1[i,*] = yfit1

  dummy_err = dblarr(n_elements(refgrid2))
  dummy_err[*] = 1.d0
  param2 = mpfitfun('fit_fringe', refgrid2, avenorm2[i,*], dummy_err, startval2, weights=weight2, maxiter=2000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit2, perror=perror, dof=dof, /quiet, /nan)
  fitparams[i,2,*] = param2
  fringefit2[i,*] = yfit2

  dummy_err = dblarr(n_elements(refgrid3))
  dummy_err[*] = 1.d0
  param3 = mpfitfun('fit_fringe', refgrid3, avenorm3[i,*], dummy_err, startval3, weights=weight3, maxiter=2000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit3, perror=perror, dof=dof, /quiet, /nan)
  fitparams[i,3,*] = param3
  fringefit3[i,*] = yfit3

  dummy_err = dblarr(n_elements(refgrid4))
  dummy_err[*] = 1.d0
  param4 = mpfitfun('fit_fringe', refgrid4, avenorm4[i,*], dummy_err, startval4, weights=weight4, maxiter=2000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit4, perror=perror, dof=dof, /quiet, /nan)
  fitparams[i,4,*] = param4
  fringefit4[i,*] = yfit4

  dummy_err = dblarr(n_elements(refgrid5))
  dummy_err[*] = 1.d0
  param5 = mpfitfun('fit_fringe', refgrid5, avenorm5[i,*], dummy_err, startval5, weights=weight5, maxiter=2000, parinfo=pi, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit5, perror=perror, dof=dof, /quiet, /nan)
  fitparams[i,5,*] = param5
  fringefit5[i,*] = yfit5


endfor

;plot fringe fits for all 24 pixels
;this plot should be always be plotted, to see if fit was successful
  window, 1, xs=1800, ys=1000
  !p.multi = [0,4,6]
  for i=0,3 do begin
    plot, (refgridW-base)*1.d6, avenormW[i,*], xst=1, charsize=2
    oplot, (refgridW-base)*1.d6, fringefitW[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid1-base)*1.d6, avenorm1[i,*], xst=1, charsize=2
    oplot, (refgrid1-base)*1.d6, fringefit1[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid2-base)*1.d6, avenorm2[i,*], xst=1, charsize=2
    oplot, (refgrid2-base)*1.d6, fringefit2[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid3-base)*1.d6, avenorm3[i,*], xst=1, charsize=2
    oplot, (refgrid3-base)*1.d6, fringefit3[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid4-base)*1.d6, avenorm4[i,*], xst=1, charsize=2
    oplot, (refgrid4-base)*1.d6, fringefit4[i,*], color=rgb(255,0,0)
  endfor
  for i=0,3 do begin
     plot, (refgrid5-base)*1.d6, avenorm5[i,*], xst=1, charsize=2
    oplot, (refgrid5-base)*1.d6, fringefit5[i,*], color=rgb(255,0,0)
  endfor
  !p.multi = [0,1,0]


;create a higher resolution refgrid for each pixel and recompute theoretical fringe package
;hr = high resolution
hrgridW = linspace(min(refgridw), max(refgridw), 200)
hrgrid1 = linspace(min(refgrid1), max(refgrid1), 1000)
hrgrid2 = linspace(min(refgrid2), max(refgrid2), 1100)
hrgrid3 = linspace(min(refgrid3), max(refgrid3), 1200)
hrgrid4 = linspace(min(refgrid4), max(refgrid4), 1300)
hrgrid5 = linspace(min(refgrid5), max(refgrid5), 1400)

hrfringeW = dblarr(4,n_elements(hrgridW))
hrfringe1 = dblarr(4,n_elements(hrgrid1))
hrfringe2 = dblarr(4,n_elements(hrgrid2))
hrfringe3 = dblarr(4,n_elements(hrgrid3))
hrfringe4 = dblarr(4,n_elements(hrgrid4))
hrfringe5 = dblarr(4,n_elements(hrgrid5))

for i=0,3 do begin

  hrfringeW[i,*] = fit_fringe(hrgridW, fitparams[i,0,*])
  hrfringe1[i,*] = fit_fringe(hrgrid1, fitparams[i,1,*])
  hrfringe2[i,*] = fit_fringe(hrgrid2, fitparams[i,2,*])
  hrfringe3[i,*] = fit_fringe(hrgrid3, fitparams[i,3,*])
  hrfringe4[i,*] = fit_fringe(hrgrid4, fitparams[i,4,*])
  hrfringe5[i,*] = fit_fringe(hrgrid5, fitparams[i,5,*])

endfor


;this is bullshit
; opddiff = dblarr(4,6)
; phshift = dblarr(4,6)
; for i=1,3 do begin
; 
;   opddiff[i,*] = fitparams[0,*,3]-fitparams[i,*,3]
;   phshift[i,*] = 2.d0*!DPI*opddiff[i,*]/fitparams[i,*,1]
;   phshift[i,*] = phshift[i,*] mod 2.d0*!DPI
;   phshift[i,*] = phshift[i,*]*(180.d0/!DPI)
; 
; endfor


fitwlen = fitparams[*,*,1]
;fitangval = phshift[1:3,*]
fitangval = fitparams[*,*,5]/!dtor

;under investigation
;some angles can be shifted by 180 degr. if this is the case its the same for all chanels
; for i=0,5 do begin
; 
;   if (fitangval[0,i] lt 0.) then fitangval[*,i] = fitangval[*,i]+180.d0
; 
; endfor

; print, ''
; print, 'Values from fitting the Fringe package'
; print, fitwlen
; print, fitangval



;============================================================================================
;the following part uses FFT and CC to determine wavelengths and phase shifts
;============================================================================================

;for FFT and CCF rename the arrays to original names from previous script, to save coding time
avenormW = normalize(hrfringeW)
avenorm1 = normalize(hrfringe1)
avenorm2 = normalize(hrfringe2)
avenorm3 = normalize(hrfringe3)
avenorm4 = normalize(hrfringe4)
avenorm5 = normalize(hrfringe5)


;FFT scans
fftW = dblarr(4,n_elements(avenormW[0,*]))
fft1 = dblarr(4,n_elements(avenorm1[0,*]))
fft2 = dblarr(4,n_elements(avenorm2[0,*]))
fft3 = dblarr(4,n_elements(avenorm3[0,*]))
fft4 = dblarr(4,n_elements(avenorm4[0,*]))
fft5 = dblarr(4,n_elements(avenorm5[0,*]))

for i=0,3 do begin

  fftW[i,*] = abs(fft(avenormW[i,*],-1))
  fft1[i,*] = abs(fft(avenorm1[i,*],-1))
  fft2[i,*] = abs(fft(avenorm2[i,*],-1))
  fft3[i,*] = abs(fft(avenorm3[i,*],-1))
  fft4[i,*] = abs(fft(avenorm4[i,*],-1))
  fft5[i,*] = abs(fft(avenorm5[i,*],-1))

endfor


;Number of the first sample that defines the tail of the fft.

;if you remove the first few elements this value has to be added later for the peak position
;if SNR allows set it to 0
; fftW = fftW[*,2:*]
; fft1 = fft1[*,2:*]
; fft2 = fft2[*,2:*]
; fft3 = fft3[*,2:*]
; fft4 = fft4[*,2:*]
; fft5 = fft5[*,2:*]
fftW[*,0:1] = 0.d0
fft1[*,0:1] = 0.d0
fft2[*,0:1] = 0.d0
fft3[*,0:1] = 0.d0
fft4[*,0:1] = 0.d0
fft5[*,0:1] = 0.d0


;only consider 1st half of FFT

fftW = fftW[*,0:round(n_elements(fftW[0,*])/2.)]
fft1 = fft1[*,0:round(n_elements(fft1[0,*])/2.)]
fft2 = fft2[*,0:round(n_elements(fft2[0,*])/2.)]
fft3 = fft3[*,0:round(n_elements(fft3[0,*])/2.)]
fft4 = fft4[*,0:round(n_elements(fft4[0,*])/2.)]
fft5 = fft5[*,0:round(n_elements(fft5[0,*])/2.)]


;GaussPeak

auxgridW = dindgen(n_elements(fftW[0,*]))
auxgrid1 = dindgen(n_elements(fft1[0,*]))
auxgrid2 = dindgen(n_elements(fft2[0,*]))
auxgrid3 = dindgen(n_elements(fft3[0,*]))
auxgrid4 = dindgen(n_elements(fft4[0,*]))
auxgrid5 = dindgen(n_elements(fft5[0,*]))

peakW = dblarr(4)
peak1 = dblarr(4)
peak2 = dblarr(4)
peak3 = dblarr(4)
peak4 = dblarr(4)
peak5 = dblarr(4)

;gauss_nterms = '3'
;     envelope, abs(xaa), 10., mnA, mxA, indA
;     envelope, abs(xaB), 10., mnB, mxB, indB
;     envelope, abs(xaC), 10., mnC, mxC, indC
;     envelope, abs(xaD), 10., mnD, mxD, indD
; 
;     startval = [1.d0, 0.d0, 500.]
;     resultA = gaussfit(lag[indA], mxA, paramA, estimates=startval, nterms=3)
;     resultB = gaussfit(lag[indB], mxB, paramB, estimates=startval, nterms=3)
;     resultC = gaussfit(lag[indC], mxC, paramC, estimates=startval, nterms=3)
;     resultD = gaussfit(lag[indD], mxD, paramD, estimates=startval, nterms=3)
dummy_err_W = dblarr(n_elements(auxgridW))
dummy_err_1 = dblarr(n_elements(auxgrid1))
dummy_err_2 = dblarr(n_elements(auxgrid2))
dummy_err_3 = dblarr(n_elements(auxgrid3))
dummy_err_4 = dblarr(n_elements(auxgrid4))
dummy_err_5 = dblarr(n_elements(auxgrid5))

startval = transpose([[0.1d0, 7.d0, 0.5d0], [0.1d0, 34.d0, 1.0d0], $
 	      [0.1d0, 37.d0, 1.0d0], [0.1d0, 40.d0, 1.0d0], $
 	      [0.1d0, 41.d0, 1.0d0], [0.1d0, 42.d0, 1.0d0]])

gauss_nterms = '3'

for i=0,3 do begin

  dummy_err = dblarr(n_elements(auxgridW))
  dummy_err[*] = 1.d0
  paramW = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', auxgridW, fftW[i,*], dummy_err, [max(fftW[i,*], maxpos), maxpos, 0.5], weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfitW, perror=perror, dof=dof, /quiet)

  dummy_err = dblarr(n_elements(auxgrid1))
  dummy_err[*] = 1.d0
  param1 = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', auxgrid1, fft1[i,*], dummy_err, [max(fft1[i,*], maxpos), maxpos, 0.5], weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit1, perror=perror, dof=dof, /quiet)

  dummy_err = dblarr(n_elements(auxgrid2))
  dummy_err[*] = 1.d0
  param2 = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', auxgrid2, fft2[i,*], dummy_err, [max(fft2[i,*], maxpos), maxpos, 0.5], weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit2, perror=perror, dof=dof, /quiet)

  dummy_err = dblarr(n_elements(auxgrid3))
  dummy_err[*] = 1.d0
  param3 = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', auxgrid3, fft3[i,*], dummy_err, [max(fft3[i,*], maxpos), maxpos, 0.5], weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit3, perror=perror, dof=dof, /quiet)

  dummy_err = dblarr(n_elements(auxgrid4))
  dummy_err[*] = 1.d0
  param4 = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', auxgrid4, fft4[i,*], dummy_err, [max(fft4[i,*], maxpos), maxpos, 0.5], weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit4, perror=perror, dof=dof, /quiet)

  dummy_err = dblarr(n_elements(auxgrid5))
  dummy_err[*] = 1.d0
  param5 = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', auxgrid5, fft5[i,*], dummy_err, [max(fft5[i,*], maxpos), maxpos, 0.5], weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit5, perror=perror, dof=dof, /quiet)

  peakW[i] = paramW[1]
  peak1[i] = param1[1]
  peak2[i] = param2[1]
  peak3[i] = param3[1]
  peak4[i] = param4[1]
  peak5[i] = param5[1]


;   resultW = gaussfit(auxgridW, fftW[i,*], paramW, estimates=[max(fftW[i,*], maxpos), maxpos, 0.5], nterms=3)
;   peakW[i] = paramW[1]
;   result1 = gaussfit(auxgrid1, fft1[i,*], param1, estimates=[max(fft1[i,*], maxpos), maxpos, 0.5], nterms=3)
;   peak1[i] = param1[1]
;   result2 = gaussfit(auxgrid2, fft2[i,*], param2, estimates=[max(fft2[i,*], maxpos), maxpos, 0.5], nterms=3)
;   peak2[i] = param2[1]
;   result3 = gaussfit(auxgrid3, fft3[i,*], param3, estimates=[max(fft3[i,*], maxpos), maxpos, 0.5], nterms=3)
;   peak3[i] = param3[1]
;   result4 = gaussfit(auxgrid4, fft4[i,*], param4, estimates=[max(fft4[i,*], maxpos), maxpos, 0.5], nterms=3)
;   peak4[i] = param4[1]
;   result5 = gaussfit(auxgrid5, fft5[i,*], param5, estimates=[max(fft5[i,*], maxpos), maxpos, 0.5], nterms=3)
;   peak5[i] = param5[1]


endfor


;stepsize = abs((max(refgrid)-min(refgrid))/double(n_elements(refgrid)))
stepsizeW = abs(median(ts_diff(hrgridW,1)))
stepsize1 = abs(median(ts_diff(hrgrid1,1)))
stepsize2 = abs(median(ts_diff(hrgrid2,1)))
stepsize3 = abs(median(ts_diff(hrgrid3,1)))
stepsize4 = abs(median(ts_diff(hrgrid4,1)))
stepsize5 = abs(median(ts_diff(hrgrid5,1)))
nsamples = double([n_elements(avenormW[0,*]),n_elements(avenorm1[0,*]),n_elements(avenorm2[0,*]), $
	  n_elements(avenorm3[0,*]),n_elements(avenorm4[0,*]),n_elements(avenorm5[0,*])])


;Peak2Wlen

dfreqW = (2.d0*!DPI/stepsizeW)/nsamples[0]
dfreq1 = (2.d0*!DPI/stepsize1)/nsamples[1]
dfreq2 = (2.d0*!DPI/stepsize2)/nsamples[2]
dfreq3 = (2.d0*!DPI/stepsize3)/nsamples[3]
dfreq4 = (2.d0*!DPI/stepsize4)/nsamples[4]
dfreq5 = (2.d0*!DPI/stepsize5)/nsamples[5]

wlenW = 2.d0*!DPI/(dfreqW*peakW-1.d0)
wlen1 = 2.d0*!DPI/(dfreq1*peak1-1.d0)
wlen2 = 2.d0*!DPI/(dfreq2*peak2-1.d0)
wlen3 = 2.d0*!DPI/(dfreq3*peak3-1.d0)
wlen4 = 2.d0*!DPI/(dfreq4*peak4-1.d0)
wlen5 = 2.d0*!DPI/(dfreq5*peak5-1.d0)

wlen = [[wlenw], [wlen1], [wlen2], [wlen3], [wlen4], [wlen5]]


;PhaseShifts

tmp = double(n_elements(avenormW[0,*]))
lagW = [dindgen(tmp-1.d0)-(tmp-1.d0), dindgen(tmp)]
tmp = double(n_elements(avenorm1[0,*]))
lag1 = [dindgen(tmp-1.d0)-(tmp-1.d0), dindgen(tmp)]
tmp = double(n_elements(avenorm2[0,*]))
lag2 = [dindgen(tmp-1.d0)-(tmp-1.d0), dindgen(tmp)]
tmp = double(n_elements(avenorm3[0,*]))
lag3 = [dindgen(tmp-1.d0)-(tmp-1.d0), dindgen(tmp)]
tmp = double(n_elements(avenorm4[0,*]))
lag4 = [dindgen(tmp-1.d0)-(tmp-1.d0), dindgen(tmp)]
tmp = double(n_elements(avenorm5[0,*]))
lag5 = [dindgen(tmp-1.d0)-(tmp-1.d0), dindgen(tmp)]


xaw = dblarr(4, n_elements(lagW))
xa1 = dblarr(4, n_elements(lag1))
xa2 = dblarr(4, n_elements(lag2))
xa3 = dblarr(4, n_elements(lag3))
xa4 = dblarr(4, n_elements(lag4))
xa5 = dblarr(4, n_elements(lag5))

for i=0,3 do begin

  xaW[i,*] = c_correlate(avenormW[0,*],avenormW[i,*], lagW, /double)
  xa1[i,*] = c_correlate(avenorm1[0,*],avenorm1[i,*], lag1, /double)
  xa2[i,*] = c_correlate(avenorm2[0,*],avenorm2[i,*], lag2, /double)
  xa3[i,*] = c_correlate(avenorm3[0,*],avenorm3[i,*], lag3, /double)
  xa4[i,*] = c_correlate(avenorm4[0,*],avenorm4[i,*], lag4, /double)
  xa5[i,*] = c_correlate(avenorm5[0,*],avenorm5[i,*], lag5, /double)

endfor

;have to make a seperate loop because I need all CCFs first

posW = dblarr(4)
pos1 = dblarr(4)
pos2 = dblarr(4)
pos3 = dblarr(4)
pos4 = dblarr(4)
pos5 = dblarr(4)
dabcdW = dblarr(4)
dabcd1 = dblarr(4)
dabcd2 = dblarr(4)
dabcd3 = dblarr(4)
dabcd4 = dblarr(4)
dabcd5 = dblarr(4)

lambdaW = mean(wlenW)
lambda1 = mean(wlen1)
lambda2 = mean(wlen2)
lambda3 = mean(wlen3)
lambda4 = mean(wlen4)
lambda5 = mean(wlen5)

for i=0,3 do begin

  if (i eq 0) then dum = max(normalize(xaW[i,*]), posWtmp) else $
    dum = max(normalize(xaW[i,0:posW[0]]), posWtmp)
  if (i eq 0) then $
    dum = max(normalize(xa1[i,*]), pos1tmp) else $
    dum = max(normalize(xa1[i,0:pos1[0]]), pos1tmp)
  if (i eq 0) then $
    dum = max(normalize(xa2[i,*]), pos2tmp) else $
    dum = max(normalize(xa2[i,0:pos2[0]]), pos2tmp)
  if (i eq 0) then $
    dum = max(normalize(xa3[i,*]), pos3tmp) else $
    dum = max(normalize(xa3[i,0:pos3[0]]), pos3tmp)
  if (i eq 0) then $
    dum = max(normalize(xa4[i,*]), pos4tmp) else $
    dum = max(normalize(xa4[i,0:pos4[0]]), pos4tmp)
  if (i eq 0) then $
    dum = max(normalize(xa5[i,*]), pos5tmp) else $
    dum = max(normalize(xa5[i,0:pos5[0]]), pos5tmp)

  posW[i] = posWtmp
  pos1[i] = pos1tmp
  pos2[i] = pos2tmp
  pos3[i] = pos3tmp
  pos4[i] = pos4tmp
  pos5[i] = pos5tmp

  dabcdW[i] = abs(lagW[posW[i]])*stepsizeW
  dabcd1[i] = abs(lag1[pos1[i]])*stepsize1
  dabcd2[i] = abs(lag2[pos2[i]])*stepsize2
  dabcd3[i] = abs(lag3[pos3[i]])*stepsize3
  dabcd4[i] = abs(lag4[pos4[i]])*stepsize4
  dabcd5[i] = abs(lag5[pos5[i]])*stepsize5

endfor


;   window, 0, xs=800, ys=1000
;   !p.multi=[0,1,6]
;   plot, lagW, normalize(xaW[0,*]), xst=1, charsize=2
;   oplot, lagW, normalize(xaW[1,*]), color=rgb(255,0,0)
;   oplot, lagW, normalize(xaW[2,*]), color=rgb(255,255,0)
;   oplot, lagW, normalize(xaW[3,*]), color=rgb(0,255,0)
;   plot, lag1, normalize(xa1[0,*]), xst=1, xr=[-180,180], charsize=2
;   oplot, lag1, normalize(xa1[1,*]), color=rgb(255,0,0)
;   oplot, lag1, normalize(xa1[2,*]), color=rgb(255,255,0)
;   oplot, lag1, normalize(xa1[3,*]), color=rgb(0,255,0)
;   plot, lag2, normalize(xa2[0,*]), xst=1, xr=[-180,180], charsize=2
;   oplot, lag2, normalize(xa2[1,*]), color=rgb(255,0,0)
;   oplot, lag2, normalize(xa2[2,*]), color=rgb(255,255,0)
;   oplot, lag2, normalize(xa2[3,*]), color=rgb(0,255,0)
;   plot, lag3, normalize(xa3[0,*]), xst=1, xr=[-180,180], charsize=2
;   oplot, lag3, normalize(xa3[1,*]), color=rgb(255,0,0)
;   oplot, lag3, normalize(xa3[2,*]), color=rgb(255,255,0)
;   oplot, lag3, normalize(xa3[3,*]), color=rgb(0,255,0)
;   plot, lag4, normalize(xa4[0,*]), xst=1, xr=[-180,180], charsize=2
;   oplot, lag4, normalize(xa4[1,*]), color=rgb(255,0,0)
;   oplot, lag4, normalize(xa4[2,*]), color=rgb(255,255,0)
;   oplot, lag4, normalize(xa4[3,*]), color=rgb(0,255,0)
;   plot, lag5, normalize(xa5[0,*]), xst=1, xr=[-180,180], charsize=2
;   oplot, lag5, normalize(xa5[1,*]), color=rgb(255,0,0)
;   oplot, lag5, normalize(xa5[2,*]), color=rgb(255,255,0)
;   oplot, lag5, normalize(xa5[3,*]), color=rgb(0,255,0)
;   !p.multi=[0,1,0]


d_abcd_radW = (dabcdW * (2.d0*!DPI))/lambdaW
d_abcd_rad1 = (dabcd1 * (2.d0*!DPI))/lambda1
d_abcd_rad2 = (dabcd2 * (2.d0*!DPI))/lambda2
d_abcd_rad3 = (dabcd3 * (2.d0*!DPI))/lambda3
d_abcd_rad4 = (dabcd4 * (2.d0*!DPI))/lambda4
d_abcd_rad5 = (dabcd5 * (2.d0*!DPI))/lambda5

angradW = [d_abcd_radW[1]-(!DPI/2.d0), d_abcd_radW[2]-(!DPI), d_abcd_radW[3]-((3.d0*!DPI)/2.d0)]
angrad1 = [d_abcd_rad1[1]-(!DPI/2.d0), d_abcd_rad1[2]-(!DPI), d_abcd_rad1[3]-((3.d0*!DPI)/2.d0)]
angrad2 = [d_abcd_rad2[1]-(!DPI/2.d0), d_abcd_rad2[2]-(!DPI), d_abcd_rad2[3]-((3.d0*!DPI)/2.d0)]
angrad3 = [d_abcd_rad3[1]-(!DPI/2.d0), d_abcd_rad3[2]-(!DPI), d_abcd_rad3[3]-((3.d0*!DPI)/2.d0)]
angrad4 = [d_abcd_rad4[1]-(!DPI/2.d0), d_abcd_rad4[2]-(!DPI), d_abcd_rad4[3]-((3.d0*!DPI)/2.d0)]
angrad5 = [d_abcd_rad5[1]-(!DPI/2.d0), d_abcd_rad5[2]-(!DPI), d_abcd_rad5[3]-((3.d0*!DPI)/2.d0)]

b0 = [angradW[0],angrad1[0],angrad2[0],angrad3[0],angrad4[0],angrad5[0]]
c0 = [angradW[1],angrad1[1],angrad2[1],angrad3[1],angrad4[1],angrad5[1]]
d0 = [angradW[2],angrad1[2],angrad2[2],angrad3[2],angrad4[2],angrad5[2]]

phserrA = sin(b0)
phserrB = 1.d0+cos(c0)
phserrC = cos(b0)+cos(d0)
phserrD = -sin(b0)-sin(d0)
phserr = [[phserrA],[phserrB],[phserrC],[phserrD]]

angvalW = [90.d0+angradW[0]/!dtor,180.d0+angradW[0]/!dtor,270.d0+angradW[2]/!dtor]
angval1 = [90.d0+angrad1[0]/!dtor,180.d0+angrad1[0]/!dtor,270.d0+angrad1[2]/!dtor]
angval2 = [90.d0+angrad2[0]/!dtor,180.d0+angrad2[0]/!dtor,270.d0+angrad2[2]/!dtor]
angval3 = [90.d0+angrad3[0]/!dtor,180.d0+angrad3[0]/!dtor,270.d0+angrad3[2]/!dtor]
angval4 = [90.d0+angrad4[0]/!dtor,180.d0+angrad4[0]/!dtor,270.d0+angrad4[2]/!dtor]
angval5 = [90.d0+angrad5[0]/!dtor,180.d0+angrad5[0]/!dtor,270.d0+angrad5[2]/!dtor]

angval = [[angvalW], [angval1], [angval2], [angval3], [angval4], [angval5]] mod 360.d0

; print, ''
; print, 'Values FFT and CC'
; print, wlen
; print, angval


;Wavelength Quality Control
qcwlen = strarr(4,6)
qcfitwlen = strarr(4,6)
for i=0,3 do begin
  for j=0,5 do begin

    if (j eq 0) then begin
      if (wlen[i,j] gt 2.2d-6 and wlen[i,j] lt 2.3d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
      if (fitwlen[i,j] gt 2.2d-6 and fitwlen[i,j] lt 2.3d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 1) then begin
      if (wlen[i,j] gt 1.95d-6 and wlen[i,j] lt 2.1d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
      if (fitwlen[i,j] gt 1.95d-6 and fitwlen[i,j] lt 2.1d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 2) then begin
      if (wlen[i,j] gt 2.05d-6 and wlen[i,j] lt 2.2d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
      if (fitwlen[i,j] gt 2.05d-6 and fitwlen[i,j] lt 2.2d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 3) then begin
      if (wlen[i,j] gt 2.15d-6 and wlen[i,j] lt 2.3d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
      if (fitwlen[i,j] gt 2.15d-6 and fitwlen[i,j] lt 2.3d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 4) then begin
      if (wlen[i,j] gt 2.3d-6 and wlen[i,j] lt 2.4d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
      if (fitwlen[i,j] gt 2.3d-6 and fitwlen[i,j] lt 2.4d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 5) then begin
      if (wlen[i,j] gt 2.35d-6 and wlen[i,j] lt 2.5d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
      if (fitwlen[i,j] gt 2.35d-6 and fitwlen[i,j] lt 2.5d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

  endfor
endfor


;Phase Shift Quality Control
qcangval = strarr(3,6)
qcfitangval = strarr(3,6)
for i=0,2 do begin
  for j=0,5 do begin

    if (i eq 0) then begin
      if (angval[i,j] ge 65.d0 and angval[i,j] le 85.d0) then qcangval[i,j] = '+1' else qcangval[i,j] = '-1'
      if (fitangval[i,j] ge 65.d0 and fitangval[i,j] le 85.d0) then qcfitangval[i,j] = '+1' else qcfitangval[i,j] = '-1'
    endif

    if (i eq 1) then begin
      if (angval[i,j] gt 155.d0 and angval[i,j] lt 175.d0) then qcangval[i,j] = '+1' else qcangval[i,j] = '-1'
      if (fitangval[i,j] gt 155.d0 and fitangval[i,j] lt 175.d0) then qcfitangval[i,j] = '+1' else qcfitangval[i,j] = '-1'
    endif

    if (i eq 2) then begin
      if (angval[i,j] gt 215.d0 and angval[i,j] lt 250.d0) then qcangval[i,j] = '+1' else qcangval[i,j] = '-1'
      if (fitangval[i,j] gt 215.d0 and fitangval[i,j] lt 250.d0) then qcfitangval[i,j] = '+1' else qcfitangval[i,j] = '-1'
    endif

  endfor
endfor



; ;output
; 
; openw, lun, resultpath+'FSUA_SkyCalibration_ABCD.dat', width=1400, /get_lun
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, phserr[i,j], format='(f15.10)'
;     endfor
;   endfor
; close, lun
; free_lun, lun
; 
; openw, lun, resultpath+'FSUA_SkyCalibration_Angles.dat', width=1400, /get_lun
;   printf, lun, 'Fit Fringe values'
;   for i=0,5 do printf, lun, fitangval[*,i], format='(3f16.10)'
;   for i=0,5 do printf, lun, qcfitangval[*,i], format='(3a10)'
;   printf, lun, ''
;   printf, lun, 'Classical Values (CC)'
;   for i=0,5 do printf, lun, angval[*,i], format='(3f15.10)'
;   for i=0,5 do printf, lun, qcangval[*,i], format='(3a10)'
; close, lun
; free_lun, lun
; 
; openw, lun, resultpath+'FSUA_SkyCalibration_Lambda.dat', width=1400, /get_lun
;   printf, lun, 'Fit Fringe values'
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, fitwlen[j,i], format='(e20.10)'
;     endfor
;   endfor
;   printf, lun, ''
;   printf, lun, 'Classical Values (CC)'
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, wlen[j,i], format='(e20.10)'
;     endfor
;   endfor
; close, lun
; free_lun, lun
; 
; openw, lun, resultpath+'FSUA_Calibration_Waves.dat', width=1400, /get_lun
;   printf, lun, 'Fit Fringe values'
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, fitwlen[j,i], qcfitwlen[j,i], format='(e20.10, a5)'
;     endfor
;   endfor
;   printf, lun, ''
;   printf, lun, 'Classical Values (CC)'
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, wlen[j,i], qcwlen[j,i], format='(e20.10, a5)'
;     endfor
;   endfor
; close, lun
; free_lun, lun
; 
; 
; 
print, ''
print, 'Wavelength and Phase values from fitting the Fringe package'
print, '         A                     B                     C                     D'
for i=0,5 do print, string(fitwlen[*,i])+'   '+string(qcfitwlen[*,i]) 
print, ''
for i=0,5 do print, string(fitangval[*,i]);+'   '+string(qcfitangval[*,i]) 

print, ''
print, 'Wavelength and Phase values from FFT and CC'
print, '         A                     B                     C                     D'
for i=0,5 do print, string(wlen[*,i])+'   '+string(qcwlen[*,i]) 
print, ''
for i=0,5 do print, string(angval[*,i])+'   '+string(qcangval[*,i]) 



stop
return
end
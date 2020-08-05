MADE SOME MODIFICATION IN MIDI_FSUA_P92RIVI/MIDI_UT

@ploterror.pro
@symcat.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@fsc_color.pro
@tvread.pro
@interpol.pro
@int_tabulated.pro
@uniq.pro
@errplot.pro
@mpfitfun.pro
@mpfit.pro
@oploterror.pro
@setdefaultvalue.pro
@tag_exist.pro
@cgquery.pro
@cgplot.pro
@setdecomposedstate.pro
@decomposedcolor.pro
@cgdefaultcolor.pro
@colorsareidentical.pro
@cgcolor.pro
@cgsnapshot.pro
@cgdefcharsize.pro
@getdecomposedstate.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@reverse.pro
@showsym.pro

;2MASS http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
; WISE mags + 2MASS von cal + sci (consider bandwidth, filterkurven, zeropoints)
;   an cals blackbody fitten
; 
;   fcorr_sci/fcorr_cal * blackbody fit
; 
; http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/311/wise&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
; 
; WISE 	http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html#PhotometricZP
;	http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/figures/sec4_3gt4.gif
;
; 2MASS http://www.ipac.caltech.edu/2mass/releases/allsky/faq.html#jansky
; http://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzermission/missionoverview/spitzertelescopehandbook/19/

function func_planck, x, p, _extra=struc

  ;based on planck.pro

  c = 299792458.d0	;m/s
  h = 6.62606957d-34	;J*s
  k = 1.3806488d-23	;J/K
  p2m = 3.0856776d16	;parsec to meter

  c1 = 2.*!DPI*h*c*c
  c2 = h*c/k

  c1 = c1*1.d11	;bring to the exponent similar to planck.pro
  c2 = c2*1.d2	;bring to the exponent similar to planck.pro

;----------------------------------------

  ;stellar radius in meter
  plx = struc.plx
  diam = struc.diam
  d = p2m/plx	;distance in meter
  R = d*tan((diam/2.d0/3600.d0)*!dtor)
  c3 = R^2./d^2.

;----------------------------------------

  val =  c2/x/p[0]
  bbflux =  c1 / ( x^5 * ( exp(val)-1. ) )
  bbflux = bbflux*1.d-8	;ergs/cm2/s/A

  bbflux = bbflux*c3

  fit = bbflux*(x/1.d-8)^2.*( 3.33564095d4)	;conversion to Jy from http://dept.astro.lsa.umich.edu/~dmaitra/conversion.html

  return, fit

end

function plot_func_planck, x, p, diam, plx

  ;based on planck.pro

  c = 299792458.d0	;m/s
  h = 6.62606957d-34	;J*s
  k = 1.3806488d-23	;J/K
  p2m = 3.0856776d16	;parsec to meter

  c1 = 2.*!DPI*h*c*c
  c2 = h*c/k

  c1 = c1*1.d11	;bring to the exponent similar to planck.pro
  c2 = c2*1.d2	;bring to the exponent similar to planck.pro

;----------------------------------------

  ;stellar radius in meter
  d = p2m/plx	;distance in meter
  R = d*tan((diam/2.d0/3600.d0)*!dtor)
  c3 = R^2./d^2.

;----------------------------------------

  val =  c2/x/p[0]
  bbflux =  c1 / ( x^5 * ( exp(val)-1. ) )
  bbflux = bbflux*1.d-8	;ergs/cm2/s/A

  bbflux = bbflux*c3

  plotfit = bbflux*(x/1.d-8)^2.*( 3.33564095d4)	;conversion to Jy from http://dept.astro.lsa.umich.edu/~dmaitra/conversion.html

  return, plotfit

end

function fint, xfilt, filt, xspec, spec

; this functions returns the integral of SPEC over FILT
;
; XFILT, XSPEC can be either wavelength or frequency, but should be the same domain
;
; NOTE: FILT is the RSR (relative system response) in normalized flux(lambda) units
;      -> SPEC should be in F_lam as well!!!
;

  ; make sure that both curves, use the same base; 
  ; it seems adequate to interpolate to the filter sampling
  return, int_tabulated(xfilt,interpol(spec,xspec,xfilt)*filt) / int_tabulated(xfilt,filt)

end


pro totflux

  obsnum = ''
;   read, 'Observation ID: ', obsnum
  obsnum = 'P92CHOQUET'

  star = ''
  read, 'Star: ', star
;star = 'HD76111'
;======================================================================================

if (star eq 'HD76111') then begin

  T = 4660.d0
  plx = 3.70d-3		;arcsec, Hipparcos
  diam = 0.411d-3	;arcsec, from II/300

endif

if (star eq 'HD80934') then begin

  T = 4660.d0
  plx = 7.40d-3		;arcsec, Hipparcos
  diam = 0.467d-3	;arcsec, from II/300

endif

struc = {plx:plx, diam:diam}


;======================================================================================

  wldir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_EWS_CustomMask/'
  resultdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/WISE_Photometry/'
  file_mkdir, resultdir


  micron = '!Mm!X'+'m'

;======================================================================================

;filter definitions

  filter = ['J', 'H', 'K', 'W1', 'W2', 'W3', 'W4']
  nf  = n_elements(filter)
  filtdata = dblarr(nf, 6)

;			        Wl mag mage Fli Fni  FniErr muiso
  readcol, 'data_'+star+'.txt', wl, d1, d2, d3, d4, d5, d6, format='d,d,d,d,d,d,d', /silent
  unit = [micron, 'mag', 'mag', 'W cm^-2 '+micron+'^-1', 'Jy', 'Jy', 'Hz']

  zeromag = [[wl], [d1], [d2], [d3], [d4], [d5], [d6]]

  ;RSR for filters
  readcol, '2MASS_RSR/J_band.txt', l, d, format='d,d', /silent
  jfilt = transpose([[l], [d]])

  readcol, '2MASS_RSR/H_band.txt', l, d, format='d,d', /silent
  hfilt = transpose([[l], [d]])

  readcol, '2MASS_RSR/K_band.txt', l, d, format='d,d', /silent
  kfilt = transpose([[l], [d]])

  readcol, 'WISE_RSR/RSR-W1.EE.txt', l, d1, d2, format='d,d,d', /silent
  w1filt = transpose([[l], [d1]])

  readcol, 'WISE_RSR/RSR-W2.EE.txt', l, d1, d2, format='d,d,d', /silent
  w2filt = transpose([[l], [d1]])

  readcol, 'WISE_RSR/RSR-W3.EE.txt', l, d1, d2, format='d,d,d', /silent
  w3filt = transpose([[l], [d1]])

  readcol, 'WISE_RSR/RSR-W4.EE.txt', l, d1, d2, format='d,d,d', /silent
  w4filt = transpose([[l], [d1]])

;======================================================================================

;plotting, RSR filters

window,0,ysize=800
!p.multi=[0,0,3]

  plot,[0],[0], xr=[1,30], yr=[0,1], /xst, /yst, title='Relative Spectral Response (RSR) curves', xtit='Wavelength / ['+micron+']', /xlog, charsize=3
    oplot, jfilt[0,*], jfilt[1,*]
    oplot, replicate(zeromag[0,0], 2), !y.crange, color=fsc_color('violet')
    oplot, hfilt[0,*], hfilt[1,*]
    oplot, replicate(zeromag[1,0], 2), !y.crange, color=fsc_color('blue')
    oplot, kfilt[0,*], kfilt[1,*]
    oplot, replicate(zeromag[2,0], 2), !y.crange, color=fsc_color('turquoise')
    oplot, w1filt[0,*], w1filt[1,*]
    oplot, replicate(zeromag[3,0], 2), !y.crange, color=fsc_color('green')
    oplot, w2filt[0,*], w2filt[1,*]
    oplot, replicate(zeromag[4,0], 2), !y.crange, color=fsc_color('yellow')
    oplot, w3filt[0,*], w3filt[1,*]
    oplot, replicate(zeromag[5,0], 2), !y.crange, color=fsc_color('orange')
    oplot, w4filt[0,*], w4filt[1,*]
    oplot, replicate(zeromag[6,0], 2), !y.crange, color=fsc_color('red')

;======================================================================================

  valuesl = (10.^(-[zeromag[0,1], zeromag[1,1], zeromag[2,1], zeromag[3,1], zeromag[4,1], zeromag[5,1], zeromag[6,1]]/2.5d0))*zeromag[*,3]
  errorsl = (10.^(-[zeromag[0,2], zeromag[1,2], zeromag[2,2], zeromag[3,2], zeromag[4,2], zeromag[5,2], zeromag[6,2]]/2.5d0))

  nsteps = 1000
  lmin = 1d & lmax = 30d ; in micron
  lam = lmin*((lmax/lmin)^(1./nsteps))^dindgen(nsteps+1)  ; logarithmic lam-scaling from 1um to 30um seems adequate

;   valuesn = (10.^(-[zeromag[0,1], zeromag[1,1], zeromag[2,1], zeromag[3,1], zeromag[4,1], zeromag[5,1], zeromag[6,1]]/2.5d0))*zeromag[*,4]
;   errorsn = (10.^(-[zeromag[0,2], zeromag[1,2], zeromag[2,2], zeromag[3,2], zeromag[4,2], zeromag[5,2], zeromag[6,2]]/2.5d0))

  valuesn = (10.^(-zeromag[*,1]/2.5d0))*zeromag[*,4]
;   errorsn = (10.^(-zeromag[*,2]/2.5d0))
  tmpe1 = (10.^(-(zeromag[*,2]/zeromag[*,1])/2.5d0))

  ;error are not symmetric because of 10^x, however take the largest of the two limits
;   tmp = (10.^(-(zeromag[*,1])/2.5d0))

  tmpe = dblarr(n_elements(zeromag[*,2]))
  tmp1 = abs((valuesn*tmpe1)-valuesn)
  tmp2 = abs((valuesn/tmpe1)-valuesn)
  for i=0,n_elements(tmpe1)-1 do tmpe[i] = max([tmp1[i], tmp2[i]])

  ;include the error of Fni, isophotic flux
  ;x = [zeromag[0,1], zeromag[1,1], zeromag[2,1], zeromag[3,1], zeromag[4,1], zeromag[5,1], zeromag[6,1]]
  valuesnerr = valuesn*sqrt((tmpe/valuesn)^2. + (zeromag[*,5]/zeromag[*,4])^2.)


  lam4 = lam^(-4.) 
  lam4 = lam4*valuesl[5]/fint(w3filt[0,*],w3filt[1,*], lam, lam4)

  plot, lam, lam4, ytit='[F_lam]', xtit='Wavelength / ['+micron+']', xr=[1,30], title='lam^-4 fit to the MIR', charsize=3, /ylog;, /xlog
  oplot, zeromag[*,0], valuesl, psym=sym(1)
  errplot, zeromag[*,0], valuesl*errorsl, valuesl/errorsl

;======================================================================================

  ;fit WISE data points with a black body
  lamCM = 299792458.d0/zeromag[*,6]*1d2
;   e1 = valuesn*errorsn
;   e2 = valuesn/errorsn
;   ndata = n_elements(e1)
;   error = dblarr(ndata)
;   for i=0,ndata-1 do error[i] = max([e1[i], e2[i]])
;   valuesnerr = abs(valuesn-error)

  weight = 1.d0/(valuesnerr/valuesn)^2.
;   error[*] = 1.
  weight[*] = 1.

  start_val = T
  params = mpfitfun('func_planck', lamCM, valuesn, valuesnerr, start_val, weights=weight, perror=perror, $
	  maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, functargs=struc, yfit=yfit, /quiet)
  print, params, perror

  ;for plotting
  plotlam = lam/1.d4
  plotfit = plot_func_planck(plotlam, params, diam, plx)
  plotfiterr1 = abs(plot_func_planck(plotlam, params+perror, diam, plx)-plotfit)
  plotfiterr2 = abs(plot_func_planck(plotlam, params-perror, diam, plx)-plotfit)
  plotfiterr = plotfiterr1
  for i=0,n_elements(plotfiterr1)-1 do plotfiterr[i] = max([plotfiterr1[i], plotfiterr2[i]])

;======================================================================================

  ;plot the interesting N-band in Jansky

  flux_jy = lam4*lam^2/299792458.d0*1.d24

  lamMIC = 299792458.d0/zeromag[*,6]*1.d6
;   plot,  lamMIC, valuesn, psym=sym(1), ytit='[Jy]', xtit='Wavelength ['+micron+']', xr=[1,30], /ylog, tit='lam^-4 (green) / Planck Fit (red)', charsize=3
;   errplot, lamMIC, valuesn*errorsn, valuesn/errorsn
;   oplot, lam, flux_jy, color=fsc_color('green')
;   oploterror, lam, plotfit, plotfiterr, color=fsc_color('red'), /nohat, errcolor=fsc_color('red')
;   oplot,  lamMIC, valuesn, psym=sym(1)
;   errplot, lamMIC, valuesn*errorsn, valuesn/errorsn
  ploterror,  lamMIC, valuesn, valuesnerr, psym=sym(1), ytit='[Jy]', xtit='Wavelength ['+micron+']', xr=[1,30], /ylog, tit='lam^-4 (green) / Planck Fit (red)', charsize=3, background=fsc_color('black')
  oploterror, lam, plotfit, plotfiterr, color=fsc_color('red'), /nohat, errcolor=fsc_color('red')
  oplot,  lamMIC, valuesn, psym=sym(1)
  oploterror,  lamMIC, valuesn, valuesnerr
  oplot, lam, flux_jy, color=fsc_color('green')

;======================================================================================

  ;get wavelength from MIDI

  file = file_search(wldir+star+'_*.corr.fits')
  file = file[0]

  st = mrdfits(file, 1, /silent)
  wlmidi = reverse(st.eff_wave)*1.d6

;   intfl = interpol(flux_jy, lam, wlmidi, /spline)

  tmpwl = wlmidi/1.d4
  fit = plot_func_planck(tmpwl, params, diam, plx)
  fiterr1 = abs(plot_func_planck(tmpwl, params+perror, diam, plx)-fit)
  fiterr2 = abs(plot_func_planck(tmpwl, params-perror, diam, plx)-fit)
  fiterr = fiterr1
  for i=0,n_elements(fiterr1)-1 do fiterr[i] = max([fiterr1[i], fiterr2[i]])


;   flux = intfl
  flux = fit
  fluxerr = fiterr

  wlmidi = wlmidi/1.d6

  fn = resultdir+star+'_WISE_Photometry.sav'
  save,filename=fn, wlmidi, flux, fluxerr
  print, 'SAVED fluxes in : '+ fn


!p.multi=[0,1,0]


stop
End
 
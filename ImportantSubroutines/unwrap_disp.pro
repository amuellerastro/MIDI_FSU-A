;from https://www.cfa.harvard.edu/~cqi/mircook.html

;------------------------------------------------------------------------------

function uti_pha_180, phas
;
; Converts phase data to be in range +-180 deg.
;
; parameters : phas -- input phase array
; result  =  modified phases
; eg. : phas = uti_pha_180(phas)
;

; next line adds 4 turns to the phase and then mods with 360.
  phas_180 = phas
  phas_180 = (1440. + phas) mod 360. 
  j = where(phas_180 gt 180.) 
  if max(j) ne -1 then phas_180[j] = phas_180[j]-360.

  return, phas_180

end

;------------------------------------------------------------------------------

function uti_pha_unwrap,pha,smooth = smooth,ramp = ramp  
;
; Unwrap phases.
;
; The data are converted to complex (with unity radius), the real
; and imaginary are smoothed and this function
; is used to track the phase through +- 180 deg.
;
; parameters : pha    -- phases to be unwrapped
;      keyword smooth -- sets number of pts to smooth over
;                        default = 3
;              ramp   -- expected phase ramp between adjacent points
;                        default = 10.
;
; result  =  -1 (error), unwrapped phase array (succesful)
;
; to generate test data :
; eg  : pha  = 20*findgen(115)
; eg. : for i = 0,n_elements(pha)-1 do begin & $
; eg. :  pha(i) = pha(i)+50*randomu(seed) & $
; eg. : endfor
; eg. : pha = uti_pha_180(pha)
;
; eg. : result = uti_pha_unwrap(pha)
;
; eg. : plot,findgen(30),pha
;
  pha_all = pha
  pha = uti_pha_180(pha)
  npts = n_elements(pha)
  ;
  ; first compute st_dev of pha from derivative of (phases - smoothed phases)
  ;
  filter = make_array(3.,/float,value = 1./3.)
  smo = pha-convol(pha,filter,/edge_truncate,/nan)
  j = where(abs(smo) lt 70.,count)
  if count gt 2 then smo = smo[j] 
  npts_smo = n_elements(smo)
  result = moment(smo-[smo(1:npts_smo-1),smo(npts_smo-1)],sdev = sdev,/nan)
  ;
  ; convert phases to re and im vis since these one can smooth
  ; w/o smoothing through discontinuities (since the complex
  ; vis is a continuous function
  ;
  uti_conv_apc,cmp,make_array(npts,/float,value = 1.),pha,/complex
  re = float(cmp) & im = imaginary(cmp)
  ;
  ; smooth the complex vis and convert back to amp,pha
  ;
  if not keyword_set(ramp) then ramp = 10.
  if not keyword_set(smooth) then smooth = (sdev/float(ramp))^2
  smooth = long(smooth)
  if (smooth mod 2L) ne 1L then smooth = smooth+1L
  if smooth gt (n_elements(pha)-2L) then smooth = n_elements(pha)-2L
  smooth = max([smooth,1L])
  filter = make_array(smooth,/float,value = 1./smooth)
  uti_conv_apc,complex(convol(re,filter,/edge_truncate,/nan), $
	convol(im,filter,/edge_truncate,/nan)),ampnf,phaf,/amp_pha
  pha_jmp = make_array(npts,/float,value = 0.)
  ;
  ; compare the smoothed phases with themselves shifted by one point
  ; ie. look for jumps
  ;
  diff = phaf-[phaf(1:npts-1),phaf(npts-1)]
  ; 
  ; first find correct lobes for smoothed phases
  ;
  ; jp  = > positive jump , jm  = > negative jump
  ;
  jp = where(diff gt 180.) & jm = where(diff lt -180.)
  if max(jp) ne -1 then begin

    for i = 0,n_elements(jp)-1 do begin  

      pha_jmp(jp[i]+1:npts-1) =  pha_jmp(jp[i]+1:npts-1) +360. 

    endfor 

  endif

  if max(jm) ne -1 then begin 

    for i = 0,n_elements(jm)-1 do begin 

      pha_jmp(jm[i]+1:npts-1) =  pha_jmp(jm[i]+1:npts-1) -360. 

    endfor 

  endif
  ;
  ;  put phase jumps into smoothed phases
  ;
  phaf = phaf+pha_jmp
  pha = pha+pha_jmp
  ;
  ; now make sure actual phases follow lobes of smoothed phases
  ;
  diff = pha-phaf
  pha_jmp = make_array(npts,/float,value = 0.)
  jp = where(diff gt 180.) & jm = where(diff lt -180.)

  if max(jp) ne -1 then begin 

    for i = 0,n_elements(jp)-1 do begin  

      pha_jmp(jp[i]) =  pha_jmp(jp[i]) -360.

    endfor 

  endif

  if max(jm) ne -1 then begin 

    for i = 0,n_elements(jm)-1 do begin 

      pha_jmp(jm[i]) =  pha_jmp(jm[i]) +360. 

    endfor 

  endif
  ;
  ; put in phase jumps and set mean value to lobe closest to zero
  ;
  pha = pha+pha_jmp
  ;
  lobe = float(round(total(pha,/nan)/max([n_elements(pha),1],/nan)/360.)) 
  phas = pha-lobe*360.
  pha = pha_all & pha = phas

  return, pha

end

;------------------------------------------------------------------------------

pro  uti_conv_apc,cmp,amp,pha,amp_pha = amp_pha,complex = complex 
;
; Converts between complex and amp/pha data  
; 
;
; parameters : cmp      -- complex array 
;              amp      -- amp array 
;              pha      -- pha array 
;              amp_pha  -- keyword to compute amp and phase
;              complex  -- keyword to compute complex
;
; eg. : uti_conv_apc,cmp,amp,pha,/complex ; convert amp,pha to complex 
; eg. : uti_conv_apc,cmp,amp,pha,/amp_pha ; convert complex to amp,pha
;
if keyword_set(amp_pha) then begin

  amp = abs(cmp) & pha = !radeg*atan(imaginary(cmp),float(cmp))

endif

if keyword_set(complex) then begin

  pha_rad = pha*!dtor
  cmp = amp*complex(cos(pha_rad),sin(pha_rad))

  endif

end

;------------------------------------------------------------------------------

pro unwrap_disp, star, miditime, resultpath, timeorig, pdelorig, gdelorig, statorig, lamorig, $
	hak=hak,box=box,xr=xr,sav=sav;, ti=ti

  if keyword_set(hak) then hakk=1 else hakk=0
  if not keyword_set(box) then box = 2001
  if (box mod 2) ne 1 then box=box+1	; make box odd for the median calculation
  if keyword_set(xr) then xr = double(xr) else xr = [0d,0d]

  time = timeorig
  pdel = pdelorig
  gdel = gdelorig
  stat = statorig
  lam = lamorig

  flag = intarr(n_elements(time))	;flag to display fringe lock

  ;backup arrays to have original array size but filled with NANs where needed
  time_bu = time
  pdel_bu = pdel
  gdel_bu = gdel
  idx1 = where(stat lt 7)
  if (idx1[0] gt -1) then begin
 
    pdel_bu[idx1] = !values.d_nan
    gdel_bu[idx1] = !values.d_nan
    flag[idx1] = !values.d_nan

  endif

  ;select only good data where FSUA locked on fringe
  idx7 = where(stat eq 7)
  if (idx7[0] gt -1) then begin
 
    time = time[idx7]
    pdel = pdel[idx7]
    gdel = gdel[idx7]
    flag[idx7] = 3

  endif

  tmp = pdel/!dtor
  up_deg = uti_pha_unwrap(tmp, smooth=box)
  up = up_deg*!dtor

  pdel = pdel/(2.d0*!DPI)*lam
  up = up/(2.d0*!DPI)*lam

  if hakk then begin

    micron = string(181B)
    window, 0, xs=1200, ys=1000, title=star+' '+miditime
    !p.multi=[0,1,2]

    plot, time/1.d6, pdel*1.d6, psym=3, $
      xr=[0,time[n_elements(time)-1]/1.d6], xst=1, yr=[-2*lam*1d6,2*lam*1d6], yst=1, $
      background=fsc_color('white'), color=fsc_color('black'), $
      xtitle='Time / [s]', ytitle='Delay ['+micron+'m]', charsize=2, $
      title='orig. PD (black) / unwrapped PD (red) / GD (blue)'
    plots, !x.crange, [lam*1d6, lam*1d6], color=fsc_color('dark gray');, linestyle=2
    plots, !x.crange, [0,0], color=fsc_color('dark gray');, linestyle=2
    plots, !x.crange, [-lam*1d6, -lam*1d6], color=fsc_color('dark gray');, linestyle=2
    oplot, time/1.d6, up*1d6, color=fsc_color('red'), psym=3
    oplot, time/1.d6, gdel*1d6, color=fsc_color('blue'), psym=3
    oplot, time/1.d6, flag, psym=3, color=fsc_color('dark green')

    plot, time/1.d6, pdel*1.d6, psym=3, $
      xr=[0,time[n_elements(time)-1]/1.d6], xst=1, yr=[-4,4], yst=1, $
      background=fsc_color('white'), color=fsc_color('black'), $
      xtitle='Time / [s]', ytitle='Delay ['+micron+'m]', charsize=2, $
      title='orig. PD (black)'

    oplot, time/1.d6, flag, psym=3, color=fsc_color('dark green')

    plots, [25,25], [-4,20], color=fsc_color('green')
    plots, [50,50], [-4,20], color=fsc_color('green')
    plots, [75,75], [-4,20], color=fsc_color('green')
    plots, [100,100], [-4,20], color=fsc_color('green')
    plots, [125,125], [-4,20], color=fsc_color('green')
    plots, [150,150], [-4,20], color=fsc_color('green')
    plots, [175,175], [-4,20], color=fsc_color('green')

    !p.multi=[0,1,0]
  endif

  set_plot, 'ps'
  device, filename=resultpath+star+'_'+miditime+'_unwrap.ps',/color,XSIZE=30, YSIZE=13, XOffset=xoffset, YOffset=yoffset
  !p.thick=4
  !x.thick=3
  !y.thick=3

    plot, time/1.d6, pdel*1.d6, psym=3, $
      xr=[0,time[n_elements(time)-1]/1.d6], xst=1, yr=[-2*lam*1d6,2*lam*1d6], yst=1, $
      background=fsc_color('white'), color=fsc_color('black'), $
      xtitle='Time / [s]', ytitle='Delay ['+micron+'m]', charsize=1.5, $
      title='orig. PD (black) / unwrapped PD (red) / GD (blue)', $
      position=[0.07,0.13,0.96,0.91]
    plots, !x.crange, [lam*1d6, lam*1d6], color=fsc_color('dark gray');, linestyle=2
    plots, !x.crange, [0,0], color=fsc_color('dark gray');, linestyle=2
    plots, !x.crange, [-lam*1d6, -lam*1d6], color=fsc_color('dark gray');, linestyle=2
    oplot, time/1.d6, up*1d6, color=fsc_color('red'), psym=3
    oplot, time/1.d6, gdel*1d6, color=fsc_color('blue'), psym=3

  device,/close
  set_plot,'x'
  !p.multi=[0,1,0]
  !p.thick=1
  !x.thick=1
  !y.thick=1


  ;fill backup arrays with corrected data but keep NANs and rename
  idx = where(finite(pdel_bu) ne 0)
  pdel_bu[idx] = up
  gdel_bu[idx] = gdel
  time = time_bu
  pdel = pdel_bu
  gdel = gdel_bu

  disp = pdel-gdel	;up is unwraped phase

  if keyword_set(sav) then begin

    fn = resultpath+star+'_'+miditime+'_unwrap.sav'
    print,'Saving final data in '+fn
    save,filename=fn,time,disp,pdel,gdel,box

  endif

end
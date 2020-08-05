@butterworth.pro
@real_part.pro

Function edge_med,arr,box

; this function reflects the edges of ARR by box 
; to get better median estimates at the edges

  if n_elements(box) gt n_elements(arr) then stop,"BOX is gt ARR, this doesn't work, stopping here..."

  arr = double(arr)
  box = long(box)

  narr = n_elements(arr)
  tarr = [arr[reverse(lindgen(box))],arr,arr[narr-1-lindgen(box)]]
  return,(median(tarr,box))[box+0:box+narr-1]
End

Function edge_smo,arr,box

; this function reflects the edges of ARR by box 
; to get better smooth estimates at the edges
;
; we use smooth ALWAYS with option /NAN

  if n_elements(box) gt n_elements(arr) then stop,"BOX is gt ARR, this doesn't work, stopping here..."

  narr = n_elements(arr)
  tarr = [arr[reverse(lindgen(box))],arr,arr[narr-1-lindgen(box)]]
  return,(smooth(tarr,box,/nan))[box+0:box+narr-1]
End

Function defgbidx,narr,bidx,box

; INPUT: NARR is the total length of the data array with 
;        a BIDX index indicating invalid points, 
;        BOX identifies the minimum length of BIDX pixel 
;         if there is data areas of lt BOX length between two BIDX areas 
;         then this good data area is also classified as bad
;         to enable reasonable BOX smoothing
;        BOX is typically the array smoothing length
;        BIDX is typically the result of BIDX = where(arr eq (bad data identifier))
;        all *IDX are indices of arr
; 
; RETURN: a structure with 
;        GIDX: the good data point index
;        BIDX: the updated bidx, the minimum length of continuous BIDX is now BOX
;        GIDXS,GIDXE: starting / ending of continuos GIDX areas
;        BIDXS,BIDXE: starting / ending of continuos BIDX areas
;

  narr = long(narr)
  bidx = long(bidx)
  box  = long(box)


  if box gt narr then stop,'BOX is larger than ARRAY, stopping here...'
  if (bidx)[0] eq (-1) then stop,'BIDX=-1 does not make sense here...'
  
  case n_elements(bidx) of 
    1: $
       begin
        bidxs = [bidx]
        bidxe = [bidx]
       end
    else:  $
      begin
        bidxe = where(ts_diff(bidx,1) lt -1)  ; this finds all bidxe, but the last
                                              ; -> this is -1 in case there is 
					      ; only one continuous strip of BIDXes
        bidxs = bidx[[0,bidxe+1]]           
        bidxe = (bidxe[0] eq (-1))?bidx[n_elements(bidx)-1]:bidx[[bidxe,(n_elements(bidx)-1)]]
      end
  endcase

  ; decrease/increase BIDX around BIDXS/E by half box size
  nbidx = bidx
  for i=0,n_elements(bidxe)-1 do $
       nbidx = [nbidx,[lindgen((bidxe[i]+floor(box/2))-(bidxs[i]-floor(box/2))+1)+(bidxs[i]-floor(box/2))]]
  nbidx = nbidx[sort(nbidx)]
  nbidx = nbidx[uniq(nbidx)]
  nbidx = nbidx[where((nbidx ge 0) and (nbidx le (narr-1)))]

  bidxe = where(ts_diff(nbidx,1) lt -1)
  bidxs = nbidx[[0,[bidxe+1]]]           
  bidxe = (bidxe[0] eq (-1))?nbidx[n_elements(nbidx)-1]:nbidx[[bidxe,(n_elements(nbidx)-1)]]

  bidxs = bidxs[sort(bidxs)]
  bidxs = bidxs[uniq(bidxs)]
  bidxe = bidxe[sort(bidxe)]
  bidxe = bidxe[uniq(bidxe)]


  ; check if now there is data groups of less/equal BOX length
  tidx =  where(([bidxs,[narr]] - [[0],bidxe]) lt (box+2))
  if (tidx)[0] ge 0 then begin
    ; fill bidxs
    for i = 0,n_elements(tidx)-1 do $
      if tidx[i] eq 0 then bidxs[0] = 0 $
        else if tidx[i] ne n_elements(bidxs) then $
         bidxs[tidx[i]] = bidxs[tidx[i]-1]

    ; fill bidxe in reverse order
    for i =n_elements(tidx)-1,0,-1 do $
      if tidx[i] eq n_elements(bidxs) then bidxe[n_elements(bidxe)-1] = narr-1$
        else if tidx[i] ne 0 then bidxe[tidx[i]-1] = bidxe[tidx[i]]
  endif

  bidxs = [bidxs[sort(bidxs)]]
  bidxs = [bidxs[uniq(bidxs)]]
  bidxe = [bidxe[sort(bidxe)]]
  bidxe = [bidxe[uniq(bidxe)]]

  ; refill IDXs
  for i=0,n_elements(bidxs)-1 do begin
    nbidx = (i eq 0)?[(lindgen(bidxe[i]-bidxs[i]+1)+bidxs[i])]:[nbidx,[(lindgen(bidxe[i]-bidxs[i]+1.)+bidxs[i])]]  
    gidx = (i eq 0)?[lindgen(bidxs[i]>1)]:[gidx,[(lindgen(bidxs[i]-bidxe[i-1]-1)+bidxe[i-1]+1)]] 
    gidxs = (i eq 0)?[0]:[gidxs,bidxe[i-1]+1]
    gidxe = (i eq 0)?[bidxs[i]-1]:[gidxe,bidxs[i]-1]
  endfor

  gidxs = [gidxs,bidxe[i-1]+1]
  gidxe = [gidxe,narr-1]
  if max(gidxs) eq narr then begin
    gidxs = gidxs[0:(n_elements(gidxs)-2)]
    gidxe = gidxe[0:(n_elements(gidxe)-2)]
  endif else gidx = [gidx,(lindgen(narr-bidxe[i-1]-1)+bidxe[i-1]+1)]

  if gidxe[0] eq -1 then begin
    gidx  = gidx[1:*]
    gidxs = gidxs[1:*]
    gidxe = gidxe[1:*]
  endif
  ; the following sorting is probably unneccesary
  gidxs = [gidxs[sort(gidxs)]]
  gidxs = [gidxs[uniq(gidxs)]]
  gidxe = [gidxe[sort(gidxe)]]
  gidxe = [gidxe[uniq(gidxe)]]

  if narr ne n_elements(gidx)+n_elements(nbidx) then stop,'We have too many indices, stopping here...'
  if n_elements(bidxs) ne n_elements(bidxe) then stop,'We have an index problem, stopping here...'
  if n_elements(gidxs) ne n_elements(gidxe) then stop,'We have an index problem, stopping here...'

  return,{GIDX:GIDX, GIDXS:GIDXS, GIDXE:GIDXE, BIDX:NBIDX, BIDXS:BIDXS, BIDXE:BIDXE}
End

Pro unwrap_disp_juwe, star, miditime, resultpath, timeorig, pdelorig, gdelorig, statorig, lamorig, $
	hak=hak,ti=ti,box=box,xr=xr,sav=sav

  time = timeorig
  pdel = pdelorig
  gdel = gdelorig
  stat = statorig
  lam = lamorig

; HAK plots a lot of debug windows
; SAV saves the final result in FILE.SAV
; TI defines a FSU time range to work on (default [0,0] -> full array)
; XR is to zoom the plots, we still calculate on TI (default [0,0] -> full array)
; BOX gives the averaging BOX in number of FSU frames (default is currently 100)
;  maybe ideally BOX should be correlated to the group delay tracking speed
;
; the data set needs to be loaded semi-manually below

  if keyword_set(hak) then hakk=1 else hakk=0
;;; configuration block ;;;

  ; set the average phase box, looking at the data, 100 [ms] seems apropriate
  ; but this is a critical value to make the phase-unwrapping 'sees' jumps
  if not keyword_set(box) then box = 100
  ; make box odd for the median calculation
  if (box mod 2) ne 1 then box=box+1
  
  if keyword_set(ti) then ti = double(ti) else ti = [0d,0d] 
  if keyword_set(xr) then xr = double(xr) else xr = [0d,0d]


  ;if (not keyword_set(ti)) or (total(ti eq [0,0]) eq 2) then ti=mima(time)
  ;tidx = where((time ge ti[0]) and (time le ti[1]))

  print,"Using BOX=",string(box),' FSU frames'
;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ntim = n_elements(time)

  pdel = pdel/(2.*!dpi)*lam
  ompdel = round(median(pdel)/lam)
  if ompdel ne 0  then stop,'Are you sure that you want to reduce data, where the median phase-delay is more than a wavelength away from 0? I am not, stopping here...'

  ; reject some datapoints and define IDXSTR

  ; here you might want to try 'lt 5' as well, 
  ; 7: fringe tracking, 5: fringe SNR below tracking, but still reasonable

  ;in case of perfect fringe track (state always 7), artificially insert an OPDC state 5
  idx7 = where(stat lt 7)
  if (idx7[0] eq -1) then stat[0] = 5

  bidx = where(stat lt 7)
  if bidx[0] gt -1 then begin
    idxstr = defgbidx(ntim,bidx,box) 
    pdel[idxstr.bidx] = !values.d_nan
    gdel[idxstr.bidx] = !values.d_nan
    bidx = 1 ; now used as flag, the indeces are in the IDXSTR
  endif else begin
    bidx = 0
  endelse


;   if hakk then begin
;     window,xs=1200,ys=1000
;     !p.multi=[0,1,4]
;   endif else begin
;     window,xs=1200,ys=250
;     !p.multi=[0,0,0]
;   endelse
  if hakk then begin

;     ;Calculate the aspect ratio of display window.
;     aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE 
;     xsize = 8.0
;       ysize = xsize * aspectRatio
;       IF ysize GT 10.5 THEN BEGIN
; 	ysize = 10.5
; 	xsize = ysize / aspectRatio
;       ENDIF
;       ; Calculate the offsets, so the output window is not off the page.
;       xoffset = (8.5 - xsize) / 2.0
;       yoffset = (11.0 - ysize) / 2.0
; 
;     set_plot, 'ps'
;     device, isolatin=1
;     device, filename=resultpath+star+'_'+miditime+'_unwrap.ps', /color,XSIZE=20, YSIZE=17, XOffset=xoffset, YOffset=yoffset

    window,xs=1200,ys=1000, title=star+' '+miditime
    !p.multi=[0,1,4]

  endif

  if hakk then begin
    micron = string(181B)
    ;!p.font=0
    ;!p.thick=4
    !x.thick=2
    !y.thick=2

    plot,time/1.d6,pdel*1.d6,/xst,/yst,xr=xr,yr=5.*[-1,1],psym=3,tit='orig. PD (black) / intra-data-unwrapped (blue) / individual-pixel-unwrapped (red)', charsize=3.4, background=fsc_color('white'), color=fsc_color('black'), $
    xtitle='Time [s]', ytitle='Delay ['+micron+'m]', yminor=2, xticklen=0.05, yticklen=0.01
    oplot,!x.crange,0*[1,1], color=fsc_color('black')
    oplot,!x.crange,-lam*1.d6*[1,1], color=fsc_color('black')
    oplot,!x.crange,lam*1.d6*[1,1], color=fsc_color('black')
  endif

  ; let's unwrap
  ; the strategy is to 
  ; FIRST: intra-data unwrap: 
  ;        define and clean-up areas with half of the phases
  ;        around pi and -pi, which belong together
  ;
  ; SECOND: unwrap individual pixel, 
  ;        where they phase-jump around the local median
  ;                       this step implies that data with phase-rms 
  ;                       gt lam/2 might be overcorrected.   
  ;                       but you probably couldn't phase track on such data 
  ;                       anyway...
  ; 
  ; THIRD: inter-pixel: 
  ; 	   unwrap where the local median jumps by more than lam/2
  ; 
  ; FOURTH: inter-block: unwrap if from one gidx-block to the next 
  ;                     there is a phase jump, 
  ;                     this jump is often NOT detected by the inter-pixel 
  ;                     unwrap, if between two gidx-blocks, the fringe track 
  ;                     was lost
  ;                     and the first retracking pixel are correctly measured 
  ;                     off the nominal position
  ;  

  ;;;;; FIRST: intra-data phase-unwrap ;;;;;;;
  ; identify areas with intra-data phase-wrapped datapoints 
  ; by comparing MEAN and MEDIAN 
  if bidx then begin
    ; estimate MEAN
    mpdel = replicate(!values.d_nan,ntim)
    ; estimate MEDIAN
    spdel = replicate(!values.d_nan,ntim)
    for i=0,n_elements(idxstr.gidxs)-1 do begin
      mpdel[idxstr.gidxs[i]:idxstr.gidxe[i]] = edge_med(pdel[idxstr.gidxs[i]:idxstr.gidxe[i]],box)
      spdel[idxstr.gidxs[i]:idxstr.gidxe[i]] = edge_smo(pdel[idxstr.gidxs[i]:idxstr.gidxe[i]],box)
    endfor
  endif else begin
    mpdel = edge_med(pdel,box)
    spdel = edge_smo(pdel,box)
  endelse

  dpdel = mpdel-spdel
  dplim = 2.5*stddev(dpdel,/nan)
  ;dplim = 7.0*stddev(dpdel,/nan)	;AMUELLER
  ;dplim = 5.d-8	;AMUELLER, experimentally derived value

  ;AMUELLER
  fixdplim = (!DPI/4.d0)*lam/(2.d0*!DPI)
  if (dplim lt fixdplim) then dplim = fixdplim

  pwidx = where(abs(dpdel) gt dplim)


  if hakk then begin
    oplot,time/1.d6,mpdel*1.d6,color=fsc_color('cyan')
    oplot,time/1.d6,spdel*1.d6,color=fsc_color('green')
    oplot,time/1.d6,dpdel*1.d6,col=fsc_color('brown'),thick=2

    oplot,!x.crange,replicate(dplim,2)*1.d6,lin=1, color=fsc_color('black')
    oplot,!x.crange,replicate(-dplim,2)*1.d6,lin=1, color=fsc_color('black')
  endif 



  if pwidx[0] ge 0 then begin
    print,'#pix found with intra-data phase wrap: ',n_elements(pwidx)
    ; increase pwidx to contiuous pieces of BOX length
    pwidxstr = (defgbidx(ntim,pwidx,box)).bidx
    ; adapt pwidxstr to the current idxstr.gidx's
    pwidxstr = pwidxstr[where(finite(pdel[pwidxstr]) eq 1)]
    ; get starts and ends
    pwidxstr = defgbidx(ntim,pwidxstr,0)

    ; let's wrap the pixel with phase gt 0 towards the pixel with phase le 0
    midx = where(pdel[pwidxstr.bidx] gt 0)
    pdel[pwidxstr.bidx[midx]] = (pdel[pwidxstr.bidx[midx]]-lam)

    if hakk then begin
      ; plot some areas of action:
      ; originally found pixel where MEDIAN-MEAN difference exceed DPLIM
      oplot,time[pwidx]/1.d6,replicate(-1.9e-6,n_elements(pwidx))*1.d6,psym=3, color=fsc_color('black')
      ; pwidx increased to BOX length
      oplot,time[pwidxstr.bidx]/1.d6,replicate(-2.0e-6,n_elements(pwidxstr.bidx))*1.d6,psym=3, color=fsc_color('black')
      ; pixel where we unwrap
      oplot,time[pwidxstr.bidx[midx]]/1.d6,replicate(-2.1e-6,n_elements(pwidxstr.bidx[midx]))*1.d6,psym=3,col=fsc_color('brown')

      ; areas of NAN
      if bidx then oplot,time[idxstr.bidx]/1.d6,replicate(2e-6,n_elements(idxstr.bidx))*1.d6,psym=3,col=fsc_color('red') $
      else oplot,time/1.d6,replicate(2e-6,ntim)*1.d6,psym=3,col=fsc_color('red')
    endif    
  endif else begin
    print,'no pix found with intra-data phase wrap.'
  endelse


  if hakk then  oplot,time/1.d6,pdel*1.d6,psym=4,color=fsc_color('blue')

  ;;;;; SECOND: pixel-unwrap  ;;;;;
  mpdel = replicate(!values.d_nan,ntim)
  dmpdel = replicate(!values.d_nan,ntim)

  ; estimate median of current PDEL
  if bidx then $
    for i=0,n_elements(idxstr.gidxs)-1 do $
      mpdel[idxstr.gidxs[i]:idxstr.gidxe[i]] = edge_med(pdel[idxstr.gidxs[i]:idxstr.gidxe[i]],box) $
        else mpdel = edge_med(pdel,box)


  if bidx then begin
    dmpdel[idxstr.gidx] $
    = floor((pdel[idxstr.gidx]-mpdel[idxstr.gidx])/lam + 0.5)
    pdel[idxstr.gidx]=pdel[idxstr.gidx]-dmpdel[idxstr.gidx]*lam
  endif else begin
    dmpdel = floor((pdel-mpdel)/lam + 0.5) 
    pdel = pdel-dmpdel*lam
  endelse
  if hakk then oplot,time/1.d6,pdel*1.d6,psym=4,color=fsc_color('red')



  ;;;;; THIRD: median-unwrap ;;;;;
  ; this sort of assumes that, after intra-data and individual-pixel unwrap, 
  ; the areas of
  ; stable phase without the need for median-unwrapping are larger than BOX
  ; so that we have a meaning full median per box


  ; estimate median of current PDEL
  ; first unwrap WITHIN gidx parts, where we intra-data unwrapped
  ; that is, between idxstr.gidx and pwidxstr.bidx
 
  ; ALL good data points are in idxstr.gidx
  ; pwidxstr.bidx describes the subset where we intra-data unwrapped

;   if bidx then $
;     mwidxstr = defgbidx(ntim,setdifference(idxstr.gidx,pwidxstr.bidx),0) $
;   else mwidxstr = defgbidx(ntim,setdifference(lindgen(ntim),pwidxstr.bidx),0)
;   ; now mwidxstr.bidx contains the GOOD data where we did NOT intra-data unwrap
  if bidx then begin
    if n_elements(pwidxstr) gt 0 then $ 
      mwidxstr = defgbidx(ntim,setdifference(idxstr.gidx,pwidxstr.bidx),0) $
    else mwidxstr = defgbidx(ntim,idxstr.gidx,0)
  endif else begin
    if n_elements(pwidxstr) gt 0 then $ 
     mwidxstr = defgbidx(ntim,setdifference(lindgen(ntim),pwidxstr.bidx),0) $
     else mwidxstr = defgbidx(ntim,lindgen(ntim),0) 
  endelse
  ; now mwidxstr.bidx contains the GOOD data where we did NOT intra-data unwrap
;--------------------

  mpdel = replicate(!values.d_nan,ntim)
  ; calc the median on continuos data good strips
  for i=0,n_elements(mwidxstr.bidxs)-1 do $
      mpdel[mwidxstr.bidxs[i]:mwidxstr.bidxe[i]] = edge_med(pdel[mwidxstr.bidxs[i]:mwidxstr.bidxe[i]],box) 
   if n_elements(pwidxstr) gt 0 then begin 
     for i=0,n_elements(pwidxstr.bidxs)-1 do $
      mpdel[pwidxstr.bidxs[i]:pwidxstr.bidxe[i]] = edge_med(pdel[pwidxstr.bidxs[i]:pwidxstr.bidxe[i]],box) 
    endif

  if hakk then begin
    oplot,time/1.d6,mpdel*1.d6, color=fsc_color('black')
    if n_elements(pwidxstr) gt 0 then $ oplot,time[pwidxstr.bidx]/1.d6,replicate(-2.4e-6,n_elements(pwidxstr.bidx))*1.d6,psym=3,col=fsc_color('red')
    oplot,time[mwidxstr.bidx]/1.d6,replicate(-2.4e-6,n_elements(mwidxstr.bidx))*1.d6,psym=3,col=fsc_color('black')
  endif

  ; now estimate locations of potential jumps due to intra-data unwrap 
  ; WITHIN gidx strips

if n_elements(pwidxstr) gt 0 then begin $

  mwidx1 = setintersection(idxstr.gidx,(pwidxstr.bidxs-1))
  mwidx2 = setintersection((idxstr.gidx-1),pwidxstr.bidxe)
  ; test
  if total(setintersection(mwidx1,mwidxstr.bidxe) ne mwidx1) gt 0 then stop,'We do not understand the indexing yet...'
  if total(setintersection(mwidx2+1,mwidxstr.bidxs) ne (mwidx2+1)) gt 0 then stop,'We do not understand the indexing yet...'

  mwidx = [mwidx1,mwidx2]
  mwidx = [mwidx[sort(mwidx)]]
  mwidx = [mwidx[uniq(mwidx)]]
  ; now mwidx points to data AFTER the jump
  mwidx = mwidx+1
  
  if hakk then for i=0,n_elements(mwidx)-1 do oplot,replicate(time[mwidx[i]],2)/1.d6,!y.crange

  if bidx then begin
    dmpdel[idxstr.gidx] $
    = shift(floor(-ts_diff(mpdel[idxstr.gidx],1,/dou)/lam + 0.5),1)
    ; only unwrap within good areas (and not between), where we intra-dataunwrapped
    dmpdel[setdifference(idxstr.gidx,mwidx)] = 0
    pdel[idxstr.gidx]=pdel[idxstr.gidx]-total(dmpdel[idxstr.gidx],/cum,/dou,/nan)*lam
  endif else begin
    dmpdel = shift(floor(-ts_diff(mpdel,1,/dou)/lam + 0.5),1)
    dmpdel[setdifference(lindgen(ntim),mwidx)] = 0
    pdel = pdel-total(dmpdel,/cum,/dou,/nan)*lam
  endelse

endif

  if hakk then begin

    ;AMUELLER
    ;there could be NAN values in pdel, plot will not work, create dummy
    idxdum = where(finite(pdel) eq 0)
    pdeldum = pdel
    if (idxdum[0] ne -1) then pdeldum[idxdum] = 0.

    plot,time/1.d6,pdel*1.d6,/xst,/yst,xr=xr,yr=[-1,1]*max(abs(pdeldum*1.d6)),psym=3,tit='unwrapped PD after intra-data-block unwrap (black)/ after gidx-block unwrap (red) / raw GD (blue)', charsize=3.4, background=fsc_color('white'), color=fsc_color('black'), xtitle='Time [s]', ytitle='Delay ['+micron+'m]', $
    yminor=2, xticklen=0.05, yticklen=0.01
    if bidx then begin
      tmp = !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0])
      oplot,time[idxstr.gidx]/1.d6,dmpdel[idxstr.gidx]*0.5+tmp,psym=2,col=fsc_color('brown') 
    endif else begin
       tmp = !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0])
       oplot,time/1.d6,dmpdel*0.5+tmp,psym=2,col=fsc_color('brown')
    endelse
    oplot,time/1.d6,gdel*1.d6,color=fsc_color('blue'),psym=3

  endif

  ; FOURTH: gidx-block unwrap
  ; we pull the median down to zero since the overall median pdel 
  ; is zero (see the ompdel-check at the beginning of this code)
  ; and since otherwise for the PDEL alone, we do not have a better reference
  ; in case of fringe drifts, the rest needs to be done 
  ; by gidx-block disp-unwrapping

  ; estimate median of current PDEL
  if bidx and (n_elements(idxstr.gidxs) ge 2) then $
    for i = 1,n_elements(idxstr.gidxs)-1 do $
      pdel[idxstr.gidxs[i]:idxstr.gidxe[i]]  $
        = pdel[idxstr.gidxs[i]:idxstr.gidxe[i]] $
           - (floor((median(pdel[idxstr.gidxs[i]:idxstr.gidxe[i]]) $
               - ompdel)/lam + 0.5)*lam)


  if hakk then oplot,time/1.d6,pdel*1.d6,psym=3,col=fsc_color('red')



  ; phase unwrap to keep the DISP centered on median(DISP)
  ; here we give up to see DISP jumps / drifts beyond lam/2
  ; essentially we unwrap each pixel if the box-median is off-centered
  ; I am not sure if that is the smartest idea, but we'll see...

  ; calculate median gdel
  mgdel = replicate(!values.d_nan,ntim)

  if bidx then $
    for i=0,n_elements(idxstr.gidxs)-1 do $
      mgdel[idxstr.gidxs[i]:idxstr.gidxe[i]] = edge_med(gdel[idxstr.gidxs[i]:idxstr.gidxe[i]],box) $
  else mgdel = edge_med(gdel,box)

  mgdisp = pdel - mgdel
;  mmgdisp = median(mgdisp)
  mmgdisp = 0
  disp = pdel - gdel

  ; plot what we have
  if hakk then begin
    plot,time/1.d6,disp*1.d6,/xst,/yst,xr=xr,yr=15.*[-1,1],psym=3,tit='orig. DISP (black; median in cyan) / median-pix-unwrapped DISP (red; and the resp. PD in blue) ', charsize=3.4, background=fsc_color('white'), color=fsc_color('black'), xtitle='Time [s]', ytitle='Delay ['+micron+'m]', yminor=2, xticklen=0.05, yticklen=0.01
    oplot,time/1.d6,mgdisp*1.d6,col=fsc_color('cyan'),psym=2,symsize=0.5
    oplot,!x.crange,mmgdisp*[1,1]*1.d6, color=fsc_color('black')
    oplot,!x.crange,mmgdisp-0.66*lam*[1,1]*1.d6,lin=1, color=fsc_color('black')
    oplot,!x.crange,mmgdisp+0.66*lam*[1,1]*1.d6,line=1, color=fsc_color('black')
  endif

  ; define stable, continuous DISP areas within the gidx blocks
  ; we define an area as STABLE if the median does not drift away 
  ; by more than lam*0.66 on scales of BOX
  ; this 0.66*lam seems to work, something like a 3*sigm definition 
  ; would be better, I also tried unbiased lam/2 but then there is 
  ; a few not realistic jumps
  ; also when disp rises more slowly then BOX, we wont detect this,
  ; but his should not happen, if we group-delay track at 10 Hz or so.
  ; Probably BOX needs to be adapted to the group delay tracking speeed.

  sumdis = replicate(!values.d_nan,ntim)

  if bidx then begin
    for i = 0, n_elements(idxstr.gidxs)-1 do begin
      dmgdisp = -ts_diff(mgdisp[idxstr.gidxs[i]:idxstr.gidxe[i]],1,/dou)
      ; SUM contains the integrated adjacent BOX mgdisp differences:
      sumdis[idxstr.gidxs[i]:idxstr.gidxe[i]] = edge_smo(dmgdisp,box)*double(box)
    endfor
  endif else begin
    dmgdisp = -ts_diff(mgdisp,1,/dou)
    sumdis = edge_smo(dmgdisp,box)*double(box)
  endelse

  if hakk then oplot,time/1.d6,sumdis*1.d6,col=fsc_color('red')
  jumps = where((abs(sumdis) gt 0.66*lam) and (finite(sumdis) eq 1))

;;;;;;;;;;;;;;;;
  if jumps[0] ge 0 then begin
    print,'#pix found with disp jumps: ',n_elements(jumps)
    ; increase jumps to contiuous pieces of BOX length
    juidxstr = (defgbidx(ntim,jumps,box)).bidx
    ; adapt juidxstr to the current idxstr.gidx's 
    juidxstr = juidxstr[where(finite(pdel[juidxstr]) eq 1)]
    ; get starts and ends
    juidxstr = defgbidx(ntim,juidxstr,0)
    jidx =1
    if hakk then begin
      ; plot some areas of action:
      ; originally found pixel where MEDIAN-MEAN difference exceed DPLIM
      oplot,time[jumps]/1.d6,replicate(-2e-6,n_elements(jumps))*1.d6,psym=3, color=fsc_color('black')
      ; jumps increased to BOX length
      oplot,time[juidxstr.bidx]/1.d6,replicate(-2.5e-6,n_elements(juidxstr.bidx))*1.d6,psym=3, color=fsc_color('black')
      ; plot areas of NAN
      if bidx then oplot,time[idxstr.bidx]/1.d6,replicate(2e-6,n_elements(idxstr.bidx))*1.d6,psym=3,col=fsc_color('red') $
      else oplot,time/1.d6,replicate(2e-6,ntim)*1.d6,psym=3,col=fsc_color('red')
    endif    
  endif else begin
    print,'no pix found with intra-data phase wrap.'
    jidx = 0
  endelse

  ; now exclude jumps from gidx (if any)

  if jidx then begin
    nbidx = (bidx eq 1)?[idxstr.bidx,juidxstr.bidx]:juidxstr.bidx
    nbidx = [nbidx[sort(nbidx)]]
    nbidx = [nbidx[uniq(nbidx)]]
    nidxstr = defgbidx(ntim,nbidx,box)
    bidx = 1
  endif else nidxstr = idxstr
  if hakk then if bidx then oplot,time[nidxstr.bidx]/1.d6,replicate(-3e-6,n_elements(nidxstr.bidx))*1.d6,psym=3,col=fsc_color('red')

  ; pull each gidx-group of DISP down to previous DISP level
  if hakk then begin
    for i = 1,n_elements(nidxstr.gidxs)-1 do begin
      oplot,[mean(time[nidxstr.gidxs[i]:nidxstr.gidxe[i]])]/1.d6,[median(disp[nidxstr.gidxs[i]:nidxstr.gidxe[i]])]*1.d6,psym=2,syms=3, color=fsc_color('black')
      print,mean(time[nidxstr.gidxs[i]:nidxstr.gidxe[i]]),median(disp[nidxstr.gidxs[i]:nidxstr.gidxe[i]])-median(disp[nidxstr.gidxs[i-1]:nidxstr.gidxe[i-1]]),(floor((median(disp[nidxstr.gidxs[i]:nidxstr.gidxe[i]]) $
             - median(disp[nidxstr.gidxs[i-1]:nidxstr.gidxe[i-1]]))/lam + 0.5)*lam)
    endfor
  endif

  if (n_elements(nidxstr.gidxs) ge 2) then begin
    for i = 1,n_elements(nidxstr.gidxs)-1 do begin
      pdel_corr = (floor((median(disp[nidxstr.gidxs[i]:nidxstr.gidxe[i]]) $
           - median(disp[nidxstr.gidxs[i-1]:nidxstr.gidxe[i-1]]))/lam + 0.5)*lam)
      disp[nidxstr.gidxs[i]:nidxstr.gidxe[i]] = $
         disp[nidxstr.gidxs[i]:nidxstr.gidxe[i]] - pdel_corr
      pdel[nidxstr.gidxs[i]:nidxstr.gidxe[i]] = $
         pdel[nidxstr.gidxs[i]:nidxstr.gidxe[i]] - pdel_corr
    endfor
  endif
  if hakk then oplot,time[nidxstr.gidx]/1.d6,disp[nidxstr.gidx]*1.d6,color=fsc_color('blue'),psym=3

  disp[nidxstr.bidx] = !values.d_nan
  pdel[nidxstr.bidx] = !values.d_nan
  gdel[nidxstr.bidx] = !values.d_nan
  ; final median
  if bidx then begin
    mdisp = replicate(!values.d_nan,ntim)
    for i=0,n_elements(nidxstr.gidxs)-1 do $
      mdisp[nidxstr.gidxs[i]:nidxstr.gidxe[i]] = edge_med(disp[nidxstr.gidxs[i]:nidxstr.gidxe[i]],box)
  endif else mdisp = edge_med(disp,box)

  ; plot the final results
  if hakk then begin	;AMUELLER
    yr = 1.d6*(mima([disp,pdel,gdel,mdisp]))
    plot,time/1.d6,disp*1.d6,/xst,/yst,xr=xr,yr=yr,psym=3,tit='final unwrapped DISP (black) / box-median DISP (brown)/ final unwrapped PD (red) / GD (blue)',ytit='Delay ['+micron+'m]',xtit='Time [s]', charsize=3.4, background=fsc_color('white'), color=fsc_color('black'), yminor=2, xticklen=0.05, yticklen=0.01
    oplot,time/1.d6,pdel*1.d6,color=fsc_color('red'),psym=3
    oplot,time/1.d6,gdel*1.d6,color=fsc_color('blue'),psym=3
    oplot,time/1.d6,mdisp*1.d6,color=fsc_color('brown'),psym=2,syms=0.5
    ; for orientation overplot few lines
    oplot,!x.crange,replicate(0,2)*1.d6, color=fsc_color('black')
    oplot,!x.crange,replicate(-lam,2)*1.d6, color=fsc_color('black')
    oplot,!x.crange,replicate(lam,2)*1.d6, color=fsc_color('black')
  endif


  if keyword_set(sav) then begin
    fn = resultpath+star+'_'+miditime+'_unwrap.sav'
    print,'Saving final data in '+fn
    save,filename=fn,time,nidxstr,disp,pdel,gdel,box,mdisp
  endif

  image = tvread(filename=resultpath+star+'_'+miditime+'_unwrap', type='tiff', /nodialog)

;   device,/close
;   set_plot,'x'

  !p.multi = [0,1,0]
;stop
End
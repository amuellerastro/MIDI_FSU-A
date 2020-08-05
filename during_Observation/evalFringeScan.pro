@get_eso_keyword.pro
@mpfit.pro
@mpfitfun.pro

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

pro evalFringeScan

path = '/home/amueller/pott/data/'
;path = '/media/disk/MIDI_FSUA/COMM16/FSUAdata/'


;don't use ls command but look for the index
file = file_search(path+'*.fits')
pos1 = strpos(file, '.fits')

num1 = strarr(n_elements(pos1))
for i=0,n_elements(pos1)-1 do num1[i] = strmid(file[i], pos1[i]-8, 8)
idx = bsort(num1, /reverse)
file = file[idx[0]]


print, 'using file '+file

;for obs_gen files
;file = 'PACMAN_OBS_GENERIC243_0011.fits'

;for gen_rec files
;file = 'PACMAN_GEN_RECORD_332_0024.fits'
;file = 'PACMAN_OBJ_SCAN_247_0003.fits'

;************************
;*  DO NOT EDIT BELOW   *
;************************

path0 = 'SimpleScan_results/'
file_mkdir, path0


;extract file name
pos1 = strpos(file, 'PAC')
name = strmid(file, strpos(file, 'PAC'), 26)

; file = path+file

;read in FSU-A data

  ;get DL sign from header
  dum = mrdfits(file, 0, hdr)
  opdsign = get_eso_keyword(hdr, 'HIERARCH ESO ISS DEL PRMOPDSIGN')
  opdsign = double(opdsign)
  ;read, 'DL sign: ', opdsign

  ;st1 = mrdfits(file, 4)
  st1 = mrdfits(file, 7)
  t1 = st1.time
  wA = st1.data1[0,*]
  wC = st1.data3[0,*]
  wB = st1.data2[0,*]
  wD = st1.data4[0,*]

  pdsnr = st1.pdsnr

  ;st2 = mrdfits(file, 8)
  st2 = mrdfits(file, 11)
  t2 = st2.time	;double sampling
  fuoff_orig = st2.fuoffset

  index_fsu_time = closest(t2,t1)
  fuoff = fuoff_orig[index_fsu_time]

;   dum = max(fuoff, idxmax)
;   fuoff = fuoff[0:idxmax]
;   wA = wA[0:idxmax]
;   wC = wC[0:idxmax]
;   wB = wB[0:idxmax]
;   wD = wD[0:idxmax]

;compute difference between white light pixel of channel A and C

  diffAC = wA-wC
  diffBD = wB-wD



;find the maximum SNR within the scan
  mmy = max(pdsnr, mmx)

;plot full scan range
  window, 0, xs=840, ys=500, title=file
  !p.multi=[0,1,2]
  plot, fuoff, diffAC, xst=1, yst=1, charsize=1.4, xtitle='Full Offset [m]', ytitle='KWA-KWC'
  oplot, fuoff[mmx-50:mmx+50], diffAC[mmx-50:mmx+50], color=fsc_color('red')
  legend, ['white light A-C'], /left, box=0, margin=0, charsize=1.5
  plot, fuoff, diffBD, xst=1, yst=1, charsize=1.4, xtitle='Full Offset [m]', ytitle='KWB-KWD'
  oplot, fuoff[mmx-50:mmx+50], diffBD[mmx-50:mmx+50], color=fsc_color('red')
  legend, ['white light B-D'], /left, box=0, margin=0, charsize=1.5
  !p.multi=[0,1,0]

;define a region of +/- fw around the max snr

  fw = 0.000010d0

  ran = where(fuoff ge fuoff[mmx]-10.d0*fw and fuoff le fuoff[mmx]+10.d0*fw)

;fit this maximum snr value as bootstrapping for futher processing

  ;--- define start parameters for fit ---
      ;c = [fuoff[mmx], 0.05d-4, 5.d6, 0.d0, 0.1d0, fuoff[mmx], 0.1d-4, mmy, 1.d0];
      c = [fuoff[mmx], 0.04d-4, 5.006385d6, -19.d0, 0.21d0, fuoff[mmx], -1.33d-6, mmy, 0.96d0]

  ;--- do least-squares fit and save param. in cc --- 
      dumerr = dblarr(n_elements(ran))
      dumerr[*] = 1.
      maxiter = 1000
      cc = mpfitfun('snrfun', fuoff[ran], pdsnr[ran], dumerr, c, weights=dumerr, $
	   maxiter=maxiter, niter=niter, status=status, bestnorm=bestnorm, $
	   perror=perror, covar=covar, yfit=yfit, dof=dof, /quiet)
      ;cc = snrfit(c,obs.foff(ran)',obs.snr(ran));
      cci = cc;


;compute OPDC SNR thresholds and ZPD offset to be applied

  snrmax = max(yfit, idx)
  fuoff_tmp = fuoff[ran]
  snrpos = fuoff_tmp[idx]
  zpdoff = snrpos*OPDsign

  openl = cc[8]
  detl = snrfun(cc[0]-1.d0*cc[1],cc)
  clol = snrfun(cc[0]-2.d0*cc[1],cc)

;plot SNR and fit
  window, 1, xs=500, ys=500, title=file
  plot, fuoff[ran], pdsnr[ran], xtitle='Full Offset [m]', ytitle='SNR', xst=1, charsize=1.5
  oplot, fuoff[ran], yfit, color=fsc_color('red'), thick=4
  plots, [snrpos,snrpos], [!y.crange[0],!y.crange[1]], linestyle=1, color=fsc_color('blue'), thick=4
  plots, [!x.crange[0],!x.crange[1]], [openl, openl], linestyle=1, color=fsc_color('yellow'), thick=4
  plots, [!x.crange[0],!x.crange[1]], [clol, clol], linestyle=1, color=fsc_color('yellow'), thick=4
  plots, [!x.crange[0],!x.crange[1]], [detl, detl], linestyle=1, color=fsc_color('yellow'), thick=4
  xyouts, 0.21, 0.88, 'ZPD offset: '+string(zpdoff)+' m', /normal, charsize=1.5

;plot detected fringe

  ;--- define plot range max(snr)+/-2*fringewidth
  ran2 = where(fuoff ge fuoff[mmx]-5.d0*fw and fuoff le fuoff[mmx]+5.d0*fw)

  window, 2, xs=840, ys=500, title=file
  !p.multi=[0,2,1]
  plot, fuoff[ran2], diffAC[ran2], xst=1, charsize=1.5, xtitle='Full Offset [m]', ytitle='KWA-KWC [counts]'
  plots, [snrpos,snrpos], [!y.crange[0],!y.crange[1]], linestyle=1, color=fsc_color('blue'), thick=4
  plot, fuoff[ran2], diffBD[ran2], xst=1, charsize=1.5, xtitle='Full Offset [m]', ytitle='KWB-KWD [counts]'
  plots, [snrpos,snrpos], [!y.crange[0],!y.crange[1]], linestyle=1, color=fsc_color('blue'), thick=4
  !p.multi=[0,1,0]


print, ''
print, '********************************'
print, 'used file   : ', file
print, 'ZPD offset  : ', zpdoff, ' m'
print, 'Detec. Level: ', detl
print, 'Close Level : ', clol
print, 'Open Level  : ', openl
print, '********************************'
print, ''


;outpu text file

  openw, lun, path0+'SimpleScan_'+name+'.txt', width=1400, /get_lun
    printf, lun, name
    printf, lun, ''
    printf, lun, 'ZPD offset  : ', zpdoff, ' m'
    printf, lun, ''
    printf, lun, 'Detec. Level: ', detl
    printf, lun, 'Close Level : ', clol
    printf, lun, 'Open Level  : ', openl
  close, lun
  free_lun, lun



;output as a PS file

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

set_plot, 'ps'
device, filename=path0+'SimpleScan_'+name+'.ps',/color,XSIZE=22, YSIZE=12, XOffset=xoffset, YOffset=yoffset
device, isolatin=1
;   !P.Font=0
  !p.thick=4
  !x.thick=3
  !y.thick=3

  !p.multi=[0,1,2]
  plot, congrid(fuoff, n_elements(fuoff)/5.d0)*1000., congrid(transpose(diffAC), n_elements(fuoff)/5.d0), xst=1, charsize=1.5, ytitle='KWA-KWC [ADU]', yminor=2, yticks=4, xticklen=0.07, yst=1, $
  xtickformat="(A1)", pos=[0.14,0.565,0.97, 0.97]
  ;oplot, fuoff[mmx-50:mmx+50], diffAC[mmx-50:mmx+50], color=fsc_color('red')
;  legend, ['white light A-C'], /left, box=0, margin=0, charsize=1.3
  plot, congrid(fuoff, n_elements(fuoff)/5.)*1000., congrid(transpose(diffBD), n_elements(fuoff)/5.), xst=1, charsize=1.5, xtitle='Full Offset [mm]', ytitle='KWB-KWD [ADU]', yminor=2, yticks=4, xticklen=0.07, yst=1, pos=[0.14,0.145,0.97, 0.56]
;   oplot, fuoff[mmx-50:mmx+50], diffBD[mmx-50:mmx+50], color=fsc_color('red')
;  legend, ['white light B-D'], /left, box=0, margin=0, charsize=1.3


  !p.multi=[0,2,1]
  plot, fuoff[ran2]*1000., diffAC[ran2], xst=1, charsize=1.5, xtitle='Full Offset [mm]', ytitle='KWA-KWC [ADU]', xminor=2, xticks=2, yticklen=0.06, yminor=2, xticklen=0.05, yst=1, $
  pos=[0.14,0.14,0.43, 0.97]
  plots, [snrpos,snrpos]*1000., [!y.crange[0],!y.crange[1]], linestyle=1, color=fsc_color('blue'), thick=4
  plot, fuoff[ran2]*1000., diffBD[ran2], xst=1, charsize=1.5, xtitle='Full Offset [mm]', ytitle='KWB-KWD [ADU]', xminor=2, xticks=2, yticklen=0.06, yminor=2, xticklen=0.05, yst=1, $
  pos=[0.63,0.14,0.92, 0.97]
  plots, [snrpos,snrpos]*1000., [!y.crange[0],!y.crange[1]], linestyle=1, color=fsc_color('blue'), thick=4


  !p.multi=[0,1,0]
  plot, fuoff[ran]*1000., pdsnr[ran], xtitle='Full Offset [mm]', ytitle='SNR', xst=1, charsize=1.5, yminor=2, pos=[0.1,0.14,0.95, 0.97], xminor=2, xticklen=0.04
  oplot, fuoff[ran]*1000., yfit, color=fsc_color('red'), thick=4
  plots, [snrpos,snrpos]*1000., [!y.crange[0],!y.crange[1]], linestyle=1, color=fsc_color('blue'), thick=4
  plots, [!x.crange[0],!x.crange[1]], [openl, openl], linestyle=1, color=fsc_color('green'), thick=4
  plots, [!x.crange[0],!x.crange[1]], [clol, clol], linestyle=1, color=fsc_color('green'), thick=4
  plots, [!x.crange[0],!x.crange[1]], [detl, detl], linestyle=1, color=fsc_color('green'), thick=4
  ;xyouts, 0.21, 0.88, 'ZPD offset: '+string(zpdoff)+' m', /normal, charsize=1.5


device,/close
set_plot,'x'

!p.multi=[0,1,0]
; !P.Font=0
!p.thick=1
!x.thick=1
!y.thick=1

;spawn, 'gv '+path0+'SimpleScan_'+name+'.ps'

stop
end

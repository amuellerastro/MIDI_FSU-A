@closest2.pro
@mean.pro
@moment.pro
@stddev.pro
@ps_background.pro
@cgcolorfill.pro
@imsl_polyregress.pro
@filepath.pro
@path_sep.pro
@imsl_long.pro
@imsl_size.pro
@imsl_polypredict.pro
@imsl_n_elements.pro
@showsym.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@caldat.pro
@mpfitfun.pro
@mpfit.pro
@fsc_color.pro
@tvread.pro
@ploterror.pro
@setdefaultvalue.pro
@cgplot.pro
@setdecomposedstate.pro
@decomposedcolor.pro
@cgdefaultcolor.pro
@colorsareidentical.pro
@cgcolor.pro
@cgsnapshot.pro
@getdecomposedstate.pro
@oploterror.pro
@tag_exist.pro
@cgquery.pro
@cgdefcharsize.pro
@symcat.pro
@legend.pro

function func_poly_deg1, x, p
  fit = p[0]+p[1]*x
  return, fit
end
function func_poly_deg2, x, p
  fit = p[0]+p[1]*x+p[2]*x^2.
  return, fit
end
function func_poly_deg3, x, p
  fit = p[0]+p[1]*x+p[2]*x^2.+p[3]*x^3.
  return, fit
end
function func_poly_deg4, x, p
  fit = p[0]+p[1]*x+p[2]*x^2.+p[3]*x^3.+p[4]*x^4.
  return, fit
end
function func_poly_deg5, x, p
  fit = p[0]+p[1]*x+p[2]*x^2.+p[3]*x^3.+p[4]*x^4.+p[5]*x^5.
  return, fit
end
function func_poly_deg6, x, p
  fit = p[0]+p[1]*x+p[2]*x^2.+p[3]*x^3.+p[4]*x^4.+p[5]*x^5.+p[6]*x^6.
  return, fit
end

function tffunc_interp3, ncal, nsci, hour_sci, hour_cal, T2_AC_cal, T2_BD_cal, T2_AC_err_cal, T2_BD_err_cal, npix

  ;fit polynom through all TFs
  deg = 3
  polfit_AC = dblarr(npix,ncal) & polfit_BD = dblarr(npix,ncal)
  fitparams_AC = dblarr(npix,deg+1) & fitparams_BD = dblarr(npix,deg+1)
  interval_AC = dblarr(npix,2,ncal) & interval_BD = dblarr(npix,2,ncal)
  T2_AC_err_sci = dblarr(npix,nsci) & T2_BD_err_sci = dblarr(npix,nsci)

  for i=0,npix-1 do begin

    ;AC polynomial fit and confidence intervals

    param = imsl_polyregress(hour_cal, reform(T2_AC_cal[i,*]), deg, weight=1.d0/reform(T2_AC_err_cal[i,*])^2., residual=residual, predict_info=infoAC, /double)
    fitparams_AC[i,*] = param
    polfit_AC[i,*] = func_poly_deg3(hour_cal, param)

    ;confidence intervals for T^2 Cal AC
    result = imsl_polypredict(infoAC, hour_cal, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, y=reform(T2_AC_cal[i,*]), /double)
    interval_AC[i,*,*] = popmean

    ;confidence intervals for T^2 Sci AC
    result = imsl_polypredict(infoAC, hour_sci, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, /double)
    T2_AC_err_sci[i,*] = abs(reform(popmean[0,*])-reform(popmean[1,*]))/2.d0


    ;BD polynomial fit and confidence intervals
    param = imsl_polyregress(hour_cal, reform(T2_BD_cal[i,*]), deg, weight=1.d0/reform(T2_BD_err_cal[i,*])^2., residual=residual, predict_info=infoBD, /double)
    fitparams_BD[i,*] = param
    polfit_BD[i,*] = func_poly_deg3(hour_cal, param)

    ;confidence intervals for T^2 Cal BD
    result = imsl_polypredict(infoBD, hour_cal, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, y=reform(T2_AC_cal[i,*]), /double)
    interval_BD[i,*,*] = popmean

    ;confidence intervals for T^2 Sci BD
    result = imsl_polypredict(infoBD, hour_sci, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, /double)
    T2_BD_err_sci[i,*] = abs(reform(popmean[0,*])-reform(popmean[1,*]))/2.d0

  endfor

  ;------------------------------------------------------------------------------

  ;get TF for time of Sci observations

  T2_AC_sci = dblarr(npix, nsci) & T2_BD_sci = dblarr(npix, nsci)

  for i=0,npix-1 do begin

    for j=0,nsci-1 do begin

      T2_AC_sci[i,j] = func_poly_deg3(hour_sci[j], fitparams_AC[i,*])
      T2_BD_sci[i,j] = func_poly_deg3(hour_sci[j], fitparams_BD[i,*])

    endfor

  endfor

  st = {T2_AC_sci:T2_AC_sci, T2_AC_err_sci:T2_AC_err_sci, T2_BD_sci:T2_BD_sci, T2_BD_err_sci:T2_BD_err_sci, polfit_AC:polfit_AC, polfit_BD:polfit_BD, interval_AC:interval_AC, interval_BD:interval_BD};, polfit_err_AC:polfit_err_AC, polfit_err_BD:polfit_err_BD

  return, st

end

function tffunc_interp2, ncal, nsci, hour_sci, hour_cal, T2_AC_cal, T2_BD_cal, T2_AC_err_cal, T2_BD_err_cal, npix

  ;fit polynom through all TFs
  deg = 2
  polfit_AC = dblarr(npix,ncal) & polfit_BD = dblarr(npix,ncal)
  fitparams_AC = dblarr(npix,deg+1) & fitparams_BD = dblarr(npix,deg+1)
  interval_AC = dblarr(npix,2,ncal) & interval_BD = dblarr(npix,2,ncal)
  T2_AC_err_sci = dblarr(npix,nsci) & T2_BD_err_sci = dblarr(npix,nsci)

  for i=0,npix-1 do begin

    ;AC polynomial fit and confidence intervals

    param = imsl_polyregress(hour_cal, reform(T2_AC_cal[i,*]), deg, weight=1.d0/reform(T2_AC_err_cal[i,*])^2., residual=residual, predict_info=infoAC, /double)
    fitparams_AC[i,*] = param
    polfit_AC[i,*] = func_poly_deg2(hour_cal, param)

    ;confidence intervals for T^2 Cal AC
    result = imsl_polypredict(infoAC, hour_cal, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, y=reform(T2_AC_cal[i,*]), /double)
    interval_AC[i,*,*] = popmean

    ;confidence intervals for T^2 Sci AC
    result = imsl_polypredict(infoAC, hour_sci, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, /double)
    T2_AC_err_sci[i,*] = abs(reform(popmean[0,*])-reform(popmean[1,*]))/2.d0


    ;BD polynomial fit and confidence intervals
    param = imsl_polyregress(hour_cal, reform(T2_BD_cal[i,*]), deg, weight=1.d0/reform(T2_BD_err_cal[i,*])^2., residual=residual, predict_info=infoBD, /double)
    fitparams_BD[i,*] = param
    polfit_BD[i,*] = func_poly_deg2(hour_cal, param)

    ;confidence intervals for T^2 Cal BD
    result = imsl_polypredict(infoBD, hour_cal, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, y=reform(T2_AC_cal[i,*]), /double)
    interval_BD[i,*,*] = popmean

    ;confidence intervals for T^2 Sci BD
    result = imsl_polypredict(infoBD, hour_sci, ci_ptw_new_samp=newsamp, ci_ptw_pop_mean=popmean, ci_scheffe=scheffe,confidence=95.d0, /double)
    T2_BD_err_sci[i,*] = abs(reform(popmean[0,*])-reform(popmean[1,*]))/2.d0

  endfor

  ;------------------------------------------------------------------------------

  ;get TF for time of Sci observations

  T2_AC_sci = dblarr(npix, nsci) & T2_BD_sci = dblarr(npix, nsci)

  for i=0,npix-1 do begin

    for j=0,nsci-1 do begin

      T2_AC_sci[i,j] = func_poly_deg2(hour_sci[j], fitparams_AC[i,*])
      T2_BD_sci[i,j] = func_poly_deg2(hour_sci[j], fitparams_BD[i,*])

    endfor

  endfor

  st = {T2_AC_sci:T2_AC_sci, T2_AC_err_sci:T2_AC_err_sci, T2_BD_sci:T2_BD_sci, T2_BD_err_sci:T2_BD_err_sci, polfit_AC:polfit_AC, polfit_BD:polfit_BD, interval_AC:interval_AC, interval_BD:interval_BD}

  return, st

end


function tffunc_ave, ncal, nsci, hour_sci, hour_cal, T2_AC_cal, T2_BD_cal, T2_AC_err_cal, T2_BD_err_cal, npix

  polfit_AC = dblarr(npix,ncal) & polfit_BD = dblarr(npix,ncal)
  interval_AC = dblarr(npix,2,ncal) & interval_BD = dblarr(npix,2,ncal)
  T2_AC_sci = dblarr(npix, nsci) & T2_BD_sci = dblarr(npix, nsci)
  T2_AC_err_sci = dblarr(npix,nsci) & T2_BD_err_sci = dblarr(npix,nsci)

  for i=0,nsci-1 do begin

    idxl = closest2(hour_sci[i], hour_cal, /lower)
    idxu = closest2(hour_sci[i], hour_cal, /upper)

    if (idxl[0] eq -1) then idxl = idxu
    if (idxu[0] eq -1) then idxu = idxl

    if (idxl eq idxu) then print, 'Science is not bracketed by 2 Calibrators!'

    for j=0,npix-1 do begin

      T2_AC_sci[j,i] = mean([T2_AC_cal[j,idxl], T2_AC_cal[j,idxu]])
      T2_BD_sci[j,i] = mean([T2_BD_cal[j,idxl], T2_BD_cal[j,idxu]])

      T2_AC_err_sci[j,i] = stddev([T2_AC_cal[j,idxl], T2_AC_cal[j,idxu]])/sqrt(2.d0)
      T2_BD_err_sci[j,i] = stddev([T2_BD_cal[j,idxl], T2_BD_cal[j,idxu]])/sqrt(2.d0)

	if (idxl eq idxu) then begin	;in case we only have a cal-sci sequence

	  T2_AC_err_sci[j,i] = T2_AC_err_cal[j,idxl]
	  T2_BD_err_sci[j,i] = T2_BD_err_cal[j,idxl]

	endif

    endfor

  endfor

  polfit_AC = T2_AC_cal
  polfit_BD = T2_BD_cal

  interval_AC[*,0,*] = T2_AC_cal-T2_AC_err_cal
  interval_AC[*,1,*] = T2_AC_cal+T2_AC_err_cal

  interval_BD[*,0,*] = T2_BD_cal-T2_BD_err_cal
  interval_BD[*,1,*] = T2_BD_cal+T2_BD_err_cal

  st = {T2_AC_sci:T2_AC_sci, T2_AC_err_sci:T2_AC_err_sci, T2_BD_sci:T2_BD_sci, T2_BD_err_sci:T2_BD_err_sci, polfit_AC:polfit_AC, polfit_BD:polfit_BD, interval_AC:interval_AC, interval_BD:interval_BD}

  return, st

end



pro FSUA_CalibrateVisibility

path = '/home/amueller/work/MIDI_FSUA/Pipeline/V_FSUA/'

;==========================================================================

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

print, ''
print, 'Available Data Sets:'
for i=0,n_elements(dataset)-1 do print, dataset[i]
read, 'Select Data Set: ', dsselect

if (dsselect eq '1') then begin
  obsnum = 'FSUAscans_20130706'	;HD155826
  caldata = path+'20130706_HD155826/allFringeScans_HD155826.sav'
  xpr = [1.5,5.5]
endif
if (dsselect eq '2') then begin
  obsnum = 'FSUAscans_20131025'  ;42Cet
  caldata = path+'20131025_42Cet/allFringeScans_42Cet.sav'
  xpr = [4.0,6.0]
endif
if (dsselect eq '3') then begin
  obsnum = 'FSUAscans_20131028'	;24Psc
  caldata = path+'20131028_24Psc/allFringeScans_24Psc.sav'
  xpr = [4.5,7.0]
endif
if (dsselect eq '4') then begin
  obsnum = 'FSUAscans_20140113_U13'	;VelaX1
  caldata = path+'20140113_VelaX1/allFringeScans_VelaX1.sav'
  xpr = [3.0,5.0]
endif
if (dsselect eq '5') then begin
  obsnum = 'FSUAscans_20140113_U14'	;VelaX1
  caldata = path+'20140113_VelaX1/allFringeScans_VelaX1.sav'
  xpr = [5.0,7.0]
endif
if (dsselect eq '6') then begin
  obsnum = 'FSUAscans_20140306_1'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
  xpr = [1.0,4.5]
endif
if (dsselect eq '7') then begin
  obsnum = 'FSUAscans_20140306_2'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
  xpr = [1.0,4.5]
endif
if (dsselect eq '8') then begin
  obsnum = 'FSUAscans_20140306_3'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
  xpr = [1.0,4.5]
endif
if (dsselect eq '9') then begin
  obsnum = 'FSUAscans_20140307_1'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
  xpr = [1.0,4.5]
endif
if (dsselect eq '10') then begin
  obsnum = 'FSUAscans_20140307_2'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
  xpr = [1.0,4.5]
endif
if (dsselect eq '11') then begin
  obsnum = 'FSUAscans_20140307_3'
  caldata = path+'20140306_P92Rivi/allFringeScans_P92Rivi.sav'
  xpr = [1.0,4.5]
endif
if (dsselect eq '12') then begin
  obsnum = 'FSUAscans_20140412_1'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
  xpr = [3.0,10.]
endif
if (dsselect eq '13') then begin
  obsnum = 'FSUAscans_20140412_2'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
  xpr = [3.0,10.]
endif
if (dsselect eq '14') then begin
  obsnum = 'FSUAscans_20140412_3'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
  xpr = [3.0,10.]
endif
if (dsselect eq '15') then begin
  obsnum = 'FSUAscans_20140412_4'
  caldata = path+'20140412_P93KRAUSU14/allFringeScans_P93KRAUSU14.sav'
  xpr = [3.0,10.]
endif
if (dsselect eq '16') then begin
  obsnum = 'FSUAscans_20140415'
  caldata = path+'20140415_P93MENU/allFringeScans_P93MENU.sav'
  xpr = [4.0,6.5]
endif
if (dsselect eq '17') then begin
  obsnum = 'FSUAscans_20140413_1'
  caldata = path+'20140413_P93KRAUSU23/allFringeScans_P93KRAUSU23.sav'
  xpr = [2.5,4.75]
endif
if (dsselect eq '18') then begin
  obsnum = 'FSUAscans_20140413_2'
  caldata = path+'20140413_P93KRAUSU23/allFringeScans_P93KRAUSU23.sav'
  xpr = [4.5,10.]
endif
if (dsselect eq '19') then begin
  obsnum = 'FSUAscans_20140413'
  caldata = path+'20140413_P93KRAUSU23/allFringeScans_P93KRAUSU23.sav'
  xpr = [2.5,10.]
endif
if (dsselect eq '20') then begin
  obsnum = 'FSUAscans_20140706_A1C1'
  caldata = path+'20140706_P93MONNIER/allFringeScans_P93MONNIER.sav'
  xpr = [4.5,8.]
endif
if (dsselect eq '21') then begin
  obsnum = 'FSUAscans_20140706_B2C1'
  caldata = path+'20140706_P93MONNIER/allFringeScans_P93MONNIER.sav'
  xpr = [7.5,10.]
endif
if (dsselect eq '22') then begin
  obsnum = 'FSUAscans_P94KOEHLERH0I1'
  caldata = path+'20141023_P94KOEHLERH0I1/allFringeScans_P94KOEHLERH0I1.sav'
  xpr = [5.,10]
endif
if (dsselect eq '23') then begin
  obsnum = 'FSUAscans_P94KOEHLERD0I1'
  caldata = path+'20141024_P94KOEHLERD0I1/allFringeScans_P94KOEHLERD0I1.sav'
  xpr = [5,10.]
endif
;==========================================================================


pxmode = ''
print, ''
print, '3: 3 pixel'
print, '5: 5 pixel'
read, 'Select Pixel mode: ', pxmode

if (pxmode ne '3' and pxmode ne '5') then begin

  print, ''
  print, 'Choice does not exist. Stop.'
  print, ''
  return

endif

if (pxmode eq '3') then npix = 4
if (pxmode eq '5') then npix = 6

;select which TF mode
tfmode = ''
print, ''
print, '1: Average'
print, '2: Interpolation 2nd degree'
print, '3: Interpolation 3rd degree'
read, 'Select TF mode: ', tfmode

if (tfmode eq '1') then idtfmode = 'Average'
if (tfmode eq '2') then idtfmode = 'InterpDeg2'
if (tfmode eq '3') then idtfmode = 'InterpDeg3'

if (tfmode ne '1' and tfmode ne '2' and tfmode ne '3') then begin

  print, ''
  print, 'Choice does not exist. Stop.'
  print, ''
  return

endif


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

datapath = path+'Vis_'+obsnum+'_'+pxmode+'pixel/CoherenceFactor/'
resultpath = path+'Vis_'+obsnum+'_'+pxmode+'pixel/TFcal_V2sci_'+idtfmode+'/'
spawn, 'mkdir -p '+resultpath

restore, caldata
if (pxmode eq '5') then begin

  ave_wl_cal = [median(lam2[*,0],/even), median(lam2[*,1],/even), median(lam2[*,2],/even), median(lam2[*,3],/even), median(lam2[*,4],/even), median(lam2[*,5],/even)]
  bw =  [median(bandwidth[*,*,0],/even),median(bandwidth[*,*,1],/even),median(bandwidth[*,*,2],/even),median(bandwidth[*,*,3],/even),median(bandwidth[*,*,4],/even),median(bandwidth[*,*,5],/even)]

endif

if (pxmode eq '3') then begin

  ave_wl_cal = [median(lam2[*,0],/even), median(lam2[*,2],/even), median(lam2[*,3],/even), median(lam2[*,4],/even)]

  bw =  [median(bandwidth[*,*,0],/even),median(bandwidth[*,*,2],/even),median(bandwidth[*,*,3],/even),median(bandwidth[*,*,4],/even)]

endif

;------------------------------------------------------------------------------

micron = '!4' + String("154B) + '!Xm'

;define arrays for Cal and Sci

nfiles = n_elements(obsid)	;total number of files
idxcal = where(starflag eq 'Cal')
idxsci = where(starflag eq 'Sci')

obsid_cal = obsid[idxcal]
scan_cal = scan[idxcal]
starflag_cal = starflag[idxcal]
diam = diam[idxcal]
obsid_sci = obsid[idxsci]
scan_sci = scan[idxsci]
starflag_sci = starflag[idxsci]


ncal = n_elements(idxcal)	;number of cal files
nsci = n_elements(idxsci)	;number of sci files

target_cal = strarr(ncal)
jd_cal = dblarr(ncal)
coh_AC_cal = dblarr(npix, ncal) & coh_BD_cal = dblarr(npix, ncal)	;coherence factor for AC, BD calbrator
coh_AC_err_cal = dblarr(npix, ncal) & coh_BD_err_cal = dblarr(npix, ncal)
V2_cal = dblarr(npix, ncal)	;theoretical visibility based on stellar diamter calbrator
T2_AC_cal = dblarr(npix, ncal) & T2_BD_cal = dblarr(npix, ncal)	;transfer function for AC and BD calbrator
T2_AC_err_cal = dblarr(npix, ncal) & T2_BD_err_cal = dblarr(npix, ncal)

target_sci = strarr(nsci)
jd_sci = dblarr(nsci)
coh_AC_sci = dblarr(npix, nsci) & coh_BD_sci = dblarr(npix, nsci)	;coherence factor for AC, BD science
coh_AC_err_sci = dblarr(npix, nsci) & coh_BD_err_sci = dblarr(npix, nsci)

fileid_cal = strarr(ncal)
fileid_sci = strarr(nsci)

;------------------------------------------------------------------------------

;read in each single CAL file to get transfer function and coherence factor, etc

if (pxmode eq '3') then begin

  for i=0,ncal-1 do begin

    ;extract file id
    fileid_cal[i] = strmid(scan_cal[i], strpos(scan_cal[i], 'SCAN_')+5, strlen(scan_cal[i]))

    filetmp = file_search(datapath+fileid_cal[i]+'_'+starflag_cal[i]+'*results.txt')

    readcol, filetmp[0], targettmp, numline=1, format='a', /silent
    target_cal[i] = targettmp[0]

    readcol, filetmp[0], dum, jdtmp, numline=1, skipline=3, format='a,d', /silent
    jd_cal[i] = jdtmp[0]+2400000.d0

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=36, format='d,d,d,d,d,d', /silent
    coh_AC_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=39, format='d,d,d,d,d,d', /silent
    coh_BD_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=42, format='d,d,d,d,d,d', /silent
    coh_AC_err_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=43, format='d,d,d,d,d,d', /silent
    coh_BD_err_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=55, format='d,d,d,d,d,d', /silent, delimiter=': '
    V2_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=57, format='d,d,d,d,d,d', /silent
    T2_AC_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=59, format='d,d,d,d,d,d', /silent
    T2_AC_err_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=61, format='d,d,d,d,d,d', /silent
    T2_BD_cal[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=63, format='d,d,d,d,d,d', /silent
    T2_BD_err_cal[*,i] = [t1, t2, t3, t4]

  endfor

;------------------------------------------------------------------------------

  ;read in each single SCI file to get coherence factor, etc
  for i=0,nsci-1 do begin

    ;extract file id
    fileid_sci[i] = strmid(scan_sci[i], strpos(scan_sci[i], 'SCAN_')+5, strlen(scan_sci[i]))

    ;read in each single file to get transfer function and coherence factor
    filetmp = file_search(datapath+fileid_sci[i]+'_'+starflag_sci[i]+'*results.txt')

    readcol, filetmp[0], targettmp, numline=1, format='a', /silent
    target_sci[i] = targettmp[0]

    readcol, filetmp[0], dum, jdtmp, numline=1, skipline=3, format='a,d', /silent
    jd_sci[i] = jdtmp[0]+2400000.d0

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=36, format='d,d,d,d,d,d', /silent
    coh_AC_sci[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=39, format='d,d,d,d,d,d', /silent
    coh_BD_sci[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=42, format='d,d,d,d,d,d', /silent
    coh_AC_err_sci[*,i] = [t1, t2, t3, t4]

    readcol, filetmp[0], t1, t2, t3, t4, numline=1, skipline=43, format='d,d,d,d,d,d', /silent
    coh_BD_err_sci[*,i] = [t1, t2, t3, t4]

  endfor

endif

if (pxmode eq '5') then begin

  for i=0,ncal-1 do begin

    ;extract file id
    fileid_cal[i] = strmid(scan_cal[i], strpos(scan_cal[i], 'SCAN_')+5, strlen(scan_cal[i]))

    filetmp = file_search(datapath+fileid_cal[i]+'_'+starflag_cal[i]+'*results.txt')

    readcol, filetmp[0], targettmp, numline=1, format='a', /silent
    target_cal[i] = targettmp[0]

    readcol, filetmp[0], dum, jdtmp, numline=1, skipline=3, format='a,d', /silent
    jd_cal[i] = jdtmp[0];+2400000.d0

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=36, format='d,d,d,d,d,d', /silent
    coh_AC_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=39, format='d,d,d,d,d,d', /silent
    coh_BD_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=42, format='d,d,d,d,d,d', /silent
    coh_AC_err_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=43, format='d,d,d,d,d,d', /silent
    coh_BD_err_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=55, format='d,d,d,d,d,d', /silent, delimiter=': '
    V2_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=57, format='d,d,d,d,d,d', /silent
    T2_AC_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=59, format='d,d,d,d,d,d', /silent
    T2_AC_err_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=61, format='d,d,d,d,d,d', /silent
    T2_BD_cal[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=63, format='d,d,d,d,d,d', /silent
    T2_BD_err_cal[*,i] = [t1, t2, t3, t4, t5, t6]

  endfor

  ;------------------------------------------------------------------------------

  ;read in each single SCI file to get coherence factor, etc
  for i=0,nsci-1 do begin

    ;extract file id
    fileid_sci[i] = strmid(scan_sci[i], strpos(scan_sci[i], 'SCAN_')+5, strlen(scan_sci[i]))

    ;read in each single file to get transfer function and coherence factor
    filetmp = file_search(datapath+fileid_sci[i]+'_'+starflag_sci[i]+'*results.txt')

    readcol, filetmp[0], targettmp, numline=1, format='a', /silent
    target_sci[i] = targettmp[0]

    readcol, filetmp[0], dum, jdtmp, numline=1, skipline=3, format='a,d', /silent
    jd_sci[i] = jdtmp[0];+2400000.d0

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=36, format='d,d,d,d,d,d', /silent
    coh_AC_sci[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=39, format='d,d,d,d,d,d', /silent
    coh_BD_sci[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=42, format='d,d,d,d,d,d', /silent
    coh_AC_err_sci[*,i] = [t1, t2, t3, t4, t5, t6]

    readcol, filetmp[0], t1, t2, t3, t4, t5, t6, numline=1, skipline=43, format='d,d,d,d,d,d', /silent
    coh_BD_err_sci[*,i] = [t1, t2, t3, t4, t5, t6]

  endfor

endif

;------------------------------------------------------------------------------

;get UT from MJD !!!2.5min difference between MJD-OBS and UTC in header!!!
caldat, jd_cal, month, day, year, hh, mm, ss
hour_cal = hh+mm/60.d0+ss/3600.d0
caldat, jd_sci, month, day, year, hh, mm, ss
hour_sci = hh+mm/60.d0+ss/3600.d0

hour = [hour_cal, hour_sci]
hour = hour[sort(hour)]

;get the date of the night
minjd = min([jd_cal, jd_sci])-0.5d0
caldat, minjd, month, day, year
if (month le 9) then month = '0'+strcompress(month,/rem) else month = strcompress(month,/rem)
if (day le 9) then day = '0'+strcompress(day,/rem) else day = strcompress(day,/rem)
date = day+'.'+month+'.'+strcompress(year,/rem)

;------------------------------------------------------------------------------

;get TF for science observations

if (tfmode eq '1') then $
  result = tffunc_ave(ncal, nsci, hour_sci, hour_cal, T2_AC_cal, T2_BD_cal, T2_AC_err_cal, T2_BD_err_cal,  npix)

if (tfmode eq '2') then $
  result = tffunc_interp2(ncal, nsci, hour_sci, hour_cal, T2_AC_cal, T2_BD_cal, T2_AC_err_cal, T2_BD_err_cal,  npix)

if (tfmode eq '3') then $
  result = tffunc_interp3(ncal, nsci, hour_sci, hour_cal, T2_AC_cal, T2_BD_cal, T2_AC_err_cal, T2_BD_err_cal, npix)


T2_AC_sci = result.T2_AC_sci
T2_AC_err_sci = result.T2_AC_err_sci
T2_BD_sci = result.T2_BD_sci
T2_BD_err_sci = result.T2_BD_err_sci
polfit_AC = result.polfit_AC
polfit_BD = result.polfit_BD
interval_AC = result.interval_AC
interval_BD = result.interval_BD


;------------------------------------------------------------------------------

;compute calibrated V^2 for science for each pixel

V2_AC_sci = dblarr(npix,nsci) & V2_BD_sci = dblarr(npix,nsci)
V2_AC_err_sci = dblarr(npix,nsci) & V2_BD_err_sci = dblarr(npix,nsci)

for i=0,npix-1 do begin

  for j=0,nsci-1 do begin

    V2_AC_sci[i,j] = coh_AC_sci[i,j]/T2_AC_sci[i,j]
    V2_BD_sci[i,j] = coh_BD_sci[i,j]/T2_BD_sci[i,j]

    V2_AC_err_sci[i,j] = V2_AC_sci[i,j]*sqrt((coh_AC_err_sci[i,j]/coh_AC_sci[i,j])^2. + (T2_AC_err_sci[i,j]/T2_AC_sci[i,j])^2.)
    V2_BD_err_sci[i,j] = V2_BD_sci[i,j]*sqrt((coh_BD_err_sci[i,j]/coh_BD_sci[i,j])^2. + (T2_BD_err_sci[i,j]/T2_BD_sci[i,j])^2.)

  endfor

endfor

;--------------------------------------------------------------------------------------------
;plot of transfer function for AC and BD

if (pxmode eq '3') then begin

  set_plot, 'ps'
  device, filename=resultpath+'TFcal_V2sci_'+idtfmode+'.ps',/color,XSIZE=45, YSIZE=30, XOffset=xoffset, YOffset=yoffset
  !p.thick=4
  !x.thick=3
  !y.thick=3
  !p.multi=[0,2,2]

  ;A-C
  xmin = floor(min(hour_cal))-0.2
  xmax = ceil(max(hour_cal))+0.2
;   ymin = floor(min(T2_AC_cal)*10.)/10.d0
;   ymax = ceil(max(T2_AC_cal)*10.)/10.d0
  ploterror, hour_cal, T2_AC_cal[0,*], T2_AC_err_cal[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), ytitle='Transfer Function T!U2!N', xtitle='UT / [hour]', title='Calibrator, AC', $; title='Calibrator '+date+' for AC', $
  charsize=1.5, xr=xpr, xst=1, yst=1, yr=[0,1.1], $	;, yr=[ymin, ymax]
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_cal, T2_AC_cal[0,*], T2_AC_err_cal[0,*], psym=sym(1), color=fsc_color('black'), symsize=1.
  oplot, hour_cal, polfit_AC[0,*], color=fsc_color('black')
  oplot, hour_cal, interval_AC[0,0,*], color=fsc_color('black'), linestyle=2
  oplot, hour_cal, interval_AC[0,1,*], color=fsc_color('black'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[1,*], T2_AC_err_cal[1,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oplot, hour_cal, polfit_AC[1,*], color=fsc_color('blue')
  oplot, hour_cal, interval_AC[1,0,*], color=fsc_color('blue'), linestyle=2
  oplot, hour_cal, interval_AC[1,1,*], color=fsc_color('blue'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[2,*], T2_AC_err_cal[2,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oplot, hour_cal, polfit_AC[2,*], color=fsc_color('green')
  oplot, hour_cal, interval_AC[2,0,*], color=fsc_color('green'), linestyle=2
  oplot, hour_cal, interval_AC[2,1,*], color=fsc_color('green'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[3,*], T2_AC_err_cal[3,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
  oplot, hour_cal, polfit_AC[3,*], color=fsc_color('orange')
  oplot, hour_cal, interval_AC[3,0,*], color=fsc_color('orange'), linestyle=2
  oplot, hour_cal, interval_AC[3,1,*], color=fsc_color('orange'), linestyle=2

  legend, ['White', '2.1'+micron, '2.2'+micron, '2.3'+micron], $
  margin=0, box=0, /left, /top, $
  color=[fsc_color('black'), fsc_color('blue'), fsc_color('green'), fsc_color('orange')], $
  psym=[sym(1),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5
  ;legend, ['95% Confidence'], margin=0, box=0, /right, /top, color=fsc_color('black'), linestyle=2, textcolor=fsc_color('black'), charsize=1.5, /clear

  ;oplot computed T^2 for science
  oploterror, hour_sci, T2_AC_sci[0,*], T2_AC_err_sci[0,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[1,*], T2_AC_err_sci[1,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[2,*], T2_AC_err_sci[2,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[3,*], T2_AC_err_sci[3,*], psym=sym(1), symsize=1., color=fsc_color('gray')

;------------------------------------------

  ;B-D
;   ymin = floor(min(T2_BD_cal)*10.)/10.d0
;   ymax = ceil(max(T2_BD_cal)*10.)/10.d0
  plot, hour_cal, T2_BD_cal[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), ytitle='Transfer Function T!U2!N', xtitle='UT / [hour]', title='Calibrator, BD', $;, title='Calibrator '+date+' for BD'
  charsize=1.5, xr=xpr, xst=1, yst=1, yr=[0,1.1], $	;, yr=[ymin, ymax]
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_cal, T2_BD_cal[0,*], T2_BD_err_cal[0,*], psym=sym(1), color=fsc_color('black'), symsize=1.
  oplot, hour_cal, polfit_BD[0,*], color=fsc_color('black')
  oplot, hour_cal, interval_BD[0,0,*], color=fsc_color('black'), linestyle=2
  oplot, hour_cal, interval_BD[0,1,*], color=fsc_color('black'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[1,*], T2_BD_err_cal[1,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oplot, hour_cal, polfit_BD[1,*], color=fsc_color('blue')
  oplot, hour_cal, interval_BD[1,0,*], color=fsc_color('blue'), linestyle=2
  oplot, hour_cal, interval_BD[1,1,*], color=fsc_color('blue'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[2,*], T2_BD_err_cal[2,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oplot, hour_cal, polfit_BD[2,*], color=fsc_color('green')
  oplot, hour_cal, interval_BD[2,0,*], color=fsc_color('green'), linestyle=2
  oplot, hour_cal, interval_BD[2,1,*], color=fsc_color('green'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[3,*], T2_BD_err_cal[3,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
  oplot, hour_cal, polfit_BD[3,*], color=fsc_color('orange')
  oplot, hour_cal, interval_BD[3,0,*], color=fsc_color('orange'), linestyle=2
  oplot, hour_cal, interval_BD[3,1,*], color=fsc_color('orange'), linestyle=2

;   legend, ['White', '2.1'+micron, '2.2'+micron, '2.3'+micron], $
;   margin=0, box=0, /left, /top, $
;   color=[fsc_color('black'), fsc_color('blue'), fsc_color('green'), fsc_color('orange')], $
;   psym=[sym(1),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5
  ;legend, ['95% Confidence'], margin=0, box=0, /right, /top, color=fsc_color('black'), linestyle=2, textcolor=fsc_color('black'), charsize=1.5, /clear

  ;oplot computed T^2 for science
  oploterror, hour_sci, T2_BD_sci[0,*], T2_BD_err_sci[0,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[1,*], T2_BD_err_sci[1,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[2,*], T2_BD_err_sci[2,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[3,*], T2_BD_err_sci[3,*], psym=sym(1), symsize=1., color=fsc_color('gray')

;--------------------------------------------------------------------------------------------
;plot of calibrated science visibility for AC and BD

  xmin = floor(min(hour_sci))-0.2
  xmax = ceil(max(hour_sci))+0.2
  ymin = floor(min(V2_AC_sci)*10.)/10.d0
  ymax = ceil(max(V2_AC_sci)*10.)/10.d0

  ;AC
  plot, hour_sci, V2_AC_sci[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), $
  ytitle='Science Calibrated Visibility V!U2!N', xtitle='UT / [hour]', title=target_sci[0]+', AC', charsize=1.5, $
  xr=xpr, xst=1, yr=[0,1.1], yst=1, $
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_sci, V2_AC_sci[0,*], V2_AC_err_sci[0,*], psym=sym(4), color=fsc_color('black'), symsize=3
;   oploterror, hour_sci, V2_AC_sci[1,*], V2_AC_err_sci[1,*], psym=sym(1), color=fsc_color('violet'), symsize=1.
  oploterror, hour_sci, V2_AC_sci[1,*], V2_AC_err_sci[1,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oploterror, hour_sci, V2_AC_sci[2,*], V2_AC_err_sci[2,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oploterror, hour_sci, V2_AC_sci[3,*], V2_AC_err_sci[3,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
;   oploterror, hour_sci, V2_AC_sci[5,*], V2_AC_err_sci[5,*], psym=sym(1), color=fsc_color('red'), symsize=1.

  plots, !x.crange, [1,1], thick=2, color=fsc_color('black')
  for i=0,n_elements(target_sci)-1 do xyouts, hour_sci[i], 1.01, target_sci[i], /data, color=fsc_color('black'), charsize=1.2

;   legend, ['White', '2.1'+micron, '2.2'+micron, '2.3'+micron], $
;   margin=0, box=0, /left, /top, $
;   color=[fsc_color('black'), fsc_color('blue'), fsc_color('green'), fsc_color('orange')], $
;   psym=[sym(1),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5

  ;BD
  ymin = floor(min(V2_BD_sci)*10.)/10.d0
  ymax = ceil(max(V2_BD_sci)*10.)/10.d0

  plot, hour_sci, V2_BD_sci[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), $
  ytitle='Science Calibrated Visibility V!U2!N', xtitle='UT / [hour]', title=target_sci[0]+', BD', charsize=1.5, $
  xr=xpr, xst=1, yr=[0,1.1], yst=1, $
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_sci, V2_BD_sci[0,*], V2_BD_err_sci[0,*], psym=sym(4), color=fsc_color('black'), symsize=3
  oploterror, hour_sci, V2_BD_sci[1,*], V2_BD_err_sci[1,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oploterror, hour_sci, V2_BD_sci[2,*], V2_BD_err_sci[2,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oploterror, hour_sci, V2_BD_sci[3,*], V2_BD_err_sci[3,*], psym=sym(1), color=fsc_color('orange'), symsize=1.

  plots, !x.crange, [1,1], thick=2, color=fsc_color('black')
  for i=0,n_elements(target_sci)-1 do xyouts, hour_sci[i], 1.01, target_sci[i], /data, color=fsc_color('black'), charsize=1.2

;   legend, ['White', '2.1'+micron, '2.2'+micron, '2.3'+micron], $
;   margin=0, box=0, /left, /top, $
;   color=[fsc_color('black'), fsc_color('blue'), fsc_color('green'), fsc_color('orange')], $
;   psym=[sym(4),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5

  device,/close
  set_plot,'x'

;   spawn, 'gv '+resultpath+'TFcal_V2sci.ps'

  !p.multi=[0,1,0]
  ; !P.Font=0
  !p.thick=1
  !x.thick=1
  !y.thick=1

endif

if (pxmode eq '5') then begin

  set_plot, 'ps'
  device, filename=resultpath+'TFcal_V2sci_'+idtfmode+'.ps',/color,XSIZE=45, YSIZE=30, XOffset=xoffset, YOffset=yoffset
  !p.thick=4
  !x.thick=3
  !y.thick=3
  !p.multi=[0,2,2]


  ;A-C
  xmin = floor(min(hour_cal))-0.2
  xmax = ceil(max(hour_cal))+0.2
;   ymin = floor(min(T2_AC_cal)*10.)/10.d0
;   ymax = ceil(max(T2_AC_cal)*10.)/10.d0
  ploterror, hour_cal, T2_AC_cal[0,*], T2_AC_err_cal[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), ytitle='Transfer Function T!U2!N', xtitle='UT / [hour]', title='Calibrator, AC', $; title='Calibrator '+date+' for AC', $
  charsize=1.5, xr=xpr, xst=1, yst=1, yr=[0,1.1], $	;, yr=[ymin, ymax]
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_cal, T2_AC_cal[0,*], T2_AC_err_cal[0,*], psym=sym(1), color=fsc_color('black'), symsize=1.
  oplot, hour_cal, polfit_AC[0,*], color=fsc_color('black')
  oplot, hour_cal, interval_AC[0,0,*], color=fsc_color('black'), linestyle=2
  oplot, hour_cal, interval_AC[0,1,*], color=fsc_color('black'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[1,*], T2_AC_err_cal[1,*], psym=sym(1), color=fsc_color('violet'), symsize=1.
  oplot, hour_cal, polfit_AC[1,*], color=fsc_color('violet')
  oplot, hour_cal, interval_AC[1,0,*], color=fsc_color('violet'), linestyle=2
  oplot, hour_cal, interval_AC[1,1,*], color=fsc_color('violet'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[2,*], T2_AC_err_cal[2,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oplot, hour_cal, polfit_AC[2,*], color=fsc_color('blue')
  oplot, hour_cal, interval_AC[2,0,*], color=fsc_color('blue'), linestyle=2
  oplot, hour_cal, interval_AC[2,1,*], color=fsc_color('blue'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[3,*], T2_AC_err_cal[3,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oplot, hour_cal, polfit_AC[3,*], color=fsc_color('green')
  oplot, hour_cal, interval_AC[3,0,*], color=fsc_color('green'), linestyle=2
  oplot, hour_cal, interval_AC[3,1,*], color=fsc_color('green'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[4,*], T2_AC_err_cal[4,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
  oplot, hour_cal, polfit_AC[4,*], color=fsc_color('orange')
  oplot, hour_cal, interval_AC[4,0,*], color=fsc_color('orange'), linestyle=2
  oplot, hour_cal, interval_AC[4,1,*], color=fsc_color('orange'), linestyle=2

  oploterror, hour_cal, T2_AC_cal[5,*], T2_AC_err_cal[5,*], psym=sym(1), color=fsc_color('red'), symsize=1.
  oplot, hour_cal, polfit_AC[5,*], color=fsc_color('red')
  oplot, hour_cal, interval_AC[5,0,*], color=fsc_color('red'), linestyle=2
  oplot, hour_cal, interval_AC[5,1,*], color=fsc_color('red'), linestyle=2

  legend, ['White', '2.0'+micron, '2.1'+micron, '2.2'+micron, '2.3'+micron, '2.4'+micron], $
  margin=-.5, box=0, /left, /top, $
  color=[fsc_color('black'), fsc_color('violet'), fsc_color('blue'), fsc_color('green'), fsc_color('orange'), fsc_color('red')], $
  psym=[sym(1),sym(1),sym(1),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5
  ;legend, ['95% Confidence'], margin=0, box=0, /right, /top, color=fsc_color('black'), linestyle=2, textcolor=fsc_color('black'), charsize=1.5, /clear

  ;oplot computed T^2 for science
  oploterror, hour_sci, T2_AC_sci[0,*], T2_AC_err_sci[0,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[1,*], T2_AC_err_sci[1,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[2,*], T2_AC_err_sci[2,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[3,*], T2_AC_err_sci[3,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[4,*], T2_AC_err_sci[4,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_AC_sci[5,*], T2_AC_err_sci[5,*], psym=sym(1), symsize=1., color=fsc_color('gray')

;------------------------------------------

  ;B-D
;   ymin = floor(min(T2_BD_cal)*10.)/10.d0
;   ymax = ceil(max(T2_BD_cal)*10.)/10.d0
  plot, hour_cal, T2_BD_cal[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), ytitle='Transfer Function T!U2!N', xtitle='UT / [hour]', title='Calibrator, BD', $; title='Calibrator '+date+' for BD', $
  charsize=1.5, xr=xpr, xst=1, yst=1, yr=[0,1.1], $	;, yr=[ymin, ymax]
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_cal, T2_BD_cal[0,*], T2_BD_err_cal[0,*], psym=sym(1), color=fsc_color('black'), symsize=1.
  oplot, hour_cal, polfit_BD[0,*], color=fsc_color('black')
  oplot, hour_cal, interval_BD[0,0,*], color=fsc_color('black'), linestyle=2
  oplot, hour_cal, interval_BD[0,1,*], color=fsc_color('black'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[1,*], T2_BD_err_cal[1,*], psym=sym(1), color=fsc_color('violet'), symsize=1.
  oplot, hour_cal, polfit_BD[1,*], color=fsc_color('violet')
  oplot, hour_cal, interval_BD[1,0,*], color=fsc_color('violet'), linestyle=2
  oplot, hour_cal, interval_BD[1,1,*], color=fsc_color('violet'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[2,*], T2_BD_err_cal[2,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oplot, hour_cal, polfit_BD[2,*], color=fsc_color('blue')
  oplot, hour_cal, interval_BD[2,0,*], color=fsc_color('blue'), linestyle=2
  oplot, hour_cal, interval_BD[2,1,*], color=fsc_color('blue'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[3,*], T2_BD_err_cal[3,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oplot, hour_cal, polfit_BD[3,*], color=fsc_color('green')
  oplot, hour_cal, interval_BD[3,0,*], color=fsc_color('green'), linestyle=2
  oplot, hour_cal, interval_BD[3,1,*], color=fsc_color('green'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[4,*], T2_BD_err_cal[4,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
  oplot, hour_cal, polfit_BD[4,*], color=fsc_color('orange')
  oplot, hour_cal, interval_BD[4,0,*], color=fsc_color('orange'), linestyle=2
  oplot, hour_cal, interval_BD[4,1,*], color=fsc_color('orange'), linestyle=2

  oploterror, hour_cal, T2_BD_cal[5,*], T2_BD_err_cal[5,*], psym=sym(1), color=fsc_color('red'), symsize=1.
  oplot, hour_cal, polfit_BD[5,*], color=fsc_color('red')
  oplot, hour_cal, interval_BD[5,0,*], color=fsc_color('red'), linestyle=2
  oplot, hour_cal, interval_BD[5,1,*], color=fsc_color('red'), linestyle=2

;   legend, ['White', '2.0'+micron, '2.1'+micron, '2.2'+micron, '2.3'+micron, '2.4'+micron], $
;   margin=0, box=0, /left, /top, $
;   color=[fsc_color('black'), fsc_color('violet'), fsc_color('blue'), fsc_color('green'), fsc_color('orange'), fsc_color('red')], $
;   psym=[sym(1),sym(1),sym(1),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5
  ;legend, ['95% Confidence'], margin=0, box=0, /right, /top, color=fsc_color('black'), linestyle=2, textcolor=fsc_color('black'), charsize=1.5, /clear

  ;oplot computed T^2 for science
  oploterror, hour_sci, T2_BD_sci[0,*], T2_BD_err_sci[0,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[1,*], T2_BD_err_sci[1,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[2,*], T2_BD_err_sci[2,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[3,*], T2_BD_err_sci[3,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[4,*], T2_BD_err_sci[4,*], psym=sym(1), symsize=1., color=fsc_color('gray')
  oploterror, hour_sci, T2_BD_sci[5,*], T2_BD_err_sci[5,*], psym=sym(1), symsize=1., color=fsc_color('gray')

;--------------------------------------------------------------------------------------------
;plot of calibrated science visibility for AC and BD

  xmin = floor(min(hour_sci))-0.2
  xmax = ceil(max(hour_sci))+0.2
  ymin = floor(min(V2_AC_sci)*10.)/10.d0
  ymax = ceil(max(V2_AC_sci)*10.)/10.d0

  ;AC
  plot, hour_sci, V2_AC_sci[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), $
  ytitle='Science Calibrated Visibility V!U2!N', xtitle='UT / [hour]', title='Science, AC', charsize=1.5, $
  xr=xpr, xst=1, yr=[0,1.1], yst=1, $
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_sci, V2_AC_sci[0,*], V2_AC_err_sci[0,*], psym=sym(4), color=fsc_color('black'), symsize=3
  oploterror, hour_sci, V2_AC_sci[1,*], V2_AC_err_sci[1,*], psym=sym(1), color=fsc_color('violet'), symsize=1.
  oploterror, hour_sci, V2_AC_sci[2,*], V2_AC_err_sci[2,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oploterror, hour_sci, V2_AC_sci[3,*], V2_AC_err_sci[3,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oploterror, hour_sci, V2_AC_sci[4,*], V2_AC_err_sci[4,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
  oploterror, hour_sci, V2_AC_sci[5,*], V2_AC_err_sci[5,*], psym=sym(1), color=fsc_color('red'), symsize=1.

  plots, !x.crange, [1,1], thick=2, color=fsc_color('black')
  for i=0,n_elements(target_sci)-1 do xyouts, hour_sci[i], 1.01, target_sci[i], /data, color=fsc_color('black'), charsize=1.2

;   legend, ['White', '2.0'+micron, '2.1'+micron, '2.2'+micron, '2.3'+micron, '2.4'+micron], $
;   margin=0, box=0, /left, /top, $
;   color=[fsc_color('black'), fsc_color('violet'), fsc_color('blue'), fsc_color('green'), fsc_color('orange'), fsc_color('red')], $
;   psym=[sym(1),sym(1),sym(1),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5

;test only where a calibrator is treated as science
; idx = where(target_cal eq 'HD152161')
; ;idx = where(target_cal eq 'HD152334')
; ; idx = where(target_cal eq 'HD146051')
; hour_sci = hour_cal[idx]
; coh_AC_sci = coh_AC_cal[*,idx]
; T2_AC_sci = dblarr(6, 4)
; V2_AC_sci = dblarr(6, 4)
; for i=0,5 do begin
;   for j=0,4-1 do begin
;     T2_AC_sci[i,j] = func_poly_deg3(hour_sci[j], fitparams_AC[i,*])
;   endfor
; endfor
; for i=0,5 do begin
;   for j=0,4-1 do begin
;     V2_AC_sci[i,j] = coh_AC_sci[i,j]/T2_AC_sci[i,j]
;   endfor
; endfor
; oplot, hour_sci, V2_AC_sci[0,*], psym=sym(1), color=fsc_color('black'), symsize=1.
; oplot, hour_sci, V2_AC_sci[1,*], psym=sym(1), color=fsc_color('violet'), symsize=1.
; oplot, hour_sci, V2_AC_sci[2,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
; oplot, hour_sci, V2_AC_sci[3,*], psym=sym(1), color=fsc_color('green'), symsize=1.
; oplot, hour_sci, V2_AC_sci[4,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
; oplot, hour_sci, V2_AC_sci[5,*], psym=sym(1), color=fsc_color('red'), symsize=1.
; stop

  ;BD
  ymin = floor(min(V2_BD_sci)*10.)/10.d0
  ymax = ceil(max(V2_BD_sci)*10.)/10.d0

  plot, hour_sci, V2_BD_sci[0,*], /nodata, $
  background=fsc_color('white'), color=fsc_color('black'), $
  ytitle='Science Calibrated Visibility V!U2!N', xtitle='UT / [hour]', title='Science, BD', charsize=1.5, $
  xr=xpr, xst=1, yr=[0,1.1], yst=1, $
  XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1

  oploterror, hour_sci, V2_BD_sci[0,*], V2_BD_err_sci[0,*], psym=sym(4), color=fsc_color('black'), symsize=3
  oploterror, hour_sci, V2_BD_sci[1,*], V2_BD_err_sci[1,*], psym=sym(1), color=fsc_color('violet'), symsize=1.
  oploterror, hour_sci, V2_BD_sci[2,*], V2_BD_err_sci[2,*], psym=sym(1), color=fsc_color('blue'), symsize=1.
  oploterror, hour_sci, V2_BD_sci[3,*], V2_BD_err_sci[3,*], psym=sym(1), color=fsc_color('green'), symsize=1.
  oploterror, hour_sci, V2_BD_sci[4,*], V2_BD_err_sci[4,*], psym=sym(1), color=fsc_color('orange'), symsize=1.
  oploterror, hour_sci, V2_BD_sci[5,*], V2_BD_err_sci[5,*], psym=sym(1), color=fsc_color('red'), symsize=1.

  plots, !x.crange, [1,1], thick=2, color=fsc_color('black')
  for i=0,n_elements(target_sci)-1 do xyouts, hour_sci[i], 1.01, target_sci[i], /data, color=fsc_color('black'), charsize=1.2

;   legend, ['White', '2.0'+micron, '2.1'+micron, '2.2'+micron, '2.3'+micron, '2.4'+micron], $
;   margin=0, box=0, /left, /top, $
;   color=[fsc_color('black'), fsc_color('violet'), fsc_color('blue'), fsc_color('green'), fsc_color('orange'), fsc_color('red')], $
;   psym=[sym(4),sym(1),sym(1),sym(1),sym(1),sym(1)], textcolor=fsc_color('black'), charsize=1.5

  device,/close
  set_plot,'x'

;   spawn, 'gv '+resultpath+'TFcal_V2sci.ps'

  !p.multi=[0,1,0]
  ; !P.Font=0
  !p.thick=1
  !x.thick=1
  !y.thick=1

endif

;------------------------------------------------------------------------------

;write out the results for the inidvidual scans

if (pxmode eq '3') then begin

  pixel = ['0','2','3','4']

  for i=0,ncal-1 do begin
  
    openw, lun, resultpath+fileid_cal[i]+'_'+starflag_cal[i]+'_'+target_cal[i]+'_TFcal.txt', width=1400, /get_lun
    printf, lun, target_cal[i]
    printf, lun, starflag_cal[i]
    printf, lun, 'JD: ', jd_cal[i], format='(a4, f16.8)'
    printf, lun, 'Hour / UT: ', hour_cal[i], format='(a11, f10.7)'
    printf, lun, ''
    printf, lun, 'Median combined wavelengths and bandwidths for all 4 pixel (White-2-3-4) / [m]'
    printf, lun, pixel[0], ave_wl_cal[0], bw[0], format='(a4,2e15.7)'
    printf, lun, pixel[1], ave_wl_cal[1], bw[1], format='(a4,2e15.7)'
    printf, lun, pixel[2], ave_wl_cal[2], bw[2], format='(a4,2e15.7)'
    printf, lun, pixel[3], ave_wl_cal[3], bw[3], format='(a4,2e15.7)'
    printf, lun, ''
    printf, lun, 'Pixel       TF^2 AC     TF_err^2 AC         TF^2 BD     TF_err^2 BD'
    for j=0,3 do printf, lun, pixel[j], T2_AC_cal[j,i], T2_AC_err_cal[j,i], T2_BD_cal[j,i], T2_BD_err_cal[j,i], format='(a3,4f16.10)'

    close, lun
    free_lun, lun

  endfor

  for i=0,nsci-1 do begin
  
    openw, lun, resultpath+fileid_sci[i]+'_'+starflag_sci[i]+'_'+target_sci[i]+'_V2sci.txt', width=1400, /get_lun
    printf, lun, target_sci[i]
    printf, lun, starflag_sci[i]
    printf, lun, 'JD: ', jd_sci[i], format='(a4, f16.8)'
    printf, lun, 'Hour / UT: ', hour_sci[i], format='(a11, f10.7)'
    printf, lun, ''
    printf, lun, 'Median combined wavelengths and bandwidths for all 6 pixel (White-2-3-4) / [m]'
    printf, lun, pixel[0], ave_wl_cal[0], bw[0], format='(a4,2e15.7)'
    printf, lun, pixel[1], ave_wl_cal[1], bw[1], format='(a4,2e15.7)'
    printf, lun, pixel[2], ave_wl_cal[2], bw[2], format='(a4,2e15.7)'
    printf, lun, pixel[3], ave_wl_cal[3], bw[3], format='(a4,2e15.7)'
    printf, lun, ''
    printf, lun, 'Pixel       TF^2 AC     TF_err^2 AC           V2 AC       V2_err AC         TF^2 BD     TF_err^2 BD           V2 BD       V2_err BD'
    for j=0,3 do printf, lun, pixel[j], T2_AC_sci[j,i], T2_AC_err_sci[j,i], V2_AC_sci[j,i], V2_AC_err_sci[j,i], T2_BD_sci[j,i], T2_BD_err_sci[j,i], V2_BD_sci[j,i], V2_BD_err_sci[j,i], format='(a3,8f16.10)'


    close, lun
    free_lun, lun

  endfor

endif

if (pxmode eq '5') then begin

  pixel = ['0','1','2','3','4','5']

  for i=0,ncal-1 do begin
  
    openw, lun, resultpath+fileid_cal[i]+'_'+starflag_cal[i]+'_'+target_cal[i]+'_TFcal.txt', width=1400, /get_lun
    printf, lun, target_cal[i]
    printf, lun, starflag_cal[i]
    printf, lun, 'JD: ', jd_cal[i], format='(a4, f16.8)'
    printf, lun, 'Hour / UT: ', hour_cal[i], format='(a11, f10.7)'
    printf, lun, ''
    printf, lun, 'Median combined wavelengths and bandwidths for all 6 pixel (White-1-2-3-4-5) / [m]'
    printf, lun, pixel[0], ave_wl_cal[0], bw[0], format='(a4,2e15.7)'
    printf, lun, pixel[1], ave_wl_cal[1], bw[1], format='(a4,2e15.7)'
    printf, lun, pixel[2], ave_wl_cal[2], bw[2], format='(a4,2e15.7)'
    printf, lun, pixel[3], ave_wl_cal[3], bw[3], format='(a4,2e15.7)'
    printf, lun, pixel[4], ave_wl_cal[4], bw[4], format='(a4,2e15.7)'
    printf, lun, pixel[5], ave_wl_cal[5], bw[5], format='(a4,2e15.7)'
    printf, lun, ''
    printf, lun, 'Pixel       TF^2 AC     TF_err^2 AC         TF^2 BD     TF_err^2 BD'
    for j=0,5 do printf, lun, pixel[j], T2_AC_cal[j,i], T2_AC_err_cal[j,i], T2_BD_cal[j,i], T2_BD_err_cal[j,i], format='(a3,4f16.10)'

    close, lun
    free_lun, lun

  endfor

  for i=0,nsci-1 do begin
  
    openw, lun, resultpath+fileid_sci[i]+'_'+starflag_sci[i]+'_'+target_sci[i]+'_V2sci.txt', width=1400, /get_lun
    printf, lun, target_sci[i]
    printf, lun, starflag_sci[i]
    printf, lun, 'JD: ', jd_sci[i], format='(a4, f16.8)'
    printf, lun, 'Hour / UT: ', hour_sci[i], format='(a11, f10.7)'
    printf, lun, ''
    printf, lun, 'Median combined wavelengths and bandwidths for all 6 pixel (White-1-2-3-4-5) / [m]'
    printf, lun, pixel[0], ave_wl_cal[0], bw[0], format='(a4,2e15.7)'
    printf, lun, pixel[1], ave_wl_cal[1], bw[1], format='(a4,2e15.7)'
    printf, lun, pixel[2], ave_wl_cal[2], bw[2], format='(a4,2e15.7)'
    printf, lun, pixel[3], ave_wl_cal[3], bw[3], format='(a4,2e15.7)'
    printf, lun, pixel[4], ave_wl_cal[4], bw[4], format='(a4,2e15.7)'
    printf, lun, pixel[5], ave_wl_cal[5], bw[5], format='(a4,2e15.7)'
    printf, lun, ''
    printf, lun, 'Pixel       TF^2 AC     TF_err^2 AC           V2 AC       V2_err AC         TF^2 BD     TF_err^2 BD           V2 BD       V2_err BD'
    for j=0,5 do printf, lun, pixel[j], T2_AC_sci[j,i], T2_AC_err_sci[j,i], V2_AC_sci[j,i], V2_AC_err_sci[j,i], T2_BD_sci[j,i], T2_BD_err_sci[j,i], V2_BD_sci[j,i], V2_BD_err_sci[j,i], format='(a3,8f16.10)'


    close, lun
    free_lun, lun

  endfor

endif

;------------------------------------------------------------------------------



;oiFitsUtils.i 6600, /src/yorick/i




stop
end
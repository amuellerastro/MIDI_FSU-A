@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@fsc_color.pro
@mima.pro
@legend.pro
@stddev.pro
@moment.pro
@showsym.pro
@sigfig.pro
@amedian.pro
@reverse.pro
@numlines.pro
@repchr.pro
@is_ieee_big.pro
@ieee_to_host.pro
@counter.pro
@fifteenb.pro

; function wrap_data, data
;   return, (((data mod 360.d0) + 180.d0) mod 360.d0) - 180.d0
; end

pro plot_MIDI_FSUA_MeasuredPredicted_CustomMask

;FOR CALIBRATORS
;**********************************************************************************************

commnum = ''
read, 'Enter number of Commissioning: ', commnum

readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', skipline=1, /silent

idx = where(comm eq commnum)
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

idxsc = where(scical eq 'Cal')
 if (idxsc[0] ne -1) then begin
  comm = comm[idxsc]
  id = id[idxsc]
  scical = scical[idxsc]
  nfl = nfl[idxsc]
  midi = midi[idxsc]
  fsu = fsu[idxsc]
  maskname = maskname[idxsc]
  n = n_elements(idxsc)
endif else begin
  print, 'No calibrator targets observed in that run.'
  stop
endelse

nobs = n_elements(idxsc)

time = strmid(midi,16,24)
night = strmid(midi, 5, 10)

EWSdir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/MIDIreduced_EWS_CustomMask/'
Koreskodir = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/MIDIreduced_Koresko_CustomMask/'

resultsdir = Koreskodir+'GD_DISP_KoreskoPredictedMeasured/'
file_mkdir, resultsdir


;**********************************************************************************************

; for i=0,nobs-1 do begin
for i=14,14 do begin

  counter, i+1, nobs, 'Observation '

  ;MIDI readin data, Disp. and GD
  st3 = mrdfits(EWSdir+id[i]+'_'+time[i]+'.groupdelay.fits', 3, /silent)
  midigd = (st3.delay)*299792458.d0*1.d6	;micrometer
  miditime = st3.time
  miditime = (miditime-miditime[0])*86400.d0

  st3 = mrdfits(EWSdir+id[i]+'_'+time[i]+'.ungroupdelay.fits', 3, /silent)
  mididisp = (st3.dispersion)	;deg


  ;FSU read in data
;   readcol, Koreskodir+id[i]+'_'+time[i]+'.KoreskoDISP', fsudisp_uw, format='d', /silent	;in deg
;   readcol, Koreskodir+id[i]+'_'+time[i]+'.KoreskoGD', fsugd, format='d', /silent	;in deg
;   readcol, Koreskodir+id[i]+'_'+time[i]+'.GDoffset', fsugdoff, format='d', /silent	;in deg
;   fsugd = (fsugd*299792458.d0)
;   fsugd = (fsugd-median(fsugd)+fsugdoff[0])*1.d6
  st3 = mrdfits(Koreskodir+id[i]+'_'+time[i]+'.groupdelay.fits', 3, /silent)
  fsugd = (st3.delay)*299792458.d0*1.d6	;micrometer
  st3 = mrdfits(Koreskodir+id[i]+'_'+time[i]+'.dispersion.fits', 3, /silent)
  fsudisp = st3.dispersion

  readcol, Koreskodir+id[i]+'_'+time[i]+'.fsuflags', fsuflag, format='d', /silent

  idx = where(fsuflag eq 1)
  midigd = midigd[idx]
  miditime = miditime[idx]
  mididisp = mididisp[idx]
  fsugd = fsugd[idx]
  fsudisp = fsudisp[idx]

  ;fsudisp = wrap_data(fsudisp_uw)

  ;convert DISP / deg to rad
  mididisp = mididisp*!dtor
  fsudisp = fsudisp*!dtor

  ;gdres = midigd-fsugd
  ;sdgdres = stddev(gdres)

  medmidigd = amedian(midigd, 51)
  gdres = medmidigd-fsugd
  sdgdres = stddev(gdres)


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
;     device, isolatin = 1
    device, filename=resultsdir+id[i]+'_'+time[i]+'_KoreskoPredMeas.ps', /color,XSIZE=23, YSIZE=17, XOffset=xoffset, YOffset=yoffset

    !p.multi=[0,1,3]
;     !p.font=0
    !p.thick=4
    !x.thick=3
    !y.thick=3

    micron = '!Mm!X'+'m'
    sigma = '!Ms!X'

    ;dispersion
    plot, /nodata, miditime, mididisp, yr=[-!DPI,!DPI], yst=1, xst=1, $	;xr=[0,max(miditime)]
      ytitle='Dispersion!C[rad]', charsize=2.5, xtickformat="(A1)", $
      color=fsc_color('black'), background=fsc_color('white'), yminor=2, $
      pos=[0.12,0.63,0.97,0.94], xminor=5, xticklen=0.08, $
      title=id[i]+' / '+night[i]+' / UTC: '+time[i]
    oplot, miditime, mididisp, color=fsc_color('blue'), psym=3
    oplot, miditime, fsudisp, color=fsc_color('red'), psym=3

    ;GD
    yr = mima([medmidigd, fsugd])
;     yr = yr*1.1
    plot, /nodata, miditime, midigd, xst=1, yr=yr, $
      ytitle='Group Delay!C['+micron+']', charsize=2.5, xtickformat="(A1)", $
      color=fsc_color('black'), background=fsc_color('white'), yminor=2, $
      pos=[0.12,0.31,0.97,0.62], xminor=5, xticklen=0.08
    oplot, miditime, midigd, color=fsc_color('dark gray'), psym=3
    oplot, miditime, medmidigd, color=fsc_color('blue'), psym=3
    oplot, miditime, fsugd, color=fsc_color('red'), psym=3
    legend, ['MIDI measured', 'FSUA predicted'], psym=[sym(1),sym(1)], color=[fsc_color('blue'), fsc_color('red')], box=0, margin=0, /right, charsize=1.3, symsize=[0.9,0.9]


    ;GD residuals
    plot, /nodata, miditime, gdres, yr=[-5,5], yst=1, xst=1, $
      xtitle='Time [second]', ytitle='GD Residuals!C['+micron+']', charsize=2.5, $
      color=fsc_color('black'), background=fsc_color('white'), yminor=1, yticks=4, $
      pos=[0.12,0.12,0.97,0.30], xminor=5, xticklen=0.1
    oplot, miditime, gdres, color=fsc_color('green'), psym=3
;    plots, !x.crange, [0,0], linestyle=1, color=fsc_color('black')
    plots, !x.crange, [1,1], linestyle=1, color=fsc_color('black')
    plots, !x.crange, [-1,-1], linestyle=1, color=fsc_color('black')
    xyouts, 0.76*!x.crange[1], 0.6*!y.crange[1], sigma+'!DGD!N: '+sigfig(sdgdres, 3)+' '+micron, charsize=1.5

    device,/close
    set_plot,'x'

    ;spawn, 'gv '+resultsdir+id[i]+'_'+time[i]+'_KoreskoPredMeas.ps'


!p.multi=[0,1,0]
!P.Font=0
!p.thick=1
!x.thick=1
!y.thick=1

; where(midigd gt median(midigd)+5.*stddev(midigd))
; where(midigd lt median(midigd)-5.*stddev(midigd))



endfor

!p.multi=[0,1,0]
!P.Font=0
!p.thick=1
!x.thick=1
!y.thick=1

stop
end
@sxaddhist.pro
@wrapphase.pro
@get_eso_keyword.pro
@sync_FSUAMIDI_Koresko.pro
@mwrfits.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@readfits.pro
@sxpar.pro
@mean.pro
@moment.pro
@fxaddpar.pro
@detabify.pro
@fxparpos.pro
@ts_diff.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@mpfitfun.pro
@mpfit.pro
@sxaddpar.pro
@counter.pro
@fifteenb.pro
@amedian.pro
@reverse.pro
@mpfitexpr.pro

; function wrap_data, data
;   return, (((data mod 360.d0) + 360.d0) mod 360.d0) - 180.d0
; end

function func_EWS_Koresko_GD, x, p, _extra=struct

  refgd = struct.midigd
  fit = x+p[0];-refgd

;   window, 0
;   plot, refgd
;   oplot, fit, color=rgb(255,0,0)

  return, fit

end



;reduces MIDI data taken with FSU-A as fringe sensor

pro reduceMIDI_Koresko_CustomMask


;**********************************************************************************************


;ewspath = '/opt/MIA+EWS-2009Dec02/c/bin/'
;ewspath = '/opt/MIA+EWS-2011Sep08/c/bin/'
; ewspath = '/opt/MIA+EWS-2011Oct22/c/bin/'
;ewsversion = 'MIA+EWS-2011Dec13'
; ewsversion = 'MIA+EWS-2013Feb02'
; ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'

obsnum = ''
read, 'Enter Observing number: ', obsnum


readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', skipline=1, /silent

idx = where(comm eq obsnum)
if (idx[0] ne -1) then begin

  comm = comm[idx]
  id = id[idx]
  scical = scical[idx]
  nfl = nfl[idx]
  midi = midi[idx]
  fsu = fsu[idx]
  maskname = maskname[idx]
  photA = photA[idx]
  photB = photB[idx]

endif else begin

  print, ''
  print, 'Observing number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

;path = '/opt/MIDI_FSUA/MIDI_FSUA_COMM'+obsnum+'/'	;masks are there
midipath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIdata/'
fsupath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/FSUAdata/'
resultsdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_Koresko_CustomMask/'
maskpath = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDI_CustomMasks/'
midiewspath  = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_EWS_CustomMask/'
file_mkdir, resultsdir
workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir

;**********************************************************************************************

prime = "'"
time = strmid(midi,16,24)

;do we want to reduce calibrator or science observations?
quest = ''
print, ''
read, 'Calibrator (c) or Science (s)?  ',quest 
print, ''

if (quest ne 's' and quest ne 'c') then begin

  print, ''
  print, 'Use one of the options.'
  stop

endif else begin

  if (quest eq 's') then begin	;SCIENCE targets
    idxsc = where(scical eq 'Sci')
    if (idxsc[0] ne -1) then begin
      comm = comm[idxsc]
      id = id[idxsc]
      scical = scical[idxsc]
      nfl = nfl[idxsc]
      midi = midi[idxsc]
      fsu = fsu[idxsc]
      maskname = maskname[idxsc]
      time = time[idxsc]
      n = n_elements(idxsc)
    endif else begin
      print, 'No science targets observed in that run.'
      stop
    endelse

  endif else begin	;CALIBRATOR targets

    idxsc = where(scical eq 'Cal')
    if (idxsc[0] ne -1) then begin
      comm = comm[idxsc]
      id = id[idxsc]
      scical = scical[idxsc]
      nfl = nfl[idxsc]
      midi = midi[idxsc]
      fsu = fsu[idxsc]
      maskname = maskname[idxsc]
      time = time[idxsc]
      n = n_elements(idxsc)
    endif else begin
      print, 'No calibrator targets observed in that run.'
      stop
    endelse

  endelse

endelse

;n = n_elements(id)

for i=0,n-1 do begin
;for i=1,1 do begin
;for i=15,15 do begin
;HD1014 03:30:57
;stop

;  counter, i+1, n, 'Processing observation '
  print, ''
  print, 'Observation '+strcompress(i+1,/rem)+' / '+strcompress(n,/rem)


  if (quest eq 's') then begin
    ;print, id[i], '  ', midi[i]
    ;get median GD position from Calibrator derived by standard EWS reduction
    readcol, resultsdir+maskname[i]+'.GDoffset', med_gd_pos, format='d', /silent
    med_gd_pos = double(med_gd_pos[0])
    ;med_gd_pos = dblarr(1)
    print, 'Using median GD position (meter) of the calibrator '+maskname[i]+': ', med_gd_pos
    print, ''
  endif


;Mask selection
  mask = maskpath+maskname[i]+'.srcmask.fits'
;   if (scical[i] eq 'Cal') then mask = maskpath+id[i]+'_'+time[i]+'.srcmask.fits'
;   if (scical[i] eq 'Sci') then begin
;     filemask = file_search(maskpath+'*.srcmask.fits', count=nmasks)
;     for j=0,nmasks-1 do print, j+1, '  ', filemask[j]
;     print, ''
;     read, 'Choose a mask from a calibrator for '+id[i]+'_'+time[i]+': ', questmask
;     print, ''
;     mask = filemask[questmask-1]
;   endif

  file = file_search(midipath+midi[i]+'*',count=nfiles)

  check99 = strmid(file[nfiles-1],strlen(file[nfiles-1])-7,2)
  if (check99 eq '99') then begin
    file = file[0:nfiles-2]
    nfiles = n_elements(file)
  endif

;--------
  if (nfiles eq 1) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 2) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 3) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 4) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 5) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 6) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 7) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 8) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 9) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 10) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 11) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 12) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'
  if (nfiles eq 13) then spawn, ewspath+'oir1dCompressData '+prime+file[0]+' '+file[1]+' '+file[2]+' '+file[3]+' '+file[4]+' '+file[5]+' '+file[6]+' '+file[7]+' '+file[8]+' '+file[9]+' '+file[10]+' '+file[11]+' '+file[12]+prime+' '+mask+' '+workingdir+id[i]+'_'+time[i]+'.compressed.fits'

  spawn, ewspath+'oirFormFringes '+workingdir+id[i]+'_'+time[i]+'.compressed.fits '+workingdir+id[i]+'_'+time[i]+'.fringes.fits'
;   if (nfl[i] le 7.) then begin
;     spawn, ewspath+'oirFormFringes '+id[i]+'_'+time[i]+'.compressed.fits '+id[i]+'_'+time[i]+'.fringes.fits -removeAverage'
;   endif else begin
;     spawn, ewspath+'oirFormFringes '+id[i]+'_'+time[i]+'.compressed.fits '+id[i]+'_'+time[i]+'.fringes.fits'
;   endelse

  spawn, ewspath+'oirRotateInsOpd '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.insopd.fits'
;   ;execute this to create a .groupdelay fits file
  spawn, ewspath+'oirGroupDelay '+workingdir+id[i]+'_'+time[i]+'.insopd.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.dummy.fits'


  ;replace MIDI GD by KoreskoGD
  openr, lun, resultsdir+id[i]+'_'+time[i]+'.KoreskoGD_orig', /get_lun
  rows = file_lines(resultsdir+id[i]+'_'+time[i]+'.KoreskoGD_orig')
  data = dblarr(1,rows)
  readf, lun, data
  tmp = data[0,*]
  kgd = reform(tmp)	;in meter, contains NAN values
  close, lun
  free_lun, lun

  openr, lun, resultsdir+id[i]+'_'+time[i]+'.KoreskoDISP_orig', /get_lun
  rows = file_lines(resultsdir+id[i]+'_'+time[i]+'.KoreskoDISP_orig')
  data = dblarr(1,rows)
  readf, lun, data
  tmp = data[0,*]
  kdisp = reform(tmp)	;in meter, contains NAN values
  close, lun
  free_lun, lun

  ;convert to seconds
  kgd_s = kgd/299792458.d0

;sync with MIDI time
  sync_FSUAMIDI_Koresko, id[i]+'_'+time[i],workingdir+id[i]+'_'+time[i]+'.compressed.fits',fsupath+fsu[i], kgd_s, kdisp, obsnum, workingdir

  ;read KoreskoGD [seconds] in and write into existent groupdelay file
  readcol, workingdir+id[i]+'_'+time[i]+'.KoreskoGD', kgd_sync, format='d', /silent

;fit GD  and Dispersion offset for bright Calibrators
  if (quest eq 'c') then begin

  ;fit GD
    kgd_sync = kgd_sync*299792458.d0	;in meter
    kgd = kgd_sync
    ;st = mrdfits(midiewspath+id[i]+'_'+time[i]+'.ORIGINALgroupdelay.fits',3, /silent)
    st = mrdfits(midiewspath+id[i]+'_'+time[i]+'.groupdelay.fits',3, /silent)
    midigd = st.delay*299792458.d0

    readcol, workingdir+id[i]+'_'+time[i]+'.fsuflags', fsuflag, format='i', /silent
    idxgood = where(fsuflag eq 1)
    midigd = midigd[idxgood]
    med_kgd = median(kgd[idxgood])
    kgd = kgd[idxgood]-med_kgd

    ;fit offset between EWS and Koresko GD
    xdum = dindgen(n_elements(midigd))
    dumerr = dblarr(n_elements(midigd))
    dumerr[*] = 1.0
    struct = {midigd:midigd}
    start_val = median(midigd)
    maxiter = 1000

    ;cut out NaN if there are any
    idx = where(finite(kgd) eq 1)
    kgd = kgd[idx]
    midigd = midigd[idx]

    ;apply the fit on a slightly median filtered MIDI GD
    med_midigd = amedian(midigd, 51)
    ;med_midigd = midigd

    med_gd_pos = MPFITfun('func_EWS_Koresko_GD', kgd, med_midigd, dumerr, start_val, weights = 1.d0, $
	    maxiter=maxiter, niter=niter, status=status, bestnorm=bestnorm, $
	    perror=perror, covar=covar, yfit=yfit, dof=dof, functargs=struct, /quiet);weights = rverr

    kgd_sync[idxgood] = ((kgd_sync[idxgood]-med_kgd)+(med_gd_pos[0]))/299792458.d0	;in seconds to write it in groupdelay.fits

    openw, lun, workingdir+id[i]+'_'+time[i]+'.GDoffset', /get_lun
      printf, lun, med_gd_pos[0], format='(e17.8)'
    close, lun
    free_lun, lun

;   ;fit dispersion
;     readcol, workingdir+id[i]+'_'+time[i]+'.KoreskoDISP', kdisp_sync, format='d', /silent
;     kdisp_sync_tmp = kdisp_sync[idxgood]	;Koresko Disp
;     st3 = mrdfits(midiewspath+id[i]+'_'+time[i]+'.ungroupdelay.fits', 3, /silent)	;MIDI EWS Disp
;     mididisp = st3.dispersion
;     mididisp_tmp = mididisp[idxgood]
;     ;both Dispersions are already in degrees, Koresko Disp is NOT wrapped
; 
;     expr = 'WrapPosNeg180(P[0]*X+ P[1])'
;     start_val = [1.d0,100.d0]
;     dumerr = dblarr(n_elements(mididisp_tmp))
;     dumerr[*] = 1d
;     maxiter = 1000
; 
;     params = mpfitexpr(expr, kdisp_sync_tmp, mididisp_tmp, dumerr, start_val, weights=dumerr, maxiter=maxiter, /quiet)
; 
;     kdisp_sync_final = WrapPosNeg180(params[0]*kdisp_sync+params[1])
; 
;     openw, lun, workingdir+id[i]+'_'+time[i]+'.DISPoffset', /get_lun
;       printf, lun, params[0], params[1], format='(2f20.13)'
;     close, lun
;     free_lun, lun

  endif

  if (quest eq 's') then  begin

    readcol, workingdir+id[i]+'_'+time[i]+'.fsuflags', fsuflag, format='d', /silent
    idxgood = where(fsuflag eq 1)
    med_kgd_sync = median(kgd_sync[idxgood])
    kgd_sync[idxgood] = (kgd_sync[idxgood]-med_kgd_sync)+(med_gd_pos[0]/299792458.d0)

  endif


  for j=0,3 do begin

    data = mrdfits(workingdir+id[i]+'_'+time[i]+'.groupdelay.dummy.fits', j, hdr,/silent)
    if (j ne 3) then begin
      mwrfits, data, workingdir+id[i]+'_'+time[i]+'.groupdelay.fits', hdr,/silent
    endif else begin
      data.delay = kgd_sync
      mwrfits, data, workingdir+id[i]+'_'+time[i]+'.groupdelay.fits', hdr, /silent
    endelse

  endfor

  spawn, 'rm '+workingdir+id[i]+'_'+time[i]+'.groupdelay.dummy.fits'

    spawn, ewspath+'oirRotateGroupDelay '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.dummy.fits'; -ngrad 3 -smooth 0.36'


;replace Dispersion table of ungroupdelay.fits with predicted phase delay by Koresko

  ;read KoreskoPD [deg] in and write into existent ungroupdelay file
;   if (quest eq 's') then readcol, workingdir+id[i]+'_'+time[i]+'.KoreskoDISP', kdisp_sync, format='d', /silent else kdisp_sync = kdisp_sync_final
  readcol, workingdir+id[i]+'_'+time[i]+'.KoreskoDISP', kdisp_sync, format='d', /silent

  for j=0,3 do begin
    data = mrdfits(workingdir+id[i]+'_'+time[i]+'.ungroupdelay.dummy.fits', j, hdr,/silent)
    if (j ne 3) then begin
	mwrfits, data, workingdir+id[i]+'_'+time[i]+'.dispersion.fits', hdr,/silent
    endif else begin
      if (quest eq 's') then kdisp_sync = WrapPosNeg180(kdisp_sync)
      ;kdisp_sync = wrap_data(kdisp_sync)
      data.dispersion = kdisp_sync
      mwrfits, data, workingdir+id[i]+'_'+time[i]+'.dispersion.fits', hdr, /silent
    endelse
  endfor
  spawn, ewspath+'oirRotateGroupDelay '+workingdir+id[i]+'_'+time[i]+'.fringes.fits '+workingdir+id[i]+'_'+time[i]+'.groupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits -phase '+workingdir+id[i]+'_'+time[i]+'.dispersion.fits'

  spawn, 'rm '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.dummy.fits'


;create artificial .flag file
  st1 = mrdfits(midiewspath+id[i]+'_'+time[i]+'.flag.fits', 0, hdr1, /silent)
  st2 = mrdfits(midiewspath+id[i]+'_'+time[i]+'.flag.fits', 1, hdr2, /silent)

  ntmp = n_elements(st2.timerang[0])


  st3 = mrdfits(workingdir+id[i]+'_'+time[i]+'.groupdelay.fits', 3, hdr3,/silent)	;need time stamps in JD
  time_tmp1 = st3.time

  readcol, workingdir+id[i]+'_'+time[i]+'.fsuflags', fsuflag, format='d', /silent

  idx0 = where(fsuflag eq 0.)
  ngood = n_elements(where(fsuflag eq 1.))

  time_tmp2 = time_tmp1[idx0]

  time_tmp3 = dblarr(2,n_elements(idx0))
  time_tmp3[0,0:n_elements(idx0)-1] = time_tmp2
  time_tmp3[1,0:n_elements(idx0)-1] = time_tmp2


  nrows = n_elements(time_tmp3[0,*])

  st3 = replicate({timerang:[0.d0,0.d0], source_id:st2[0].source_id, telescope:[st2[0].telescope[0],st2[0].telescope[1]], chans:[st2[0].chans[0],st2[0].chans[1]], reason:''},nrows)

  st3[*].timerang[0] = transpose(time_tmp3[0,*])
  st3[*].timerang[1] = transpose(time_tmp3[1,*])

  reason_tmp = strarr(nrows)
  reason_tmp[*] = 'JUMP IN GROUP DELAY'
;  reason_tmp[n_elements(reason_tmp)-1] = 'END of OBS'

  st3[*].reason = reason_tmp

  ;change number of NGOOD in header
  sxaddpar, hdr1, 'NGOOD', strcompress(ngood, /rem)

  sxaddpar, hdr2, 'EXTNAME', 'FLAG'

  mwrfits, st1, workingdir+id[i]+'_'+time[i]+'.flag.fits', hdr1, /silent
  mwrfits, st3, workingdir+id[i]+'_'+time[i]+'.flag.fits', hdr2, /silent

  spawn, ewspath+'oirAverageVis '+workingdir+id[i]+'_'+time[i]+'.ungroupdelay.fits '+workingdir+id[i]+'_'+time[i]+'.flag.fits '+workingdir+id[i]+'_'+time[i]+'.corr.fits'


  spawn, 'mv '+workingdir+id[i]+'* '+resultsdir+'/'
;stop
endfor

spawn, 'rm -r '+workingdir

;stop
end

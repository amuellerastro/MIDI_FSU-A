@stddev.pro
@reverse.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@uniq.pro
@mrdfits.pro
@fxposit.pro
@fxmove.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@linspace.pro
@mean.pro
@moment.pro
@fsc_color.pro
@oploterror.pro
@setdefaultvalue.pro
@tag_exist.pro
@cgquery.pro
@cgplot.pro
@setdecomposedstate.pro
@decomposedcolor.pro
@cgdefaultcolor.pro
@getdecomposedstate.pro
@colorsareidentical.pro
@cgcolor.pro
@cgdefcharsize.pro
@mwrfits.pro
@fxaddpar.pro
@detabify.pro
@sxaddhist.pro
@fxparpos.pro



pro avergae_faint_photometry_CustomMask

;------------------------------------

win = 0.2	;average of wl+/-win

;------------------------------------

; ewsversion = 'MIA+EWS-2013May18'
;ewsversion = 'MIA+EWS-2013Nov19'
ewsversion = 'MIA+EWS-2014Feb16'

ewspath = '/home/amueller/src/'+ewsversion+'/c/bin/'

obsnum = ''
obsnum = 'P92CHOQUET'
; read, 'Enter Observation ID: ', obsnum


readcol, 'observation.txt', comm, id, scical, nfl, midi, photA, photB, fsu, maskname, photo, format='a,a,a,d,a,a,a,a,a,a', skipline=1, /silent

idx = where(comm eq obsnum)
if (idx[0] ne -1) then begin

  comm = comm[idx]
  id = id[idx]
  scical = scical[idx]
  nfl = nfl[idx]
  midi = midi[idx]
  fsu = fsu[idx]
  maskname = maskname[idx]
  photo = photo[idx]
  photA = photA[idx]
  photB = photB[idx]

endif else begin

  print, ''
  print, 'Commissioning number does not match with current data set in observation.txt.'
  print, ''
  return

endelse

photdir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/Photometry_CustomMask/'
dir = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obsnum+'/MIDIreduced_Koresko_CustomMask/'
workingdir = '/home/amueller/work/MIDI_FSUA/temp_WorkingDir/'
file_mkdir, workingdir

time = strmid(midi,16,24)

idsingle = id[sort(id)]
idsingle = idsingle[uniq(idsingle)]
nstars = n_elements(idsingle)

p = dblarr(n_elements(id), 6, 171)	;masked photometry of each obs
pe = dblarr(n_elements(id),6, 171)	;corresponding error
vis = dblarr(n_elements(id), 171)
vise = vis
phi = vis
phie = vis

; avep = dblarr(nstars, 171)	;average photometry
; avepe = dblarr(nstars, 171)	;error

;================================================================================

;read in data

for xx=0,n_elements(id)-1 do begin

  file = file_search(photdir+id[xx]+'_'+time[xx]+'.photometry.fits', count=nfiles)
    st = mrdfits(file, 2, /silent)
    p[xx,0,*] = st[0].data1
    p[xx,1,*] = st[1].data1
    p[xx,2,*] = st[2].data1
    p[xx,3,*] = st[3].data1
    p[xx,4,*] = st[4].data1
    p[xx,5,*] = st[5].data1
    pe[xx,0,*] = st[6].data1
    pe[xx,1,*] = st[7].data1
    pe[xx,2,*] = st[8].data1
    pe[xx,3,*] = st[9].data1
    pe[xx,4,*] = st[10].data1
    pe[xx,5,*] = st[11].data1

  ;for wavelength only
  if (xx eq 0) then begin

    st = mrdfits(dir+id[xx]+'_'+time[xx]+'.corr.fits', 1, /silent)
    wl = st.eff_wave*1.d6

  endif

  st = mrdfits(dir+id[xx]+'_'+time[xx]+'.corr.fits', 3, /silent)
  vis[xx,*] = st.visamp
  vise[xx,*] = st.visamperr
  phi[xx,*] = st.visphi
  phie[xx,*] = st.visphierr

endfor

;================================================================================

;define averaging inside an observation

n = round((13.5-7.)/win)
avewl = linspace(7., 13.5, n)

mp = dblarr(n_elements(id), 6, n_elements(avewl))	;binned photometry
mpe = mp
visamp = dblarr(n_elements(id), n_elements(avewl))
visamperr = dblarr(n_elements(id), n_elements(avewl))
visphi = dblarr(n_elements(id), n_elements(avewl))
visphierr = dblarr(n_elements(id), n_elements(avewl))
nbin = dblarr(n_elements(id), n)	;number of elements in one bin for each spectrum
sigmashort = dblarr(n_elements(id), 6, n_elements(avewl))

for xx=0,n_elements(id)-1 do begin

  for i=0,n_elements(avewl)-1 do begin

    idx = where(wl ge avewl[i]-win and wl le avewl[i]+win)

    nbin[xx, i] = n_elements(idx)

    mp[xx,0,i] = mean(p[xx,0,idx])
    mp[xx,1,i] = mean(p[xx,1,idx])
    mp[xx,2,i] = mean(p[xx,2,idx])
    mp[xx,3,i] = mean(p[xx,3,idx])
    mp[xx,4,i] = mean(p[xx,4,idx])
    mp[xx,5,i] = mean(p[xx,5,idx])

    mpe[xx,0,i] = mean(pe[xx,0,idx])
    mpe[xx,1,i] = mean(pe[xx,1,idx])
    mpe[xx,2,i] = mean(pe[xx,2,idx])
    mpe[xx,3,i] = mean(pe[xx,3,idx])
    mpe[xx,4,i] = mean(pe[xx,4,idx])
    mpe[xx,5,i] = mean(pe[xx,5,idx])


    sigmashort[xx,0,i] = mpe[xx,0,i]
    sigmashort[xx,1,i] = mpe[xx,1,i]
    sigmashort[xx,2,i] = mpe[xx,2,i]
    sigmashort[xx,3,i] = mpe[xx,3,i]
    sigmashort[xx,4,i] = mpe[xx,4,i]
    sigmashort[xx,5,i] = mpe[xx,5,i]


    ;weighted mean for Fcorr and its error

;     visamp[xx,i] = mean(vis[xx,idx])
    visamp[xx,i] = (total((1./vise[xx,idx]^2.)*vis[xx,idx]))/total(1./vise[xx,idx]^2.)	;weighted mean
    visamperr[xx,i] = sqrt(n_elements(vis[xx,idx])/total(1./vise[xx,idx]^2.))

    visphi[xx,i] = (total((1./phie[xx,idx]^2.)*phi[xx,idx]))/total(1./phie[xx,idx]^2.)	;weighted mean
    visphierr[xx,i] = sqrt(n_elements(phi[xx,idx])/total(1./phie[xx,idx]^2.))

;     visphi[xx,i] = mean(phi[xx,idx])
;     tmp = (moment(phie[xx,idx]))[1]
;     tmp = tmp;/(nbin[xx,i])
;     visphierr[xx,i] = tmp


  endfor

endfor

;================================================================================

;compute sigma long

sigmalong = dblarr(nstars, 6, n_elements(avewl))

for xx=0,nstars-1 do begin

  idx = where(id eq idsingle[xx])

  ;compute sigma long
  for i=0,n_elements(avewl)-1 do begin

;     var = sigmashort[idx,0,i]
;     sigmalong[xx,0,i] = mean(var)-mean(sqrt(var))^2.
;     var = sigmashort[idx,1,i]
;     sigmalong[xx,1,i] = mean(var)-mean(sqrt(var))^2.
;     var = sigmashort[idx,2,i]
;     sigmalong[xx,2,i] = mean(var)-mean(sqrt(var))^2.
;     var = sigmashort[idx,3,i]
;     sigmalong[xx,3,i] = mean(var)-mean(sqrt(var))^2.
;     var = sigmashort[idx,4,i]
;     sigmalong[xx,4,i] = mean(var)-mean(sqrt(var))^2.
;     var = sigmashort[idx,5,i]
;     sigmalong[xx,5,i] = mean(var)-mean(sqrt(var))^2.

    sigmalong[xx,0,i] = stddev(mp[idx,0,i])
    sigmalong[xx,1,i] = stddev(mp[idx,1,i])
    sigmalong[xx,2,i] = stddev(mp[idx,2,i])
    sigmalong[xx,3,i] = stddev(mp[idx,3,i])
    sigmalong[xx,4,i] = stddev(mp[idx,4,i])
    sigmalong[xx,5,i] = stddev(mp[idx,5,i])


  endfor


endfor

;================================================================================

;compute final error, Eq. 5.3 in Burtscher, 2011

mphot = dblarr(nstars, 6, n)
mphoterr = dblarr(nstars, 6, n)

for xx=0,nstars-1 do begin

  idx = where(id eq idsingle[xx])

  for i=0,n_elements(avewl)-1 do begin

    mphot[xx,0,i] = mean(mp[idx,0,i])
    mphoterr[xx,0,i] = sqrt(((mean(sigmashort[idx,0,i]^2.)/nbin[0,i]) + sigmalong[xx,0,i]^2.)/double(n_elements(idx)))
    mphot[xx,1,i] = mean(mp[idx,1,i])
    mphoterr[xx,1,i] = sqrt(((mean(sigmashort[idx,1,i]^2.)/nbin[0,i]) + sigmalong[xx,1,i]^2.)/double(n_elements(idx)))
    mphot[xx,2,i] = mean(mp[idx,2,i])
    mphoterr[xx,2,i] = sqrt(((mean(sigmashort[idx,2,i]^2.)/nbin[0,i]) + sigmalong[xx,2,i]^2.)/double(n_elements(idx)))
    mphot[xx,3,i] = mean(mp[idx,3,i])
    mphoterr[xx,3,i] = sqrt(((mean(sigmashort[idx,3,i]^2.)/nbin[0,i]) + sigmalong[xx,3,i]^2.)/double(n_elements(idx)))
    mphot[xx,4,i] = mean(mp[idx,4,i])
    mphoterr[xx,4,i] = sqrt(((mean(sigmashort[idx,4,i]^2.)/nbin[0,i]) + sigmalong[xx,4,i]^2.)/double(n_elements(idx)))
    mphot[xx,5,i] = mean(mp[idx,5,i])
    mphoterr[xx,5,i] = sqrt(((mean(sigmashort[idx,5,i]^2.)/nbin[0,i]) + sigmalong[xx,5,i]^2.)/double(n_elements(idx)))

  endfor

endfor


mwl = avewl
; mphot = mphot
; mphoterr = mphoterr

;================================================================================


micron = '!Mm!X'+'m'

;output masked A+B flux only

for i=0,nstars-1 do begin

  set_plot, 'ps'
  device, isolatin=1
  device, filename=photdir+idsingle[i]+'_AveragePhotometry.ps', /color,XSIZE=25, YSIZE=15, XOffset=xoffset, YOffset=yoffset
  !p.font=0
  ; !p.multi=[0,1,nstars]

  idx = where(idsingle[i] eq id)

  yrmax = max(p[idx,5,*])
  plot, wl, p[idx[0],5,*], /nodata, charsize=1.5, xr=[8,13], xst=1, yr=[0,yrmax], title=idsingle[i], $
    xtitle='Wavelength / ['+micron+']', ytitle='Masked Flux Geometric Mean of A&B / [ADU/s]'
  for j=0,n_elements(idx)-1 do begin

    oplot, wl, p[idx[j],5,*]

  endfor
  oploterror, mwl, mphot[i,5,*], mphoterr[i,5,*], color=fsc_color('red'), /nohat, thick=3, errthick=3


  device,/close
  set_plot,'x'

;   spawn, 'gv '+photdir+idsingle[i]+'_AveragePhotometry.ps'


  openw, lun, photdir+idsingle[i]+'_AveragePhotometry.txt', width=1400, /get_lun

    printf, lun, 'Wavelength   Average Flux    Error'
    for j=0,n_elements(mwl)-1 do printf, lun, mwl[j], mphot[i,5,j], mphoterr[i,5,j], format='(f10.6, 2f15.8)'

  close, lun
  free_lun, lun

endfor

!p.multi=[0,1,0]

;================================================================================

;photometry is now averaged and has less than 171 elements
;create a temporary photometry file as well as a temporary correlated flux file needed for oirRedCal to compute intrumental visibilities

for xx=0,nstars-1 do begin

 idx = where(id eq idsingle[xx])

  ;---------------------------------------------
  ;.photometry.fits

  for i=0,n_elements(idx)-1 do begin

    for j=0,2 do begin

      st = mrdfits(photdir+id[idx[i]]+'_'+time[idx[i]]+'.photometry.fits', j, hdr, /silent)

      if (j ne 2) then begin

	if (j eq 0) then begin	;change header keyword

	  if (strmatch(hdr[98], '*HIERARCH ESO DET WIN1 NX*') eq 1) then hdr[98] = 'HIERARCH ESO DET WIN1 NX     =          '+strcompress(uint(n),/rem)+' / # of pixels along X' else stop

	  if (strmatch(hdr[102], '*HIERARCH ESO DET WIN2 NX*') eq 1) then hdr[102] = 'HIERARCH ESO DET WIN2 NX     =          '+strcompress(uint(n),/rem)+' / # of pixels along X' else stop

	endif

      if (j eq 1) then begin

	st.naxis[0] = uint(n)

      endif

	mwrfits, st, workingdir+id[idx[i]]+'_'+time[idx[i]]+'.photometry.fits', hdr, /silent

      endif else begin

	if (strmatch(hdr[51], '*TDIM11*') eq 1) then hdr[51] = 'TDIM11  = ('+strcompress(uint(n),/rem)+',1)            / Column dimensionality' else stop

	prime = "'"
 	if (strmatch(hdr[35], '*TFORM11*') eq 1) then hdr[35] = 'TFORM11 = '+prime+strcompress(uint(n), /rem)+'E'+prime+'                / data format of field: 4-byte REAL' else stop

	data1 = fltarr(n)
	tmp = {frame:st[0].frame,time:st[0].time,exptime:st[0].exptime,opt_train:st[0].opt_train,reference:st[0].reference,opd:st[0].opd,localopd:st[0].localopd,offset:st[0].offset,rotation:st[0].rotation,stepping_phase:st[0].stepping_phase,data1:data1,target1:st[0].target1,tartyp1:st[0].tartyp1,ins_train:st[0].ins_train}
; 
	stnew = replicate(tmp, 12)

 	stnew[0].data1 = reverse(reform(mphot[xx,0,*]))
 	stnew[1].data1 = reverse(reform(mphot[xx,1,*]))
 	stnew[2].data1 = reverse(reform(mphot[xx,2,*]))
 	stnew[3].data1 = reverse(reform(mphot[xx,3,*]))
 	stnew[4].data1 = reverse(reform(mphot[xx,4,*]))
 	stnew[5].data1 = reverse(reform(mphot[xx,5,*]))
 	stnew[6].data1 = reverse(reform(mphoterr[xx,0,*]))
 	stnew[7].data1 = reverse(reform(mphoterr[xx,1,*]))
 	stnew[8].data1 = reverse(reform(mphoterr[xx,2,*]))
 	stnew[9].data1 = reverse(reform(mphoterr[xx,3,*]))
 	stnew[10].data1 = reverse(reform(mphoterr[xx,4,*]))
 	stnew[11].data1 = reverse(reform(mphoterr[xx,5,*]))

;  	stnew[0].data1 = reverse(reform(mphot[xx,5,*]))
;  	stnew[1].data1 = reverse(reform(mphot[xx,5,*]))
;  	stnew[2].data1 = reverse(reform(mphot[xx,5,*]))
;  	stnew[3].data1 = reverse(reform(mphot[xx,5,*]))
;  	stnew[4].data1 = reverse(reform(mphot[xx,5,*]))
;  	stnew[5].data1 = reverse(reform(mphot[xx,5,*]))
;  	stnew[6].data1 = reverse(reform(mphoterr[xx,5,*]))
;  	stnew[7].data1 = reverse(reform(mphoterr[xx,5,*]))
;  	stnew[8].data1 = reverse(reform(mphoterr[xx,5,*]))
;  	stnew[9].data1 = reverse(reform(mphoterr[xx,5,*]))
;  	stnew[10].data1 = reverse(reform(mphoterr[xx,5,*]))
;  	stnew[11].data1 = reverse(reform(mphoterr[xx,5,*]))

	mwrfits, stnew, workingdir+id[idx[i]]+'_'+time[idx[i]]+'.photometry.fits', hdr, /silent;, /lscale;, /no_types;, /silent
; 	writefits, workingdir+id[idx[i]]+'_'+time[idx[i]]+'.photometry.fits', stnew, hdr, /append

      endelse

    endfor

  endfor

  ;---------------------------------------------
  ;.corr.fits

  for i=0,n_elements(idx)-1 do begin

    for j=0,3 do begin

      st = mrdfits(dir+id[idx[i]]+'_'+time[idx[i]]+'.corr.fits', j, hdr, /silent)

      if (j ne 3) then begin

	if (j eq 0) then begin	;change header keyword

	  if (strmatch(hdr[98], '*HIERARCH ESO DET WIN1 NX*') eq 1) then hdr[98] = 'HIERARCH ESO DET WIN1 NX     =          '+strcompress(uint(n),/rem)+' / # of pixels along X' else stop

	  if (strmatch(hdr[102], '*HIERARCH ESO DET WIN2 NX*') eq 1) then hdr[102] = 'HIERARCH ESO DET WIN2 NX     =          '+strcompress(uint(n),/rem)+' / # of pixels along X' else stop

	endif

	if (j eq 1) then begin	;replace wavelength

	  tmp = st.eff_band
	  eff_band = tmp[0:n-1]
	  sttmp = {eff_wave:0., eff_band:0.}
	  st = replicate(sttmp, n)
	  st.eff_wave = reverse(mwl/1.d6)
	  st.eff_band = eff_band

	  if (strmatch(hdr[18], '*NWAVE*') eq 1) then hdr[18] = 'NWAVE   =                  '+strcompress(uint(n),/rem)+' / Number of Wavelengths in table' else stop

	endif

	mwrfits, st, workingdir+id[idx[i]]+'_'+time[idx[i]]+'.corr.fits', hdr, /silent

      endif else begin

; 	flag = st.flag
; 	flag = flag[0:n-1]
	flag = strarr(n)
	flag[*] = ('F')

	stnew = {target_id:st.target_id,time:st.time,mjd:st.mjd,int_time:st.int_time,visamp:reverse(reform(visamp[idx[i],*])),visamperr:reverse(reform(visamperr[idx[i],*])),visphi:reverse(reform(visphi[idx[i],*])),visphierr:reverse(reform(visphierr[idx[i],*])),ucoord:st.ucoord,vcoord:st.vcoord,sta_index:st.sta_index,flag:flag}

	if (strmatch(hdr[38], '*TFORM12*') eq 1) then hdr[38] = 'TFORM12 = '+strcompress(uint(n), /rem)+'L'+'               / data format of field: 1-byte LOGICAL' else stop

	mwrfits, stnew, workingdir+id[idx[i]]+'_'+time[idx[i]]+'.corr.fits', hdr, logical_cols=[12], /silent

      endelse

    endfor

  endfor

endfor

;================================================================================

;oirRedcal to get instrumental visibilities

for i=0,n_elements(id)-1 do begin

  spawn, ewspath+'oirRedCal '+workingdir+id[i]+'_'+time[i]

endfor


;================================================================================


; 
; window, 0
; x=mrdfits(workingdir+'HD76111_03:37:09.redcal.fits', 3)
; plot, x.visamp
; ploterror, x.visamp, x.visamperr
; 
; window, 1
; x=mrdfits(workingdir+'HD76111_03:34:51.redcal.fits', 3)
; plot, x.visamp
; ploterror, x.visamp, x.visamperr

;================================================================================

;move files to resultsdir and rename customized .corr and .photometry

bindir = dir+'unbinned_corr/'
file_mkdir, bindir
bindirphot = photdir+'unbinned_photometry/'
file_mkdir, bindirphot

file = file_search(dir+'*.corr.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+bindir+'.'

file = file_search(photdir+'*photometry.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+bindirphot+'.'


file = file_search(workingdir+'*redcal.fits', count=nfiles)
for i=0,nfiles-1 do spawn, 'mv '+file[i]+' '+dir+'.'

for i=0,n_elements(id)-1 do begin

  ;spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'.corr.fits '+workingdir+id[i]+'_'+time[i]+'.binned.corr.fits'
;   spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'.binned.corr.fits '+dir+'.'
  spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'.corr.fits '+dir+'.'
;   spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'.photometry.fits '+workingdir+id[i]+'_'+time[i]+'.binned.photometry.fits'
;   spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'.binned.photometry.fits '+photdir+'.'
  spawn, 'mv '+workingdir+id[i]+'_'+time[i]+'.photometry.fits '+photdir+'.'

endfor


spawn, 'rm -r '+workingdir

stop
end
;if the normalized scan has negative values / values close to zero there will be big artefacts and the fringe will be filtered out...

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
@mean.pro
@moment.pro
@rgb.pro
@linspace.pro
@c_correlate.pro
@sixlin.pro
@proceeding_text.pro
@resistant_mean.pro
@poly.pro
@stddev.pro

;=================================================================================================
;no improvement if used instead auf Gauss function
; function fit_snr, x, p
;   fit = p[7]*( exp(-0.5d0*((x-p[0])/p[1])^2.) + p[4]*exp(-0.5d0*((x-p[5])/p[6])^2.)*sin(p[2]*x-p[3]) )+p[8]
;   return, fit
; end

function fit_sinus, x, p

  fit = p[0]*sin((2.d0*!DPI*x/p[1])+p[2])+p[3]

  return, fit

end

;=================================================================================================

function fit_fringe, x, p

  ;p0: I0, p1: l0, p2: dl0, p3: doffset, p4: Ioffset, p5: Phase

  lcoh = (p[1]^2.-(p[2]^2.)/4.d0)/p[2]
  eta0 = 1.d0

  k0 = 2.d0*!DPI/p[1]
  sinc = (sin(!DPI*(x-p[3])/lcoh))/(!DPI*(x-p[3])/lcoh)
  idxnan = where(finite(sinc) ne 1)	;in case there is a NAN when fringe going through 0 OPD
  if (idxnan[0] ne -1) then sinc[idxnan] = 1.d0
  ;Iint = (2.d0*p[0]*p[2]*eta0*sinc*sin(k0*(x-p[3])))+p[4]
  Iint = (2.d0*p[0]*p[2]*eta0*sinc*sin(k0*(x-p[3])+p[5]))+p[4]

;   gauss = p[0]*exp(-0.5d0*((x-p[3])/(lcoh*2.d0))^2.)
;   Iint = (2.d0*p[0]*p[2]*eta0*gauss*sin(k0*(x-p[3])+p[5]))+p[4]

  return, Iint

end

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

;=================================================================================================

function normalize, x

  norm = x
  idx = where(finite(x) eq 1)
  norm[idx] = (x[idx]-min(x[idx]))/(max(x[idx])-min(x[idx]))

  return, norm

end

;=================================================================================================

function snrfun, x, p
  fit = p[7]*( exp(-0.5d0*((x-p[0])/p[1])^2.) + p[4]*exp(-0.5d0*((x-p[5])/p[6])^2.)*sin(p[2]*x-p[3]) )+p[8]
  return, fit
end

;=================================================================================================

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

;=================================================================================================

function loadFSU, path, file, fsuid, opdcid

;read in FSU data
  if (fsuid eq 'A') then obs = mrdfits(path+file,7, /silent)
  ;if (fsuid eq 'B') then obs = mrdfits(path+file,8, /silent)

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

  ;remove non interesting states, where OPDC was not in 20 or 21
  idx = where(state eq 20 or state eq 21)
  time = time[idx]
  snr = snr[idx]
  gd = gd[idx]
  A = A[*,idx]
  B = B[*,idx]
  C = C[*,idx]
  D = D[*,idx]
  darkwhite = darkwhite[*,idx]
  fuoffset = fuoffset[idx]
  rtoffset = rtoffset[idx]
  state = state[idx]
  dlpos = dlpos[idx]

;remove 1st and last scan as they might not cross the whole fringe package
  ;remove 1st scan
  diff = ts_diff(state,1)
  idx = where(diff ne 0.d0)
  time = time[idx[0]+1:*]
  snr = snr[idx[0]+1:*]
  gd = gd[idx[0]+1:*]
  A = A[*,idx[0]+1:*]
  B = B[*,idx[0]+1:*]
  C = C[*,idx[0]+1:*]
  D = D[*,idx[0]+1:*]
  darkwhite = darkwhite[*,idx[0]+1:*]
  fuoffset = fuoffset[idx[0]+1:*]
  rtoffset = rtoffset[idx[0]+1:*]
  state = state[idx[0]+1:*]
  dlpos = dlpos[idx[0]+1:*]

  ;remove 2nd scan
  diff = ts_diff(state,1)
  idx = where(diff ne 0.d0)
  idx = idx[n_elements(idx)-1]	;ladt index
  time = time[0:idx]
  snr = snr[0:idx]
  gd = gd[0:idx]
  A = A[*,0:idx]
  B = B[*,0:idx]
  C = C[*,0:idx]
  D = D[*,0:idx]
  darkwhite = darkwhite[*,0:idx]
  fuoffset = fuoffset[0:idx]
  rtoffset = rtoffset[0:idx]
  state = state[0:idx]
  dlpos = dlpos[0:idx]

  data = {time:time, snr:snr, gd:gd, A:A, B:B, C:C, D:D, darkwhite:darkwhite, fuoffset:fuoffset, rtoffset:rtoffset, state:state, dlpos:dlpos}

  return, data

end

;=================================================================================================

function getFSUA_Dark, pathdf, darkfile

  dark = dblarr(4,6)

  st_d1 = mrdfits(pathdf+darkfile,6, /silent)

  idxd1 = where(st_d1.time ne 0.d0)


  for i=0,3 do begin	;loop over quadrants

    for j=0,5 do begin	;loop over pixel

      if (strmatch(darkfile, 'PACMAN*') eq 1) then begin

	if (i eq 0) then dark[i,j] = median(st_d1[idxd1].data1[j])
	if (i eq 1) then dark[i,j] = median(st_d1[idxd1].data2[j])
	if (i eq 2) then dark[i,j] = median(st_d1[idxd1].data3[j])
	if (i eq 3) then dark[i,j] = median(st_d1[idxd1].data4[j])

      endif else begin

	;use median instead of mean as there can be straylight if brighter objects are observed
	;situation should be improve if telescope chopping is implemented instead of using the ACUs
	if (i eq 0) then dark[i,j] = median(st_d1[idxd1].data1[j]-st_d1[idxd1].data5[0])
	if (i eq 1) then dark[i,j] = median(st_d1[idxd1].data2[j]-st_d1[idxd1].data5[1])
	if (i eq 2) then dark[i,j] = median(st_d1[idxd1].data3[j]-st_d1[idxd1].data5[2])
	if (i eq 3) then dark[i,j] = median(st_d1[idxd1].data4[j]-st_d1[idxd1].data5[3])

      endelse

    endfor

  endfor

  return, dark

end


;=================================================================================================

function getFSUA_DarkFlat, rawdata, pathdf, darkfile, flatfile1, flatfile2

  dark = dblarr(6,4)	;4 quadrants containing 6 pixels where 1st pixel is white light pixel
  flat = dblarr(6,4)	;final flat
  flat1 = dblarr(6,4)	;flat beam1
  flat2 = dblarr(6,4)	;flat beam2

  st_d1 = mrdfits(pathdf+darkfile,6, /silent)
  st_f1 = mrdfits(pathdf+flatfile1,6, /silent)
  st_f2 = mrdfits(pathdf+flatfile2,6, /silent)

  idxd1 = where(st_d1.time ne 0.d0)
  idxf1 = where(st_f1.time ne 0.d0)
  idxf2 = where(st_f2.time ne 0.d0)

; ;for normalized flat field
;   idxff = closest(st_f2[idxf2].time, st_f2[idxf1].time)
;   ff1a = st_f1[idxf1].data1[0]
;   ff1b = st_f1[idxf1].data2[0]
;   ff1c = st_f1[idxf1].data3[0]
;   ff1d = st_f1[idxf1].data4[0]
;   ff2a = st_f2[idxf2[idxff]].data1[0]
;   ff2b = st_f2[idxf2[idxff]].data2[0]
;   ff2c = st_f2[idxf2[idxff]].data3[0]
;   ff2d = st_f2[idxf2[idxff]].data4[0]

;   S0 = dblarr(4, n_elements(idxff))


;compute dark and flat and normalized flat exposure

  for i=0,3 do begin	;loop over quadrants

    for j=0,5 do begin	;loop over channels

      if (i eq 0) then begin
	dark[j,i] = median(st_d1[idxd1].data1[j])
	flat1[j,i] = median(st_f1[idxf1].data1[j])
	flat2[j,i] = median(st_f2[idxf2].data1[j])
      endif

      if (i eq 1) then begin 
	dark[j,i] = median(st_d1[idxd1].data2[j])
	flat1[j,i] = median(st_f1[idxf1].data2[j])
	flat2[j,i] = median(st_f2[idxf2].data2[j])
      endif

      if (i eq 2) then begin 
	dark[j,i] = median(st_d1[idxd1].data3[j])
	flat1[j,i] = median(st_f1[idxf1].data3[j])
	flat2[j,i] = median(st_f2[idxf2].data3[j])
      endif

      if (i eq 3) then begin
	dark[j,i] = median(st_d1[idxd1].data4[j])
	flat1[j,i] = median(st_f1[idxf1].data4[j])
	flat2[j,i] = median(st_f2[idxf2].data4[j])
      endif

    endfor

  endfor

  flat = flat1 + flat2	;combine flats of beam 1 and 2 to final flat

;   ;normalized FF
;   S0 = dblarr(4)
;   for i=0,3 do  S0[i] = ((flat1[0,i]-dark[0,i])+(flat2[0,i]-dark[0,i]))/(flat[0,i]-2.d0*dark[0,i])

  rawA = rawdata.A
  rawB = rawdata.B
  rawC = rawdata.C
  rawD = rawdata.D

  rawA[0,*] = rawA[0,*]-rawdata.darkwhite[0,*]
  rawB[0,*] = rawB[0,*]-rawdata.darkwhite[1,*]
  rawC[0,*] = rawC[0,*]-rawdata.darkwhite[2,*]
  rawD[0,*] = rawD[0,*]-rawdata.darkwhite[3,*]

  A = dblarr(6,n_elements(rawdata.A[0,*])) & B = A & C = A & D = A
  PD = A
  for i=0,5 do begin

      A[i,*] = (rawA[i,*] - dark[i,0]) / (flat[i,0] - 2.d0*dark[i,0])
      B[i,*] = (rawB[i,*] - dark[i,1]) / (flat[i,1] - 2.d0*dark[i,1])
      C[i,*] = (rawC[i,*] - dark[i,2]) / (flat[i,2] - 2.d0*dark[i,2])
      D[i,*] = (rawD[i,*] - dark[i,3]) / (flat[i,3] - 2.d0*dark[i,3])

  endfor

  data = {A:A, B:B, C:C, D:D}

  return, data

end

;=================================================================================================

function func_phot_offset, x, p

  fit = x*p[0]
  return, fit

end

;=================================================================================================


pro calibrate_FSUA_SKY_ffscan_filter_old, dodark=dodark, nodark=nodark

;stupid keyword check, enforce the usage of one option
check = 0
if keyword_set(dodark) then check = 1
if keyword_set(nodark) then check = 2

if (check eq 0) then begin

  print, 'Use /dodark or /nodark option! NOdark preferred.'
  retall

endif

;define path and fsu 2nd fringe scan file
;   path = ''
;   print, ''
;   read, 'Enter path, e.g. /media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_P90RATZKA/FSUAdata/FringeScan: ', path
;   if (path eq '') then begin
;     print, 'You have to provide a path.'
;     return
;   endif
;   path = path+'/'
  path = '/home/amueller/work/MIDI_FSUA/Pipeline/Calibration/'

  print,''
  print,'Following FSUA FRINGE SCAN files found: '
  tmp = file_search(path+'PACMAN_OBJ_SCAN_*.fits', count=nf)
  if (tmp[0] eq '') then begin

    print, ''
    print, 'No FSUA SECOND FRINGE SCAN files found'
    stop

  endif else begin

    for i=0,nf-1 do begin
      print, uint(i+1), ' ', tmp[i]
    endfor

    print, ''

    read,'Enter which file you want to use: ', usefiles
;    usefiles = 1
;     usefiles = 1
    print, ''

    file = tmp[usefiles-1]

  endelse

;dirctory and dark and flat files
pathdf = path+'SkyCalibration/'

; darkfile = 'FSUA_SKY_DARK_132_0013.fits'
; darkfile = 'FSUA_SKY_DARK_333_0001.fits'
; flatfile1 = 'FSUA_SKY_FLAT_333_0001.fits'
; flatfile2 = 'FSUA_SKY_FLAT_333_0002.fits'
; darkfile = 'FSUA_SKY_DARK_350_0008.fits'
; flatfile1 = 'FSUA_SKY_FLAT_350_0015.fits'
; flatfile2 = 'FSUA_SKY_FLAT_350_0016.fits'

;=================================================================================================

;remove full path in front of file
  pos1 = strpos(file, 'PACMAN_OBJ_SCAN')
  pos2 = strpos(file, '.fits')
  file = strmid(file, pos1, pos2-pos1+5)


;extract file ID
  pos1 = strpos(file, 'SCAN_')
  pos2 = strpos(file, '.fits')
  fileid = strmid(file, pos1+5, pos2-pos1-5)


;define path for results
  resultpath = path+'SkyCal_'+fileid+'/'
  spawn, 'mkdir '+resultpath

;extract keywords present in OBJ_SCAN files
dum = readfits(path+file, hdr, /silent)
  base = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN BASE'))	;center position
  channel = (get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN CHANNEL'))
  nslew = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN NSLEW'))	;number of slews
  samprate = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN SAMPRATE'));number of points/micron
  stroke = double(get_eso_keyword(hdr, 'HIERARCH ESO INS OPDSCAN STROKE'))	;meter
  mjd = double(get_eso_keyword(hdr, 'MJD-OBS'))
  airmass = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS AIRM START'))
  altitude = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS ALT'))
  azimut = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS AZ'))
  seeing = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS AMBI FWHM START'))
  tau0 = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS AMBI TAU0 START'))*1.d3	;ms
  opl = abs(double(get_eso_keyword(hdr, 'HIERARCH ESO DEL DLT1 OPL START'))-double(get_eso_keyword(hdr, 'HIERARCH ESO DEL DLT2 OPL START')))
  bl = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PBL12 START'))
  object = get_eso_keyword(hdr, 'HIERARCH ESO OBS TARG NAME')
  station = get_eso_keyword(hdr, 'HIERARCH ESO ISS CONF STATION1')+'-'+get_eso_keyword(hdr, 'HIERARCH ESO ISS CONF STATION2')

print, '     sampling       stroke    nslew         seeing             tau'
print, round(samprate), round(stroke*1.d6), round(nslew), seeing, tau0

;=================================================================================================

;LoadFSUA
fsuid = 'A'
opdcid = 'opdc'
rawdata = loadFSU(path, file, fsuid, opdcid)

;=================================================================================================

;LoadDarkFlat
; datadf = getFSUA_DarkFlat(rawdata, pathdf, darkfile, flatfile1, flatfile2)
;use photometry from scans instead

if keyword_set(dodark) then dark = getFSUA_Dark(path, darkfile)

;=================================================================================================

;SelectScans
diff = ts_diff(rawdata.state,1)
idx1 = where(diff eq -1.d0)
idx2 = where(diff eq 1.d0)
idxscan = [idx1,idx2]
idxscan = idxscan[sort(idxscan)]

; idxscan = idxscan[0:35]
; print, ''
; print, 'REMOVE LINE 433'
; print, ''

nscans = n_elements(idxscan)-1

if (total(idxscan) eq -2.) then begin
  print, 'Error. No Fringe Scans Found.'
  return
endif


max_snr = dblarr(nscans)
snr_pos = dblarr(nscans)
nval = dblarr(nscans)	;number of values per scan

for i=0,nscans-1 do begin

  snr = rawdata.snr[idxscan[i]+1:idxscan[i+1]]
  rtoffset = rawdata.rtoffset[idxscan[i]+1:idxscan[i+1]]
  nval[i] = n_elements(snr)

  dum = max(snr, idxsnr)
  startval = [5.d0, rtoffset[idxsnr], 1.d-5, 1.d0]

;   A = rawdata.A[*,idxscan[i]+1:idxscan[i+1]]
;   B = rawdata.B[*,idxscan[i]+1:idxscan[i+1]]
;   C = rawdata.C[*,idxscan[i]+1:idxscan[i+1]]
;   D = rawdata.D[*,idxscan[i]+1:idxscan[i+1]]
;   window, 0, xs=1500, ys=1000
;   !p.multi=[0,4,6]
;   for i=0,5 do begin
;     plot, A[i,*]
;     plot, B[i,*]
;     plot, C[i,*]
;     plot, D[i,*]
;   endfor
;   !p.multi=[0,1,0]
;   window, 1, xs=200, ys=200
;   plot, snr, yr=[0,10]
;   hak

  gauss_nterms = '4'
  dummy_err = dblarr(n_elements(rtoffset))
  dummy_err[*] = 1.d0
  fitparams_gauss = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', rtoffset, smooth(snr, 3), dummy_err, startval, weights=snr, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
  ;sigma = perror*sqrt(bestnorm/dof)
  max_snr[i] = fitparams_gauss[0]+fitparams_gauss[3]
  snr_pos[i] = fitparams_gauss[1]

;   window, 0
;   plot, rtoffset, smooth(snr, 3), xst=1
;   oplot, rtoffset, yfit, color=fsc_color('red')
;   hak

  proceeding_text,loop=(nscans), i=i, prompt='> Fit SNR                        '+string(i+1,form='(I4)')

endfor

time = dblarr(nscans, max(nval))
snr = dblarr(nscans, max(nval))
gd = dblarr(nscans, max(nval))
rawA = dblarr(nscans, 6, max(nval))
rawB = dblarr(nscans, 6, max(nval))
rawC = dblarr(nscans, 6, max(nval))
rawD = dblarr(nscans, 6, max(nval))
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

;=================================================================================================

;shorten arrays because there can be zeros in the arrays because of unequal number of values of the individual scans

  idx = where(time eq 0.d0)
  if (idx[0] ne -1) then begin

    for i=0,nscans-1 do begin

      idx = where(time[i,*] eq 0.d0)
      if (idx[0] ne -1) then begin

	time = time[*,0:idx-1]
	snr = snr[*,0:idx-1]
	gd = gd[*,0:idx-1]
	fuoffset = fuoffset[*,0:idx-1]
	rtoffset = rtoffset[*,0:idx-1]
	state = state[*,0:idx-1]
	rawA = rawA[*,*,0:idx-1]
	rawB = rawB[*,*,0:idx-1]
	rawC = rawC[*,*,0:idx-1]
	rawD = rawD[*,*,0:idx-1]

      endif

    endfor

  endif

  ndata = n_elements(rawA[0,0,*])

; 
; ;check for 2 consecutive scans
; idx1 = where(time[0,*] eq 0.d0)
; idx2 = where(time[1,*] eq 0.d0)
; 
; if (idx1[0] ne -1) then begin
; 
;   time = time[0:idx1-1]
;   snr = snr[0:idx1-1]
;   gd = gd[0:idx1-1]
;   fuoffset = fuoffset[0:idx1-1]
;   rtoffset = rtoffset[0:idx1-1]
;   state = state[*,0:idx1-1]
;   rawA = rawA[*,*,0:idx1-1]
;   rawB = rawB[*,*,0:idx1-1]
;   rawC = rawC[*,*,0:idx1-1]
;   rawD = rawD[*,*,0:idx1-1]
; 
; endif
; 
; if (idx2[0] ne -1) then begin
; 
;   time = time[*,0:idx2-1]
;   snr = snr[*,0:idx2-1]
;   gd = gd[*,0:idx2-1]
;   fuoffset = fuoffset[*,0:idx2-1]
;   rtoffset = rtoffset[*,0:idx2-1]
;   state = state[*,0:idx2-1]
;   rawA = rawA[*,*,0:idx2-1]
;   rawB = rawB[*,*,0:idx2-1]
;   rawC = rawC[*,*,0:idx2-1]
;   rawD = rawD[*,*,0:idx2-1]
; 
; endif

;=================================================================================================

;remove scans based on deviation from maxSNR position

resistant_mean, snr_pos, 2.0, robmeansnrpos, sigma, NumRej, GoodVec=GoodInd, /double

time = time[GoodInd,*]
snr = snr[GoodInd, *]
snr_pos = snr_pos[GoodInd]
max_snr = max_snr[GoodInd]
gd = gd[GoodInd, *]
fuoffset = fuoffset[GoodInd, *]
rtoffset = rtoffset[GoodInd, *]
state = state[GoodInd, *]
rawA = rawA[GoodInd, *, *]
rawB = rawB[GoodInd, *, *]
rawC = rawC[GoodInd, *, *]
rawD = rawD[GoodInd, *, *]
nscans = n_elements(GoodInd)

;=================================================================================================

;get photometry from each scan and 'normalize' each scan by dividing the real time photometry
;subtract dark and divide by flat
;cut out roughly the fringes in order to be able to correct for a residual offset between A-C and B-D photometry
cut = [2.0d-05, 3.5d-05, 4.0d-05, 5.d-05, 5.5d-05, 7.d-05]

A = dblarr(nscans, 6, n_elements(rawA[0,0,*])) & B=A & C=A & D=A
rawphotAC = dblarr(nscans, 6, n_elements(rawA[0,0,*]))
rawphotBD = dblarr(nscans, 6, n_elements(rawA[0,0,*]))
rawphotA = dblarr(nscans, 6, n_elements(rawA[0,0,*])) & rawphotB = rawphotA & rawphotC = rawphotA & rawphotD = rawphotA
photA = dblarr(nscans, 6, n_elements(rawA[0,0,*])) & photB = photA & photC = photA & photD = photA
offAC = dblarr(nscans, 6)
offBD = dblarr(nscans, 6)

maxiter = 1000.

for i=0,nscans-1 do begin

  for j=0,5 do begin

    idx = where(rtoffset[i,*] le snr_pos[i]-cut[j] or rtoffset[i,*] ge snr_pos[i]+cut[j])
    dumerr = n_elements(idx)
    dumerr[*] = 1.d0
    start_val = [1.d0]
    off1 = MPFITfun('func_phot_offset', rawC[i,j,idx], rawA[i,j,idx], dumerr, start_val, weights = 1.d0, $
	    maxiter=maxiter, niter=niter, status=status, bestnorm=bestnorm, $
	    perror=perror, covar=covar, yfit=yfit, dof=dof, /quiet)
    offAC[i,j] = off1[0]

    off2 = MPFITfun('func_phot_offset', rawD[i,j,idx], rawB[i,j,idx], dumerr, start_val, weights = 1.d0, $
	    maxiter=maxiter, niter=niter, status=status, bestnorm=bestnorm, $
	    perror=perror, covar=covar, yfit=yfit, dof=dof, /quiet)
    offBD[i,j] = off2[0]

    rawphotAC[i,j,*] = rawA[i,j,*]+rawC[i,j,*]*offAC[i,j]
    rawphotBD[i,j,*] = rawB[i,j,*]+rawD[i,j,*]*offBD[i,j]

    ;raw photometry = flatfield
    rawphotA[i,j,*] = rawphotAC[i,j,*]/2.d0
    rawphotC[i,j,*] = rawphotAC[i,j,*]/2.d0/offAC[i,j]
    rawphotB[i,j,*] = rawphotBD[i,j,*]/2.d0
    rawphotD[i,j,*] = rawphotBD[i,j,*]/2.d0/offBD[i,j]

;     ;plot fringe scan and real time photometry
;     window, 0, xs=800, ys=1000
;     !p.multi=[0,1,4]
;     plot, rawA[i,j,*], xst=1
;     oplot, rawphotA[i,j,*], color=rgb(255,0,0)
;     plot, rawB[i,j,*], xst=1
;     oplot, rawphotB[i,j,*], color=rgb(255,0,0)
;     plot, rawC[i,j,*], xst=1
;     oplot, rawphotC[i,j,*], color=rgb(255,0,0)
;     plot, rawD[i,j,*], xst=1
;     oplot, rawphotD[i,j,*], color=rgb(255,0,0)
;     hak
;     !p.multi=[0,1,0]

    if keyword_set(dodark) then begin

      photA[i,j,*] = rawphotA[i,j,*] - 2.d0*dark[0,j]
      photB[i,j,*] = rawphotB[i,j,*] - 2.d0*dark[1,j]
      photC[i,j,*] = rawphotC[i,j,*] - 2.d0*dark[2,j]
      photD[i,j,*] = rawphotD[i,j,*] - 2.d0*dark[3,j]

      A[i,j,*] = (rawA[i,j,*]-dark[0,j])/photA[i,j,*]
      B[i,j,*] = (rawB[i,j,*]-dark[1,j])/photB[i,j,*]
      C[i,j,*] = (rawC[i,j,*]-dark[2,j])/photC[i,j,*]
      D[i,j,*] = (rawD[i,j,*]-dark[3,j])/photD[i,j,*]

    endif

    if keyword_set(nodark) then begin

      photA[i,j,*] = rawphotA[i,j,*]
      photB[i,j,*] = rawphotB[i,j,*]
      photC[i,j,*] = rawphotC[i,j,*]
      photD[i,j,*] = rawphotD[i,j,*]

      A[i,j,*] = rawA[i,j,*]/photA[i,j,*]
      B[i,j,*] = rawB[i,j,*]/photB[i,j,*]
      C[i,j,*] = rawC[i,j,*]/photC[i,j,*]
      D[i,j,*] = rawD[i,j,*]/photD[i,j,*]

    endif

  endfor

  proceeding_text,loop=(nscans), i=i, prompt='> Photometry                     '+string(i+1,form='(I4)')

endfor


;use OPD values of first scan as reference
grid = rtoffset
; refgrid = reform(grid[0,*])

;=================================================================================================

;filter scans
step = abs(median(ts_diff(transpose(grid[0,*]),1)))
ndata = double(n_elements(A[0,0,*]))

filtA = dblarr(nscans, 6, ndata)
filtB = filtA & filtC = filtA & filtD = filtA

dfreq = 1.d0/(step*ndata)
auxgrid = dindgen(ndata)
wavenum = dfreq*auxgrid

;values are measured manually from a high SNR 
;		white			1		2			3				4			5
winwn = [[3.8d5, 5.3d5], [466937.42d0, 530228.56d0], [444433.d0, 500692.d0], [420000.28d0, 475376.d0], [399426.d0, 451466.d0], [382549.88d0, 433182.d0]]

for i=0,nscans-1 do begin

  for j=0,5 do begin

    idx1 = where(wavenum gt winwn[0,j] and wavenum lt winwn[1,j])	;white light
    idx2 = n_elements(wavenum)-idx1
    idx = [[idx1],[idx2]]
    filter = dblarr(ndata)
    filter[idx] = 1.d0

    tmpA = fft(A[i,j,*],-1)
    tmpB = fft(B[i,j,*],-1)
    tmpC = fft(C[i,j,*],-1)
    tmpD = fft(D[i,j,*],-1)

    tmpA2 = tmpA*filter
    tmpB2 = tmpB*filter
    tmpC2 = tmpC*filter
    tmpD2 = tmpD*filter

    filtA[i,j,*] = normalize(fft(tmpA2,1))
    filtB[i,j,*] = normalize(fft(tmpB2,1))
    filtC[i,j,*] = normalize(fft(tmpC2,1))
    filtD[i,j,*] = normalize(fft(tmpD2,1))

  endfor

  proceeding_text,loop=(nscans), i=i, prompt='> Filter scans                   '+string(i+1,form='(I4)')

endfor

;=================================================================================================

;recompute SNR on filtered scans

;Phase shift errors are median values from previous calibrations
;pse = array: channel, pixel / 4,6
pse = [[-0.33825935d0, 1.9410529d0, 1.3795450d0, 1.4108432d0], [-0.33967445d0, 1.9405431d0, 1.5009768d0, 1.3153251d0], [-0.35784562d0, 1.9337808d0, 1.5444702d0, 1.2025047d0], [-0.38910486d0, 1.9211935d0, 1.4401357d0, 1.3150547d0], [-0.41229578d0, 1.9110500d0, 1.3632777d0, 1.3905010d0], [-0.39482077d0, 1.9187582d0, 1.3212059d0, 1.4338347d0]]


csc = dblarr(6)
for i=0,5 do csc[i] = 1.d0/(pse[1,i]*pse[2,i]-pse[0,i]*pse[3,i])

X = dblarr(nscans, 6, ndata) & Y = X & filtSNR = X
q=x & w=x & fs = x

for i=0,nscans-1 do begin

  for j=0,5 do begin

    X[i,j,*] = ((pse[2,j]*(filtA[i,j,*]-filtC[i,j,*])) - (pse[0,j]*(filtB[i,j,*]-filtD[i,j,*])));*csc[j]
    Y[i,j,*] = ((pse[1,j]*(filtB[i,j,*]-filtD[i,j,*])) - (pse[3,j]*(filtA[i,j,*]-filtC[i,j,*])));*csc[j]
    filtSNR[i,j,*] = sqrt(X[i,j,*]^2. + Y[i,j,*]^2.)

  endfor

  proceeding_text,loop=(nscans), i=i, prompt='> Recompute SNR of filtered scans'+string(i+1,form='(I4)')

endfor

;=================================================================================================
;fit new SNR computed from filtered scans

filt_max_snr = dblarr(nscans)
filt_snr_pos = dblarr(nscans)

for i=0,nscans-1 do begin

  dum = max(filtSNR[i,0,*], idxsnr)
  startval = [0.5, rtoffset[i,idxsnr], 1.d-5, 0.d0]

  gauss_nterms = '4'
  dummy_err = dblarr(ndata)
  dummy_err[*] = 1.d0
  fitparams_gauss = mpfitfun('fit_gauss_'+strcompress(gauss_nterms, /rem)+'terms', rtoffset[i,*], filtSNR[i,0,*], dummy_err, startval, weights=filtSNR[i,0,*], maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
  ;sigma = perror*sqrt(bestnorm/dof)
  filt_max_snr[i] = fitparams_gauss[0]+fitparams_gauss[3]
  filt_snr_pos[i] = fitparams_gauss[1]

;   window, 0
;   plot, rtoffset[i,*], filtSNR[i,0,*], xst=1
;   oplot, rtoffset[i,*], yfit, color=fsc_color('red')
;   hak

  proceeding_text,loop=(nscans), i=i, prompt='> Fit SNR of filtered scans      '+string(i+1,form='(I4)')

endfor


;=================================================================================================

;select only high SNR scans

resistant_mean, filt_max_snr, 1.0, robmeansnr, sigmaSNR, NumRejSNR, GoodVec=GoodIndSNR, /double
;we don't want reject high SNR data
BadIndSNR = intarr(nscans)
for i=0,nscans-1 do begin

  idx = where(GoodIndSNR eq i)
  if (idx[0] eq -1) then BadIndSNR[i] = i+1	;+1 because 0 can be missing as well

endfor

idx = where(BadIndSNR ne 0)
BadIndSNR = BadIndSNR[idx]-1

idx = where(filt_max_snr[BadIndSNR] lt median(filt_max_snr))
if (idx[0] ne -1) then BadIndSNR = BadIndSNR[idx]

GoodIndSNR = indgen(nscans)
for i=0,nscans-1 do begin

  idx = where(GoodIndSNR[i] eq BadIndSNR)
  if (idx[0] ne -1) then GoodIndSNR[i] = -99

endfor

idx = where(GoodIndSNR ne -99)
if (idx[0] ne -1) then GoodIndSNR = GoodIndSNR[idx]

filt_snr_pos = filt_snr_pos[GoodIndSNR]
filt_max_snr = filt_max_snr[GoodIndSNR]
filtA = filtA[GoodIndSNR, *, *]
filtB = filtB[GoodIndSNR, *, *]
filtC = filtC[GoodIndSNR, *, *]
filtD = filtD[GoodIndSNR, *, *]
grid = grid[GoodIndSNR, *]
photA = photA[GoodIndSNR,*,*]
photB = photB[GoodIndSNR,*,*]
photC = photC[GoodIndSNR,*,*]
photD = photD[GoodIndSNR,*,*]
gd = gd[GoodIndSNR, *]
nscans = n_elements(GoodIndSNR)


;remove outliers based on maxSNR position (problem with e.g. FSUA_SECOND_FRINGE_SCAN_332_0003.fits)
resistant_mean, filt_snr_pos, 2.5, robmeansnr, sigmaSNR, NumRejSNR, GoodVec=GoodIndSNR, /double

filt_snr_pos = filt_snr_pos[GoodIndSNR]
filt_max_snr = filt_max_snr[GoodIndSNR]
filtSNR = filtSNR[GoodIndSNR, *,*]
filtA = filtA[GoodIndSNR, *, *]
filtB = filtB[GoodIndSNR, *, *]
filtC = filtC[GoodIndSNR, *, *]
filtD = filtD[GoodIndSNR, *, *]
grid = grid[GoodIndSNR, *]
photA = photA[GoodIndSNR,*,*]
photB = photB[GoodIndSNR,*,*]
photC = photC[GoodIndSNR,*,*]
photD = photD[GoodIndSNR,*,*]
gd = gd[GoodIndSNR, *]
state = state[GoodIndSNR,*]
nscans = n_elements(GoodIndSNR)


filtABCD = dblarr(nscans, 4, 6, ndata)
filtABCD[*,0,*, *] = filtA
filtABCD[*,1,*, *] = filtB
filtABCD[*,2,*, *] = filtC
filtABCD[*,3,*, *] = filtD

photABCD = dblarr(nscans, 4, 6, ndata)
photABCD[*,0,*,*] = photA
photABCD[*,1,*,*] = photB
photABCD[*,2,*,*] = photC
photABCD[*,3,*,*] = photD


;=================================================================================================

;fit fringe packages

fitparams = dblarr(nscans,4,6,6)	;channel-pixel-params
fringefit = dblarr(nscans,4,6,ndata)

  ;p0: Intensity, p1: wavelength, p2: bandwidth, p3: ZOPD offset, p4: Ioffset, p5: Phase
startval = [[500000.d0, 2.25d-6, 5.d-7, base, 0.5d0, 0.d0],$	;base only dummy, use filt_snr_pos
	    [1000000.d0, 2.0d-6, 1.4d-7, base, 0.5d0, 0.d0],$
	    [1000000.d0, 2.1d-6, 1.4d-7, base, 0.5d0, 0.d0],$
	    [1000000.d0, 2.2d-6, 1.4d-7, base, 0.5d0, 0.d0],$
	    [1000000.d0, 2.3d-6, 1.4d-7, base, 0.5d0, 0.d0],$
	    [1000000.d0, 2.4d-6, 1.4d-7, base, 0.5d0, 0.d0]]

;phasestart = ([0.d0, 90.d0, 180.d0, 270.d0])*!dtor
;phasestart = ([0.d0, 50.d0, 150.d0, 215.d0])*!dtor
;phasestart = ([0.d0, 0.d0, 0.d0, 0.d0])*!dtor

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
; pi[0].limits[0] = 0.d0	;Intensity has to be greater than 0
; pi[0].limited[0] = 1
pi[1].limits[0] = 1.6d-6	;wavelength limited between 1.8 and 2.7 micrometer
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
; pi[5].limited[0] = 1
; pi[5].limited[1] = 1


dummy_err = dblarr(ndata)
dummy_err[*] = 1.d0
weight = dblarr(4,6,ndata)

for i=0,nscans-1 do begin	;scans

  startval[3,*] = filt_snr_pos[i]

  for j=0,3 do begin	;FSU channels

    for k=0,5 do begin	;spectral pixel

      lcoh = (startval[1,k]^2.-(startval[2,k]^2.)/4.d0)/startval[2,k]
      weight[j,k,*] = (abs((sin(!DPI*(grid[i,*]-startval[3,k])/lcoh))/(!DPI*(grid[i,*]-startval[3,k])/lcoh)))^2.

      tmpparam = mpfitfun('fit_fringe', grid[i,*], filtABCD[i,j,k,*], dummy_err, startval[*,k], weights=weight[j,k,*], maxiter=100, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet, /nan, parinfo=pi)
      fitparams[i,j,k,*] = tmpparam
      fringefit[i,j,k,*] = yfit

    endfor

;     ;plot all spectra of one channel
;     window, 0, xs=1000, ys=900
;     !p.multi=[0,1,6]
;     for k=0,5 do begin
;       plot, grid[i,*], filtABCD[i,j,k,*], xst=1
;       oplot, grid[i,*], weight[j,k,*], color=rgb(0,255,0)
;       oplot, grid[i,*], fringefit[i,j,k,*], color=rgb(255,0,0)
;     endfor
;     !p.multi=[0,1,5]
;     hak

  endfor

  proceeding_text,loop=(nscans), i=i, prompt='> Fit fringe package             '+string(i+1,form='(I4)')

endfor

wlen = dblarr(4,6) & wlenerr = wlen

for i=0,3 do begin

  for j=0,5 do begin

    resistant_mean, fitparams[*,i,j,1], 1.0, robmean, sigma, NumRej, GoodVec=GoodInd, /double
    wlen[i,j] = median(fitparams[GoodInd,i,j,1])
    wlenerr[i,j] = stddev(fitparams[GoodInd,i,j,1])

  endfor

endfor

lambda = dblarr(6)
for i=0,5 do lambda[i] = median(wlen[*,i])	;median waelength of spectral pixel of all channels

;=================================================================================================

;create high-resolution fringes based on fits

hrgrid = dblarr(nscans, 2000)
ndata = n_elements(hrgrid[0,*])
hrABCD = dblarr(nscans, 4, 6, ndata)	;computed fringe package

for i=0,nscans-1 do begin

  hrgrid[i,*] = linspace(filt_snr_pos[i]-100.d-6, filt_snr_pos[i]+100.d-6, ndata)

  for j=0,3 do begin

    for k=0,5 do begin

      hrABCD[i,j,k,*] = normalize(fit_fringe(hrgrid[i,*], reform(fitparams[i,j,k,*])))

    endfor

  endfor

;   window, 0, xs=1500, ys=500	;white light pixel
;   plot, hrgrid[i,*], hrABCD[i,0,0,*], xst=1, xr=[median(filt_snr_pos)-12.d-6,median(filt_snr_pos)+12.d-6]
;   oplot, hrgrid[i,*], hrABCD[i,1,0,*], color=fsc_color('red')
;   oplot, hrgrid[i,*], hrABCD[i,2,0,*], color=fsc_color('yellow')
;   oplot, hrgrid[i,*], hrABCD[i,3,0,*], color=fsc_color('green')
;   legend, ['A','B','C','D'], color=[fsc_color('white'),fsc_color('red'),fsc_color('yellow'),fsc_color('green')], margin=0, box=0, /left, linestyle=0, charsize=1.5
;   hak

  proceeding_text,loop=(nscans), i=i, prompt='> Create high resolution fringes '+string(i+1,form='(I4)')

endfor

;=================================================================================================

;get phases via cross correlation of high-resolution fringe

stepsize = abs(median(ts_diff(reform(hrgrid[0,*]),1)))	;stepsize same for all scans
lag = [dindgen(ndata-1.d0)-(ndata-1.d0), dindgen(ndata)]
cc = dblarr(nscans,3,6,n_elements(lag))	;only '3' bacause CC between AB,AC,AD
ind_maxposcc = intarr(nscans,3,6,1)
dabcd = dblarr(nscans,3,6,1) & dabcd_rad = dabcd
angrad = dblarr(nscans,3,6,1) & angdegall = angrad
phserrall = dblarr(nscans,4,6,1)

for i=0,nscans-1 do begin

  for j=0,2 do begin

    for k=0,5 do begin

      cc[i,j,k,*] = c_correlate(hrABCD[i,0,k,*], hrABCD[i,j+1,k,*], lag, /double)
      dum = max(cc[i,j,k,*], tmp)
      ind_maxposcc[i,j,k,*] = tmp
      dabcd[i,j,k,*] = abs(lag[ind_maxposcc[i,j,k,*]])*stepsize
      dabcd_rad[i,j,k,*] = (dabcd[i,j,k,*]*2.d0*!DPI)/wlen[j,k]

      if (j eq 0) then angrad[i,j,k,*] = dabcd_rad[i,j,k,0]-(!DPI/2.d0)
      if (j eq 1) then angrad[i,j,k,*] = dabcd_rad[i,j,k,0]-!DPI
      if (j eq 2) then angrad[i,j,k,*] = dabcd_rad[i,j,k,0]-(3.d0*!DPI/2.d0)

    endfor

    if (j eq 0) then begin
      phserrall[i,j,*,*] = sin(angrad[i,0,*,*])
      angdegall[i,j,*,*] = 90.d0+angrad[i,j,*,*]/!dtor
    endif
    if (j eq 1) then begin
      phserrall[i,j,*,*] = 1.d0+cos(angrad[i,1,*,*])
      angdegall[i,j,*,*] = 180.d0+angrad[i,j,*,*]/!dtor
    endif
    if (j eq 2) then begin
      phserrall[i,j,*,*] = cos(angrad[i,0,*,*])+cos(angrad[i,2,*,*])
      angdegall[i,j,*,*] = 270.d0+angrad[i,j,*,*]/!dtor
    endif

  endfor

  phserrall[i,3,*,*] = -sin(angrad[i,0,*,*])-sin(angrad[i,2,*,*])


  proceeding_text,loop=(nscans), i=i, prompt='> Compute angular phase shifts   '+string(i+1,form='(I4)')

endfor

angdegall = angdegall mod 360.d0

for i=0,nscans-1 do begin

  for j=0,2 do begin

    for k=0,5 do begin

      if (j eq 0) then if (angdegall[i,j,k] gt 180.d0) then angdegall[i,j,k] = 360.d0-angdegall[i,j,k]
      if (j eq 1) then if (angdegall[i,j,k] gt 270.d0) then angdegall[i,j,k] = 360.d0-angdegall[i,j,k]
      if (j eq 2) then if (angdegall[i,j,k] lt 180.d0) then angdegall[i,j,k] = 360.d0-angdegall[i,j,k]

    endfor

  endfor

endfor

;=================================================================================================

;average phase shifts and phserr

angdeg = dblarr(3,6) & angdegerr = angdeg 
phserr = dblarr(4,6) & phserrerr = phserr

for i=0,2 do begin
  for j=0,5 do begin

    resistant_mean, angdegall[*,i,j], 1.0, robmean, sigma, NumRej, GoodVec=GoodInd, /double
    angdeg[i,j] = median(angdegall[GoodInd,i,j])
    angdegerr[i,j] = stddev(angdegall[GoodInd,i,j])

  endfor
endfor

for i=0,3 do begin
  for j=0,5 do begin

    resistant_mean, phserrall[*,i,j], 1.0, robmean, sigma, NumRej, GoodVec=GoodInd, /double
    phserr[i,j] = median(phserrall[GoodInd,i,j])
    phserrerr[i,j] = stddev(phserrall[GoodInd,i,j])

  endfor
endfor


; ;get phases via cross correlation of observed fringes
; 
; stepsize = abs(median(ts_diff(reform(grid[0,*]),1)))	;stepsize same for all scans
; lag = [dindgen(ndata-1.d0)-(ndata-1.d0), dindgen(ndata)]
; cc = dblarr(nscans,3,6,n_elements(lag))	;only '3' bacause CC between AB,AC,AD
; ind_maxposcc = intarr(nscans,3,6,1)
; dabcd = dblarr(nscans,3,6,1) & dabcd_rad = dabcd
; angrad = dblarr(nscans,3,6,1) & angdeg = angrad
; phserr = dblarr(nscans,3,6,1)
; 
; for i=0,nscans-1 do begin
; 
;   for j=0,2 do begin
; 
;     for k=0,5 do begin
; 
;       cc[i,j,k,*] = c_correlate(filtABCD[i,0,k,*], filtABCD[i,j+1,k,*], lag, /double)
;       dum = max(cc[i,j,k,*], tmp)
;       ind_maxposcc[i,j,k,*] = tmp
;       dabcd[i,j,k,*] = abs(lag[ind_maxposcc[i,j,k,*]])*stepsize
;       dabcd_rad[i,j,k,*] = (dabcd[i,j,k,*]*2.d0*!DPI)/wlen[j,k]
; 
;       if (j eq 0) then angrad[i,j,k,*] = dabcd_rad[i,j,k,0]-(!DPI/2.d0)
;       if (j eq 1) then angrad[i,j,k,*] = dabcd_rad[i,j,k,0]-!DPI
;       if (j eq 2) then angrad[i,j,k,*] = dabcd_rad[i,j,k,0]-(3.d0*!DPI/2.d0)
; 
;     endfor
; 
;     if (j eq 0) then begin
;       phserr[i,j,*,*] = sin(angrad[i,0,*,*])
;       angdeg[i,j,*,*] = 90.d0+angrad[i,j,*,*]/!dtor
;     endif
;     if (j eq 1) then begin
;       phserr[i,j,*,*] = 1.d0+cos(angrad[i,1,*,*])
;       angdeg[i,j,*,*] = 180.d0+angrad[i,j,*,*]/!dtor
;     endif
;     if (j eq 2) then begin
;       phserr[i,j,*,*] = cos(angrad[i,0,*,*])+cos(angrad[i,2,*,*])
;       angdeg[i,j,*,*] = 270.d0+angrad[i,j,*,*]/!dtor
;     endif
;     if (j eq 3) then phserr[i,j,*,*] = -sin(angrad[i,0,*,*])-sin(angrad[i,2,*,*])
; 
;   endfor
; 
;   proceeding_text,loop=(nscans), i=i, prompt='> Compute angular phase shifts   '+string(i+1,form='(I4)')
; 
; endfor
; 
; angdeg = angdeg mod 360.d0

;=================================================================================================


;Wavelength Quality Control
qcwlen = strarr(4,6)
qcfitwlen = strarr(4,6)
for i=0,3 do begin
  for j=0,5 do begin

    if (j eq 0) then begin
      if (wlen[i,j] gt 2.2d-6 and wlen[i,j] lt 2.3d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
;       if (fitwlen[i,j] gt 2.2d-6 and fitwlen[i,j] lt 2.3d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 1) then begin
      if (wlen[i,j] gt 1.95d-6 and wlen[i,j] lt 2.1d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
;       if (fitwlen[i,j] gt 1.95d-6 and fitwlen[i,j] lt 2.1d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 2) then begin
      if (wlen[i,j] gt 2.05d-6 and wlen[i,j] lt 2.2d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
;       if (fitwlen[i,j] gt 2.05d-6 and fitwlen[i,j] lt 2.2d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 3) then begin
      if (wlen[i,j] gt 2.15d-6 and wlen[i,j] lt 2.3d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
;       if (fitwlen[i,j] gt 2.15d-6 and fitwlen[i,j] lt 2.3d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 4) then begin
      if (wlen[i,j] gt 2.3d-6 and wlen[i,j] lt 2.4d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
;       if (fitwlen[i,j] gt 2.3d-6 and fitwlen[i,j] lt 2.4d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

    if (j eq 5) then begin
      if (wlen[i,j] gt 2.35d-6 and wlen[i,j] lt 2.5d-6) then qcwlen[i,j] = '+1' else qcwlen[i,j] = '-1'
;       if (fitwlen[i,j] gt 2.35d-6 and fitwlen[i,j] lt 2.5d-6) then qcfitwlen[i,j] = '+1' else qcfitwlen[i,j] = '-1'
    endif

  endfor
endfor


;Phase Shift Quality Control
qcangdeg = strarr(3,6)
qcfitangdeg = strarr(3,6)
for i=0,2 do begin
  for j=0,5 do begin

    if (i eq 0) then begin
      if (angdeg[i,j] ge 65.d0 and angdeg[i,j] le 85.d0) then qcangdeg[i,j] = '+1' else qcangdeg[i,j] = '-1'
;       if (fitangdeg[i,j] ge 65.d0 and fitangdeg[i,j] le 85.d0) then qcfitangdeg[i,j] = '+1' else qcfitangdeg[i,j] = '-1'
    endif

    if (i eq 1) then begin
      if (angdeg[i,j] gt 155.d0 and angdeg[i,j] lt 205.d0) then qcangdeg[i,j] = '+1' else qcangdeg[i,j] = '-1'
;       if (fitangdeg[i,j] gt 155.d0 and fitangdeg[i,j] lt 175.d0) then qcfitangdeg[i,j] = '+1' else qcfitangdeg[i,j] = '-1'
    endif

    if (i eq 2) then begin
      if (angdeg[i,j] gt 215.d0 and angdeg[i,j] lt 250.d0) then qcangdeg[i,j] = '+1' else qcangdeg[i,j] = '-1'
;       if (fitangdeg[i,j] gt 215.d0 and fitangdeg[i,j] lt 250.d0) then qcfitangdeg[i,j] = '+1' else qcfitangdeg[i,j] = '-1'
    endif

  endfor
endfor



;output

openw, lun, resultpath+'FSUA_SkyCalibration_ABCD.dat', width=1400, /get_lun
;   printf, lun, 'Fit Fringe values'
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, fitphserr[i,j], format='(f15.10)'
;     endfor
;   endfor
;   printf, lun, ''
;   printf, lun, 'Classical Values (CC)'
  for i=0,5 do begin
    for j=0,3 do begin
      printf, lun, phserr[j,i], format='(f15.10)'
    endfor
  endfor
close, lun
free_lun, lun

openw, lun, resultpath+'FSUA_SkyCalibration_Angles.dat', width=1400, /get_lun
;   printf, lun, 'Fit Fringe values'
;   for i=0,5 do printf, lun, fitangval[*,i], format='(4f16.10)'
;   ;for i=0,5 do printf, lun, qcfitangval[*,i], format='(3a10)'
;   printf, lun, ''
;   printf, lun, 'Classical Values (CC)'
  for i=0,5 do printf, lun, angdeg[*,i], format='(3f15.10)'
  for i=0,5 do printf, lun, qcangdeg[*,i], format='(3a10)'
close, lun
free_lun, lun

openw, lun, resultpath+'FSUA_SkyCalibration_Lambda.dat', width=1400, /get_lun
;   printf, lun, 'Fit Fringe values'
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, fitwlen[j,i], format='(e20.10)'
;     endfor
;   endfor
;   printf, lun, ''
;   printf, lun, 'Classical Values (FFT)'
  for i=0,5 do begin
    for j=0,3 do begin
      printf, lun, wlen[j,i], format='(e20.10)'
    endfor
  endfor
close, lun
free_lun, lun

openw, lun, resultpath+'FSUA_Calibration_Waves.dat', width=1400, /get_lun
;   printf, lun, 'Fit Fringe values'
;   for i=0,5 do begin
;     for j=0,3 do begin
;       printf, lun, fitwlen[j,i], qcfitwlen[j,i], format='(e20.10, a5)'
;     endfor
;   endfor
;   printf, lun, ''
;   printf, lun, 'Classical Values (FFT)'
  for i=0,5 do begin
    for j=0,3 do begin
      printf, lun, wlen[j,i], qcwlen[j,i], format='(e20.10, a5)'
    endfor
  endfor
close, lun
free_lun, lun



print, ''
print, 'Wavelength and Phase values from fitting the Fringe package'
print, '         A                     B                     C                     D'
for i=0,5 do print, string(wlen[*,i])+'   '+string(qcfitwlen[*,i]) 
print, ''
for i=0,5 do print, string(angdeg[*,i]);+'   '+string(qcfitangval[*,i]) 

; print, ''
; print, 'Wavelength and Phase values from FFT and CC'
; print, '         A                     B                     C                     D'
; for i=0,5 do print, string(wlen[*,i])+'   '+string(qcwlen[*,i]) 
; print, ''
; for i=0,5 do print, string(angval[*,i])+'   '+string(qcangval[*,i]) 


stop
return
end
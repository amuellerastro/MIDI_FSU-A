; /home/amueller/work/MIDI_FSUA/DualFeed/prpipe/xpprojBL.m:      clat = cos(lat_vlti);
; /home/amueller/work/MIDI_FSUA/DualFeed/prpipe/xpprojBL.asv:      clat = cos(lat_vlti);
; /home/amueller/work/MIDI_FSUA/MIDI_FSUA_COMMP88/24Psc/24Psc/Sci_BinaryVis.pro:  clat = cos(lat_vlti)
; /home/amueller/work/MIDI_FSUA/MIDI_FSUA_COMMP88/24Psc/Cal/Cal_Vis.pro:  clat = cos(lat_vlti)
; /home/amueller/work/MIDI_FSUA/Pipeline/Conditions/get_V_projbaseline.pro:  clat = cos(lat_vlti)
; /home/amueller/work/MIDI_FSUA/Pipeline/ImportantSubroutines/get_baseline.pro:  clat = cos(lat_vlti)
; /home/amueller/work/MIDI_FSUA/Pipeline/V_FSUA/vsci.pro:  clat = cos(lat_vlti)
; /home/amueller/work/MIDI_FSUA/Pipeline/V_FSUA/vcal.pro:  clat = cos(lat_vlti)
; /home/amueller/work/MIDI_FSUA/WaterVapor/Akeson/calc_delay_array.pro:    clat = cos(lat_vlti)
; /home/amueller/work/MIDI_FSUA/WaterVapor/Akeson/calc_delay.pro:    clat = cos(lat_vlti)


@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@strnumber.pro
@readfits.pro
@sxpar.pro
@get_eso_keyword.pro
@caldat.pro
@ten.pro
@calcpos.pro
@addpm.pro
@get_visibility.pro
@int_tabulated.pro
@uniq.pro

pro vcal

readcol, 'FSUA_FringeScans.txt', obsid, scan, bg, ff1, ff2, starflag, diam, format='a,a,a,a,a,a,d', /silent

;-----------------------------

obsnum = ''
; read, 'Enter ObsID: ', obsnum
obsnum = 'FSUAscans_20130706' 

; wl = 2.25d-6
wl = [2.1915746d-06, 2.0002667d-06, 2.0959865d-06, 2.2190671d-06, 2.3369855d-06, 2.4290174d-06]	;actual median combined values from all 16 scan files of HD155826 using calibrate_all2ndScan.pro


;==================================================================

idx = where(obsid eq obsnum)

  obsid = obsid[idx]
  scan = scan[idx]
  bg = bg[idx]
  ff1 = ff1[idx]
  ff2 = ff2[idx]
  starflag = starflag[idx]
  diam = diam[idx]

idx = where(starflag eq 'Cal')

  obsid = obsid[idx]
  scan = scan[idx]
  bg = bg[idx]
  ff1 = ff1[idx]
  ff2 = ff2[idx]
  starflag = starflag[idx]
  diam = diam[idx]

nfiles = n_elements(obsid)
target = strarr(nfiles)
vis2 = dblarr(nfiles) & bl = vis2 & padeg = vis2 & u = vis2 & v = vis2 & w = vis2 & ha = vis2
vis2_px3 = dblarr(nfiles) & vis2_all = dblarr(nfiles, 6)
vis2_theoWhite = dblarr(nfiles)
utc = dblarr(nfiles)
fileid = strarr(nfiles)

for i=0,nfiles-1 do begin

  file = obsid[i]+'/'+scan[i]+'.fits'

  ;remove full path in front of file
  pos1 = strpos(file, 'FSUA_SECOND_FRINGE')
  pos2 = strpos(file, '.fits')
  file = strmid(file, pos1, pos2-pos1+5)

  ;extract file ID
  pos1 = strpos(file, 'SCAN_')
  pos2 = strpos(file, '.fits')
  fileid[i] = strmid(file, pos1+5, pos2-pos1-5)

  ;extract keywords present in OBJ_SCAN files

  dummy = readfits(obsid[i]+'/'+file, header, /silent)

  jd = double(get_eso_keyword(header, 'MJD-OBS'))+0.5d0+2400000.d0
;   object = strtrim(sxpar(header,'OBJECT'),2)
;   pos3 = strpos(object,'_')
;   id1 = strmid(object,0,pos3)
;   id2 = strmid(object,pos3+1,strlen(object)-pos3)
;   dummy = [id1,id2]
;   object = strjoin(dummy)
  tmp = get_eso_keyword(header, 'HIERARCH ESO OBS TARG NAME')
  target[i] = tmp

  tel1 = get_eso_keyword(header,'HIERARCH ESO ISS CONF T1NAME')	;id of 1st telescope
  tel2 = get_eso_keyword(header,'HIERARCH ESO ISS CONF T1NAME')	;id of 2nd telescope
  station1 = get_eso_keyword(header,'HIERARCH ESO ISS CONF STATION1')
  station2 = get_eso_keyword(header,'HIERARCH ESO ISS CONF STATION2')

  radeg = double(strtrim(sxpar(header,'RA'),2))	;degrees
  dec = double(strtrim(sxpar(header,'DEC'),2))	;degrees
  lstsec = double(strtrim(sxpar(header,'LST'),2))	;seconds

  ;get the middle of the observation
  tmid = double(get_eso_keyword(header, 'EXPTIME'))/2.d0
  lstsec = lstsec+tmid
;   if (lstsec gt 86400.d0) then lstsec = lstsec-86400.d0	;needed?
  utc[i] = double((get_eso_keyword(header, 'UTC')+tmid)/3600.d0)

  t1x = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T1X'))
  t1y = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T1Y'))
  t1z = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T1Z'))
  t2x = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T2X'))
  t2y = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T2Y'))
  t2z = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T2Z'))
  obs_long = double(get_eso_keyword(header,'HIERARCH ESO ISS GEOLON'))	;longitude VLTI
  lat_vlti = double(get_eso_keyword(header,'HIERARCH ESO ISS GEOLAT'))	;latitude VLTI

  dec = dec*!dtor
  lst = (lstsec)*15.d0/(3600.d0*!radeg)
  ra = radeg*!dtor
  ha[i] = lst - ra

  ; if (ha lt 0.) then ha = ha + 360.0d0
;   if (ha gt 180.d0) then ha = ha - 360.0d0
;   if (ha[i] gt 2.*!DPI or ha[i] lt -2.*!DPI) then stop	;something to do here?
  if (ha[i] lt -!DPI) then ha[i] = ha[i]+2.*!DPI
  if (ha[i] gt !DPI) then ha[i] = ha[i]-2.*!DPI

  x = t1x - t2x
  y = t1y - t2y
  z = t1z - t2z

  lat_vlti  = lat_vlti/!radeg

;   ihar = ha*!dtor 	;/ 57.29577951d0
;   decr = dec*!dtor	; / 57.29577951d0

  cha = cos(ha[i])
  sha = sin(ha[i])
  cdec = cos(dec)
  sdec = sin(dec)
  clat = cos(lat_vlti)
  slat = sin(lat_vlti)

  ;1st baseline T1 - T2

  ; compute u,v,w
  u[i] = x*cha - y*slat*sha + z*clat*sha
  v[i] = x*sha*sdec + y*(slat*cha*sdec + clat*cdec) - $
	  z*(clat*cha*sdec - slat*cdec)
  w[i] = -x*sha*cdec - y*(slat*cha*cdec - clat*sdec) + $
		z*(clat*cha*cdec + slat*sdec)

  SinAlt =   (clat*cha*cdec + slat*sdec)
  alt = asin(SinAlt) ; <<<  in Radians!
  pa = atan(y, x) ; also in radians
  para = atan(sha, (slat/clat)*cdec-sdec*cha)

  alt = alt*!radeg ; to degrees

  bl[i] = sqrt((u[i])^2.+(v[i])^2.)
  padeg[i] = (atan(u[i], v[i]))*!radeg

  if (padeg[i] lt 0. and padeg[i] ge -180.d0) then padeg[i] = padeg[i] + 180.d0 else stop

  vis2_px3[i] = (get_visibility(bl[i], wl[3], diam[i]/1.d3, 0.))^2.	;for quick check compute vis^2 using middle of K-band, assuming 3rd pixel of FSUA

  for j=0,5 do vis2_all[i,j] = (get_visibility(bl[i], wl[j], diam[i]/1.d3, 0.))^2.

  dum = dblarr(5)
  dum[*] = 1.d0
  vis2_theoWhite[i] = int_tabulated(wl[1:*],vis2_all[i,1:*]*wl[1:*]^2.,/double)/int_tabulated(wl[1:*], wl[1:*]^2)	;theoretical Visibility due to bandpass smearing, Kervella+2004

endfor
print, ''
print, '**************************************************************************************************************************'
print, 'Target     FileID             UTC             V^2           proj.BL           PA             U               V              HA'
for i=0,nfiles-1 do print, target[i], '  ', fileid[i], utc[i], vis2_px3[i], bl[i], padeg[i], u[i], v[i], ha[i]/!dtor/15.

print, ''
print, 'theoretical V^2 for White Light pixel'
print, vis2_theoWhite

print, '**************************************************************************************************************************'
print, ''

stop
end
function get_baseline, path, id

  if (strmatch(id, '*fits') eq 1) then file = file_search(path+id,count=nfiles) $
    else  file = file_search(path+id+'*.fits',count=nfiles)

  ;get time in seconds in the middle of the observation since start
;   dummy1 = mrdfits(file[0], 9, /silent)
;   dummy2 = mrdfits(file[nfiles-1], 9, /silent)
;   tmid = 86400.d0*(abs(dummy2[n_elements(dummy2.time)-1].time - dummy1[0].time))/2.d0

  file = file[0]
  dummy = readfits(file, header, /silent)
  catg = get_eso_keyword(header,'HIERARCH ESO DPR CATG')
  t1xdummy = get_eso_keyword(header,'HIERARCH ESO ISS CONF T1X')
  mode = get_eso_keyword(header,'HIERARCH ESO DET NRTS MODE')

  object = strtrim(sxpar(header,'OBJECT'),2)
  pos3 = strpos(object,'_')
  id1 = strmid(object,0,pos3)
  id2 = strmid(object,pos3+1,strlen(object)-pos3)
  dummy = [id1,id2]
  object = strjoin(dummy)


  tel1 = get_eso_keyword(header,'HIERARCH ESO ISS CONF T1NAME')	;id of 1st telescope
  tel2 = get_eso_keyword(header,'HIERARCH ESO ISS CONF T1NAME')	;id of 2nd telescope
  station1 = get_eso_keyword(header,'HIERARCH ESO ISS CONF STATION1')
  station2 = get_eso_keyword(header,'HIERARCH ESO ISS CONF STATION2')

  radeg = double(strtrim(sxpar(header,'RA'),2))	;degrees
  dec = double(strtrim(sxpar(header,'DEC'),2))	;degrees
  lstsec = double(strtrim(sxpar(header,'LST'),2))	;seconds
  t1x = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T1X'))
  t1y = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T1Y'))
  t1z = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T1Z'))
  t2x = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T2X'))
  t2y = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T2Y'))
  t2z = double(get_eso_keyword(header,'HIERARCH ESO ISS CONF T2Z'))

  dec = dec*!dtor

;   lst = (lstsec+tmid)*15.d0/3600.d0
  lst = (lstsec)*15.d0/(3600.d0*!radeg)
  obs_long = double(get_eso_keyword(header,'HIERARCH ESO ISS GEOLON'))	;longitude VLTI
  lat_vlti = double(get_eso_keyword(header,'HIERARCH ESO ISS GEOLAT'))	;latitude VLTI
  ;lat_vlti  = -0.429838786d0
  ra = radeg*!dtor
  ha = lst - ra

  ; if (ha lt 0.) then ha = ha + 360.0d0
;   if (ha gt 180.d0) then ha = ha - 360.0d0
  if (ha lt -!DPI) then ha = ha+2.*!DPI
  if (ha gt !DPI) then ha = ha-2.*!DPI

  x = t1x - t2x
  y = t1y - t2y
  z = t1z - t2z

  lat_vlti  = lat_vlti/!radeg

;   ihar = ha*!dtor 	;/ 57.29577951d0
;   decr = dec*!dtor	; / 57.29577951d0

  cha = cos(ha)
  sha = sin(ha)
  cdec = cos(dec)
  sdec = sin(dec)
  clat = cos(lat_vlti)
  slat = sin(lat_vlti)

  ;1st baseline T1 - T2

  ; compute u,v,w
  u = x*cha - y*slat*sha + z*clat*sha
  v = x*sha*sdec + y*(slat*cha*sdec + clat*cdec) - $
	  z*(clat*cha*sdec - slat*cdec)
  w = -x*sha*cdec - y*(slat*cha*cdec - clat*sdec) + $
		z*(clat*cha*cdec + slat*sdec)


  SinAlt =   (clat*cha*cdec + slat*sdec)
  alt = asin(SinAlt) ; <<<  in Radians!
  pa = atan(y, x) ; also in radians
  para = atan(sha, (slat/clat)*cdec-sdec*cha)


  ;getuvwa,ha, dec, xyz, uvw,para=para,alt=a
  alt = alt*!radeg ; to degrees

  bl = sqrt((u)^2.+(v)^2.)
  padeg = (atan(u, v))*!radeg
;   if (padeg lt 0.) then padeg = padeg + 360.d0
  if (padeg lt 0. and padeg ge -180.d0) then padeg = padeg + 180.d0; else stop


  return, {bl:bl, padeg:padeg, u:u, v:v, w:w, ha:ha}

end
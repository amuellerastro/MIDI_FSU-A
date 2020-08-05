@counter.pro
@fifteenb.pro

pro get_V_projbaseline

commnum = ''
read, 'Enter Observation ID: ', commnum


midipath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/MIDIdata/'
fsupath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/FSUAdata/'
path = '/home/amueller/work/MIDI_FSUA/MIDI_FSUA_'+commnum+'/'


readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', cn, id, scical, nf, mf, photA, photB, pf, maskf, phot, format = 'a,a,a,d,a,a,a,a,a,a', /silent
idx = where(cn eq commnum)
id = id[idx]
nf = nf[idx]
scical = scical[idx]
mf = mf[idx]
pf = pf[idx]
maskf = maskf[idx]
phot = phot[idx]
photA = photA[idx]
photB = photB[idx]

outfile = 'projBaseline_'+commnum+'.txt'

openw, lun, path+outfile, width=1400, /get_lun;, /append
printf, lun, '                        File              Object  Station        U         V        W   proj. BL       PA'
printf, lun, '                                                               [m]       [m]       [m]       [m]     [deg]'

for i=0,n_elements(mf)-1 do begin

  counter, i+1, n_elements(mf), 'Processing observation '

  file = file_search(midipath+mf[i]+'*.fits',count=nfiles)

  ;get time in seconds in the middle of the observation since start
  dummy1 = mrdfits(file[0], 9, /silent)
  dummy2 = mrdfits(file[nfiles-1], 9, /silent)
  tmid = 86400.d0*(abs(dummy2[n_elements(dummy2.time)-1].time - dummy1[0].time))/2.d0

  file = file[0]

  dummy = readfits(file,header, /silent)
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

;   lst = (lstsec+tmid)*360.d0/86400.d0
  lst = (lstsec)*15.d0/(3600.d0*!radeg)
  obs_long = double(get_eso_keyword(header,'HIERARCH ESO ISS GEOLON'))	;longitude VLTI
  lat_vlti = double(get_eso_keyword(header,'HIERARCH ESO ISS GEOLAT'))	;latitude VLTI

  ra = radeg*!dtor
  ha = lst - ra

  ; if (ha lt 0.) then ha = ha + 360.0d0
;   if (ha gt 180.) then ha = ha - 360.0d0
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

  bl = sqrt((u)^2+(v)^2)
  padeg = (atan(u, v))*!radeg
;  if (pa lt 0.) then pa = pa + 180.d0
;   if (padeg lt 0.) then padeg = padeg + 360.d0
  if (padeg lt 0. and padeg ge -180.d0) then padeg = padeg + 180.d0; else stop


  dum = '-'
  printf, lun, mf[i],id[i],station1,dum,station2,u,v,w,bl,padeg,format='(a28,a20,a5,a1,a2,5f10.2)'

  ;     PRINT,'lst ',lst,' ',lst/15.d0
  ;     PRINT,'HA [deg/hr] ',ha,' ',ha/15.
  ;     PRINT,'xyz [m] ',x,y,z,sqrt((x)^2+(y)^2+z^2),format='(a,4f10.5)'
  ;     PRINT,'uvw [m] ',u, v, w,sqrt((u)^2+(v)^2+(w)^2),sqrt((u)^2+(v)^2),format='(a, 5f10.5)'
  ;     PRINT,'proj. BL [m] ',bl
  ;     PRINT,'para angle [deg] ',para*!radeg
  ;     PRINT,'PA [deg]',pa
  ;     PRINT,'Alt. [deg] ',alt
  ;     PRINT,'airmass ',1./sin(alt/!radeg)

endfor



close, lun
free_lun, lun


;compute K and N band Visibilities
readcol,path+outfile,ofile,object,station,u,v,w,pbl,pa,format='a,a,a,d,d,d,d,d', /silent
readcol,path+'stars_'+commnum+'.txt',id,ra,dec,vmag,jmag,hmag,kmag,f10,diam,format='a,a,a,d,d,d,d,d,d', /silent


conv = 180.d0*3600.d0/!DPI
wlfsu = 2.2*1.0D-6	;K-band
wlmidi = 10.0*1.0D-6	;N-band


openw, lun, path+'Vis_projBaseline_COMM'+commnum+'.txt', width=1400, /get_lun
printf, lun, '                        File              Object Station  proj. BL       PA    V(2.2um)  V(10um)'
printf, lun, '                                                              [m]     [deg]'

for i=0,n_elements(object)-1 do begin

  index = where(object[i] eq id)
  if (index(0) ne -1) then begin

  diameter = diam[index]

  if (diameter ne 99) then begin

    diameter = diameter/1000.

    sffsu = pbl[i]/wlfsu/conv	;spatial frequency
    sfmidi = pbl[i]/wlmidi/conv
    visfsu = abs(2.d0*BESELJ(!DPI*sffsu*diameter,1,/DOUBLE)/(!DPI*sffsu*diameter))
    vismidi = abs(2.d0*BESELJ(!DPI*sfmidi*diameter,1,/DOUBLE)/(!DPI*sfmidi*diameter))

    v2fsu = visfsu*visfsu
    v2midi = vismidi*vismidi

  endif else begin

    visfsu = !values.d_nan
    vismidi = !values.d_nan
    v2fsu = !values.d_nan
    v2midi = !values.d_nan

  endelse

  printf, lun,ofile[i],object[i],station[i],pbl[i],pa[i],visfsu,vismidi,format='(a28,a20,a7,4f10.2)'

  endif

endfor

close, lun
free_lun, lun

end
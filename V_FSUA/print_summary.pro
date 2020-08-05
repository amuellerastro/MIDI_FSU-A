pro print_summary

obsnum = 'FSUAscans_P94KOEHLERD0I1'
; obsnum = 'FSUAscans_P94KOEHLERH0I1'

dir = 'Vis_'+obsnum+'_5pixel'
visdir = dir+'/TFcal_V2sci_InterpDeg2/'

;=============================================================================

readcol, 'FSUA_FringeScans.txt', obsid, scan, bg, ff1, ff2, starflag, diam, format='a,a,a,a,a,a,d', /silent

idx = where(obsid eq obsnum and starflag eq 'Sci')

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

nvis = n_elements(obsid)

file = file_search(visdir+'*_Sci_*.txt', count=nfiles)


pbl = dblarr(nvis)
pa = dblarr(nvis)
jd = dblarr(nvis)
v0 = dblarr(nvis) & v0e = v0

for i=0,nvis-1 do begin

  idx = strmatch(file, '*'+strmid(scan[i], 24,8)+'*')
  idx1 = where(idx eq 1)

  readcol, file[idx1], t1,t2,t3,t4,t5,t6,t7,t8,t9, format='d,d,d,d,d,d,d,d,d', /silent, skipline=14, numline=1

  v0[i] = t4
  v0e[i] = t5

  ;------------------------------------------

  st = mrdfits(obsnum+'/'+scan[i]+'.fits', 0, hdr, /silent)

  pbl1 = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PBL12 START'))
  pbl2 = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PBL12 END'))
  pbl[i] = mean([pbl1, pbl2])

  pa1 = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PBLA12 START'))
  pa2 = double(get_eso_keyword(hdr, 'HIERARCH ESO ISS PBLA12 END'))
  pa[i] = mean([pa1, pa2])

  mjd = double(get_eso_keyword(hdr, 'MJD-OBS'))
  inttime = ((double(get_eso_keyword(hdr, 'EXPTIME')))/86400.d0)/2.d0

  jd[i] = mjd+inttime+2400000.5d0


endfor


openw, lun, 'summary_'+obsnum+'.txt', width=1400, /get_lun

  printf, lun, '                JD   BL / m  PA / deg  V^2    dV^2'
  for i=0,nvis-1 do printf, lun, jd[i], pbl[i], pa[i], v0[i], v0e[i], format='(f18.8, 2f9.2, 2f7.3)'

close, lun
free_lun, lun


stop
end
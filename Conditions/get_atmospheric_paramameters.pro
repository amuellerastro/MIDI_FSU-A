@counter.pro
@fifteenb.pro
@get_eso_keyword.pro


pro get_atmospheric_paramameters

commnum = ''
read, 'Enter Observation ID: ', commnum

path = '/home/amueller/work/MIDI_FSUA/MIDI_FSUA_'+commnum+'/'
midipath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/MIDIdata/'
fsupath = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/FSUAdata/'

readcol, '/home/amueller/work/MIDI_FSUA/observation.txt', cn, id, scical, nf, mf, photA, photB, pf, maskf, phot, format = 'a,a,a,d,a,a,a,a,a,a', /silent
idx = where(cn eq commnum)
id = id[idx]
nf = nf[idx]
mf = mf[idx]
scical = scical[idx]
pf = pf[idx]
maskf = maskf[idx]
phot = phot[idx]
photA = photA[idx]
photB = photB[idx]

openw,lun,path+'conditions_'+commnum+'.txt', width=1400,/get_lun
printf, lun, '                        File              Star         AM  seeing [arcsec] tau0 [ms]   rel. Humid. [%]'

for i=0,n_elements(id)-1 do begin

  ;counter, i+1, n_elements(id), 'Processing observation '

  file = file_search(midipath+mf[i]+'*.fits',count=nfiles)
  file = file[0]

  dummy = readfits(file,hdr,/silent)
  star = sxpar(hdr, 'OBJECT')
  airmass_start = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AIRM START'))
  airmass_end = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AIRM END'))
  seeing_start = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI FWHM START'))
  seeing_end = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI FWHM END'))
  tc_start = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI TAU0 START'))	;in [sec]
  tc_end = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI TAU0 END'))	;in [sec]
  rh = double(get_eso_keyword(hdr,'HIERARCH ESO ISS AMBI RHUM'))	;relative humidity in %

  if (seeing_start ge 0. and seeing_end ge 0.) then seeing = (double(seeing_start)+double(seeing_end))/2.
  if (airmass_start ge 0. and airmass_end ge 0.) then airmass = (double(airmass_start)+double(airmass_end))/2.
  if (tc_start ge 0. and tc_end ge 0.) then tc = ((double(tc_start)+double(tc_end))/2.)*1000.

  if (seeing_start le 0. and seeing_end ge 0.) then seeing = seeing_end 
  if (seeing_start ge 0. and seeing_end le 0.) then seeing = seeing_start
  if (seeing_start ge 0. and seeing_end ge 0.) then seeing = mean([seeing_start,seeing_end])
  if (seeing_start le 0. and seeing_end le 0.) then seeing = 999.

  if (seeing eq 999.) then begin
    dummy = readfits(fsupath+pf[i], hdrfsu, /silent)
    seeing_start = get_eso_keyword(hdrfsu,'HIERARCH ESO ISS AMBI FWHM START')
    seeing_end = get_eso_keyword(hdrfsu,'HIERARCH ESO ISS AMBI FWHM END')
    if (seeing_start le 0. and seeing_end ge 0.) then seeing = seeing_end 
    if (seeing_start ge 0. and seeing_end le 0.) then seeing = seeing_start
    if (seeing_start ge 0. and seeing_end ge 0.) then seeing = mean([seeing_start,seeing_end])
    if (seeing_start le 0. and seeing_end le 0.) then seeing = 999.
  endif


  if (airmass_start le 0. and airmass_end ge 0.) then airmass = airmass_end 
  if (airmass_start ge 0. and airmass_end le 0.) then airmass = airmass_start
  if (airmass_start ge 0. and airmass_end ge 0.) then airmass = mean([airmass_start,airmass_end])
  if (airmass_start le 0. and airmass_end le 0.) then airmass = 999.

  if (tc_start le 0. and tc_end ge 0.) then tc = tc_end*1000.
  if (tc_start ge 0. and tc_end le 0.) then tc = tc_start*1000.
  if (tc_start ge 0. and tc_end ge 0.) then tc = mean([tc_start,tc_end])*1000.
  if (tc_start le 0. and tc_end le 0.) then tc = 999.

  printf, lun, mf[i], id[i], airmass, seeing, tc, rh, format='(a28,a20,f11.3,2f11.2,I10)'

endfor

close, lun
free_lun, lun

;stop
end
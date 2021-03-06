@readfits.pro
@sxpar.pro
@gettok.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mwrfits.pro
@fxaddpar.pro
@fxmove.pro
@mrd_skip.pro
@match.pro
@mrd_struct.pro
@detabify.pro
@fxparpos.pro
@sxaddhist.pro


pro merge_MIDI_FSUA_observations


;works only for two observations for the moment
;===================================================
commnum = 'P94RATZKA'
indmidi = ['MIDI.2015-01-21T09:16:32', 'MIDI.2015-01-21T09:23:49']
indfsu = ['FSUA_OBS_021_0015.fits', 'FSUA_OBS_021_0016.fits']
finalmidi = indmidi[0]
finalfsu = indfsu[0]
; ;===================================================
; commnum = 'P94RATZKA'
; indmidi = ['MIDI.2015-01-21T08:42:10', 'MIDI.2015-01-21T08:47:55']
; indfsu = ['FSUA_OBS_021_0012.fits', 'FSUA_OBS_021_0013.fits']
; finalmidi = indmidi[0]
; finalfsu = indfsu[0]
; ;===================================================
; commnum = 'P94RATZKA'
; indmidi = ['MIDI.2015-01-21T08:02:01', 'MIDI.2015-01-21T08:06:33']
; indfsu = ['FSUA_OBS_021_0009.fits', 'FSUA_OBS_021_0010.fits']
; finalmidi = indmidi[0]
; finalfsu = indfsu[0]
; ;===================================================
; commnum = 'P94RATZKA'
; indmidi = ['MIDI.2015-01-21T07:26:55', 'MIDI.2015-01-21T07:36:27']
; indfsu = ['FSUA_OBS_021_0006.fits', 'FSUA_OBS_021_0007.fits']
; finalmidi = indmidi[0]
; finalfsu = indfsu[0]
; ;===================================================
; commnum = 'P94RATZKA'
; indmidi = ['MIDI.2015-01-21T06:29:50', 'MIDI.2015-01-21T06:45:56']
; indfsu = ['FSUA_OBS_021_0002.fits', 'FSUA_OBS_021_0004.fits']
; finalmidi = indmidi[0]
; finalfsu = indfsu[0]
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T07:30:19', 'MIDI.2015-01-07T07:34:08']
; finalmidi = indmidi[0]
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T06:48:19', 'MIDI.2015-01-07T07:03:12']
; finalmidi = indmidi[0]
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T06:27:46', 'MIDI.2015-01-07T06:33:42']
; finalmidi = indmidi[0]
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T05:41:10', 'MIDI.2015-01-07T05:45:23']
; finalmidi = indmidi[0]
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T05:02:02', 'MIDI.2015-01-07T05:16:44']
; finalmidi = indmidi[0]
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T04:41:00', 'MIDI.2015-01-07T04:46:32']
; finalmidi = indmidi[0]
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T03:58:54', 'MIDI.2015-01-07T04:03:13']
; finalmidi = 'MIDI.2015-01-07T03:58:54'
; ;===================================================
; commnum = 'P94KRAUS'
; indmidi = ['MIDI.2015-01-07T01:17:51', 'MIDI.2015-01-07T01:49:40']
; finalmidi = 'MIDI.2015-01-07T01:17:51'
; ;===================================================
; commnum = 'P93KRAUSU23'
; 
; ;individual files
; indfsu = ['FSUA_OBS_104_0018.fits', 'FSUA_OBS_104_0019.fits']
; indmidi = ['MIDI.2014-04-14T09:16:52', 'MIDI.2014-04-14T09:29:20']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_104_0018.fits'
; finalmidi = 'MIDI.2014-04-14T09:16:52'
;===================================================
; ;===================================================
; commnum = 'P93KRAUSU23'
; 
; ;individual files
; indfsu = ['FSUA_OBS_104_0013.fits', 'FSUA_OBS_104_0014.fits']
; indmidi = ['MIDI.2014-04-14T07:31:36', 'MIDI.2014-04-14T07:35:35']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_104_0013.fits'
; finalmidi = 'MIDI.2014-04-14T07:31:36'
; ;===================================================
; ;===================================================
; commnum = 'P93KRAUSU23'
; 
; ;individual files
; indfsu = ['FSUA_OBS_104_0010.fits', 'FSUA_OBS_104_0011.fits']
; indmidi = ['MIDI.2014-04-14T06:41:43', 'MIDI.2014-04-14T06:52:28']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_104_0010.fits'
; finalmidi = 'MIDI.2014-04-14T06:41:43'
; ;===================================================
; ;===================================================
; commnum = 'P93KRAUSU23'
; 
; ;individual files
; indfsu = ['FSUA_OBS_104_0002.fits', 'FSUA_OBS_104_0003.fits']
; indmidi = ['MIDI.2014-04-14T03:19:08', 'MIDI.2014-04-14T03:25:38']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_104_0002.fits'
; finalmidi = 'MIDI.2014-04-14T03:19:08'
; ;===================================================
; ;===================================================
; commnum = 'P93KRAUSU14'
; 
; ;individual files
; indfsu = ['FSUA_OBS_103_0012.fits', 'FSUA_OBS_103_0013.fits']
; indmidi = ['MIDI.2014-04-13T08:44:00', 'MIDI.2014-04-13T08:56:23']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_103_0012.fits'
; finalmidi = 'MIDI.2014-04-13T08:44:00'
; ;===================================================
; ;===================================================
; commnum = 'P93KRAUSU14'
; 
; ;individual files
; indfsu = ['FSUA_OBS_103_0009.fits', 'FSUA_OBS_103_0010.fits']
; indmidi = ['MIDI.2014-04-13T07:46:00', 'MIDI.2014-04-13T07:57:00']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_103_0009.fits'
; finalmidi = 'MIDI.2014-04-13T07:46:00'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_067_0006.fits', 'FSUA_OBS_067_0007.fits']
; indmidi = ['MIDI.2014-03-08T03:23:23', 'MIDI.2014-03-08T03:27:21']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_067_0006.fits'
; finalmidi = 'MIDI.2014-03-08T03:23:23'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0016.fits', 'FSUA_OBS_066_0020.fits']
; indmidi = ['MIDI.2014-03-07T04:13:58', 'MIDI.2014-03-07T04:38:29']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0016.fits'
; finalmidi = 'MIDI.2014-03-07T04:13:58'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0016.fits', 'FSUA_OBS_066_0019.fits']
; indmidi = ['MIDI.2014-03-07T04:13:58', 'MIDI.2014-03-07T04:25:56']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0016.fits'
; finalmidi = 'MIDI.2014-03-07T04:13:58'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0016.fits', 'FSUA_OBS_066_0018.fits']
; indmidi = ['MIDI.2014-03-07T04:13:58', 'MIDI.2014-03-07T04:21:55']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0016.fits'
; finalmidi = 'MIDI.2014-03-07T04:13:58'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0016.fits', 'FSUA_OBS_066_0017.fits']
; indmidi = ['MIDI.2014-03-07T04:13:58', 'MIDI.2014-03-07T04:17:55']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0016.fits'
; finalmidi = 'MIDI.2014-03-07T04:13:58'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0014.fits', 'FSUA_OBS_066_0015.fits']
; indmidi = ['MIDI.2014-03-07T03:50:28', 'MIDI.2014-03-07T03:54:17']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0014.fits'
; finalmidi = 'MIDI.2014-03-07T03:50:28'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0012.fits', 'FSUA_OBS_066_0013.fits']
; indmidi = ['MIDI.2014-03-07T03:30:49', 'MIDI.2014-03-07T03:34:39']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0012.fits'
; finalmidi = 'MIDI.2014-03-07T03:30:49'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0009.fits', 'FSUA_OBS_066_0010.fits']
; indmidi = ['MIDI.2014-03-07T02:42:11', 'MIDI.2014-03-07T02:48:31']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0009.fits'
; finalmidi = 'MIDI.2014-03-07T02:42:11'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0002.fits', 'FSUA_OBS_066_0008.fits']
; indmidi = ['MIDI.2014-03-07T01:37:15', 'MIDI.2014-03-07T02:01:49']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0002.fits'
; finalmidi = 'MIDI.2014-03-07T01:37:15'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0002.fits', 'FSUA_OBS_066_0007.fits']
; indmidi = ['MIDI.2014-03-07T01:37:15', 'MIDI.2014-03-07T01:57:54']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0002.fits'
; finalmidi = 'MIDI.2014-03-07T01:37:15'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0002.fits', 'FSUA_OBS_066_0006.fits']
; indmidi = ['MIDI.2014-03-07T01:37:15', 'MIDI.2014-03-07T01:54:05']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0002.fits'
; finalmidi = 'MIDI.2014-03-07T01:37:15'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0002.fits', 'FSUA_OBS_066_0005.fits']
; indmidi = ['MIDI.2014-03-07T01:37:15', 'MIDI.2014-03-07T01:49:39']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0002.fits'
; finalmidi = 'MIDI.2014-03-07T01:37:15'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0002.fits', 'FSUA_OBS_066_0004.fits']
; indmidi = ['MIDI.2014-03-07T01:37:15', 'MIDI.2014-03-07T01:45:29']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0002.fits'
; finalmidi = 'MIDI.2014-03-07T01:37:15'
; ;===================================================
; ;===================================================
; commnum = 'P92RIVI'
; 
; ;individual files
; indfsu = ['FSUA_OBS_066_0002.fits', 'FSUA_OBS_066_0003.fits']
; indmidi = ['MIDI.2014-03-07T01:37:15', 'MIDI.2014-03-07T01:41:16']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_066_0002.fits'
; finalmidi = 'MIDI.2014-03-07T01:37:15'
; ;===================================================
; ;===================================================
; commnum = 'P92KRAUS1'
; 
; ;individual files
; indfsu = ['FSUA_OBS_353_0010.fits', 'FSUA_OBS_353_0011.fits']
; indmidi = ['MIDI.2013-12-19T04:17:46', 'MIDI.2013-12-19T04:24:12']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_353_0010.fits'
; finalmidi = 'MIDI.2013-12-19T04:17:46'
; ;===================================================
; ;===================================================
; commnum = 'P92KRAUS2'
; 
; ;individual files
; indfsu = ['FSUA_OBS_355_0006.fits', 'FSUA_OBS_355_0007.fits']
; indmidi = ['MIDI.2013-12-21T02:14:16', 'MIDI.2013-12-21T02:28:21']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_355_0006.fits'
; finalmidi = 'MIDI.2013-12-21T02:14:16'
; ;===================================================
; ;===================================================
; commnum = 'P92KRAUS2'
; 
; ;individual files
; indfsu = ['FSUA_OBS_355_0004.fits', 'FSUA_OBS_355_0005.fits']
; indmidi = ['MIDI.2013-12-21T01:44:22', 'MIDI.2013-12-21T01:54:33']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_355_0004.fits'
; finalmidi = 'MIDI.2013-12-21T01:44:22'
; ;===================================================
; ;===================================================
; commnum = 'P92KRAUS1'
; 
; ;individual files
; indfsu = ['FSUA_OBS_353_0013.fits', 'FSUA_OBS_353_0014.fits']
; indmidi = ['MIDI.2013-12-19T05:01:37', 'MIDI.2013-12-19T05:13:22']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_353_0013.fits'
; finalmidi = 'MIDI.2013-12-19T05:01:37'
; ;===================================================
; ;===================================================
; commnum = 'P92KRAUS1'
; 
; ;individual files
; indfsu = ['FSUA_OBS_353_0004.fits', 'FSUA_OBS_353_0005.fits']
; indmidi = ['MIDI.2013-12-19T02:10:12', 'MIDI.2013-12-19T02:19:39']
; 
; ;final file names of the merged data files
; finalfsu = 'FSUA_OBS_353_0004.fits'
; finalmidi = 'MIDI.2013-12-19T02:10:12'
; ;===================================================



path = '/media/amueller/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/'

midipath = path+'MIDIdata/'
fsupath = path+'FSUAdata/'

backupmidi = midipath+'not_merged/'
spawn, 'mkdir '+backupmidi
backupfsu = fsupath+'not_merged/'
spawn, 'mkdir '+backupfsu


;backup the individual files before merging
nobs = n_elements(indfsu)
for i=0,nobs-1 do begin

  spawn, 'mv '+fsupath+indfsu[i]+' '+backupfsu
  if (i eq 0) then spawn, 'cp '+midipath+indmidi[i]+'* '+backupmidi $
    else spawn, 'mv '+midipath+indmidi[i]+'* '+backupmidi

endfor


;==========================================================================================

;FSU
;===

dum = readfits(backupfsu+indfsu[0], hdr1, /silent)
dum = readfits(backupfsu+indfsu[1], hdr2, /silent)

fsu_mjd1 = sxpar(hdr1,'MJD-OBS')
fsu_utc_dum = sxpar(hdr1,'DATE-OBS')
fsu_utc = dblarr(3)
fsu_utc[0] = double(strmid(fsu_utc_dum,11,2))	;UTC hrs
fsu_utc[1] = double(strmid(fsu_utc_dum,14,2))	;UTC min
fsu_utc[2] = double(strmid(fsu_utc_dum,17,strlen(fsu_utc_dum)-17))	;UTC sec
;recalculate the start MJD more accurate
jd1 = floor(fsu_mjd1) + fsu_utc[0]/24.d0 + fsu_utc[1]/(24.d0*60.d0) + fsu_utc[2]/(24.d0*3600.d0)

fsu_mjd2 = sxpar(hdr2,'MJD-OBS')
fsu_utc_dum = sxpar(hdr2,'DATE-OBS')
fsu_utc = dblarr(3)
fsu_utc[0] = double(strmid(fsu_utc_dum,11,2))	;UTC hrs
fsu_utc[1] = double(strmid(fsu_utc_dum,14,2))	;UTC min
fsu_utc[2] = double(strmid(fsu_utc_dum,17,strlen(fsu_utc_dum)-17))	;UTC sec
;recalculate the start MJD more accurate
jd2 = floor(fsu_mjd2) + fsu_utc[0]/24.d0 + fsu_utc[1]/(24.d0*60.d0) + fsu_utc[2]/(24.d0*3600.d0)

diff = (jd2 - jd1)*86400.d0*1.d6	;time in microsec when 2nd observation started

for j=0,11 do begin

  data1 = mrdfits(backupfsu+indfsu[0], j, hdr1, /silent)
  data2 = mrdfits(backupfsu+indfsu[1], j, hdr2, /silent)

  if (j ne 6 and j ne 7 and j ne 8 and j ne 11) then begin 

    mwrfits, data1, fsupath+finalfsu, hdr1, /silent

  endif else begin

    data2.time = data2.time+diff
    data = [data1,data2]

    mwrfits, data, fsupath+finalfsu, hdr1, /silent

  endelse

endfor

;==========================================================================================
; 
;MIDI
;====

;merging means basically just renaming the files as time stamps are in JD

;extract the last number
reffile = file_search(backupmidi+indmidi[0]+'*', count=ln)

;extract last frame number of reference file
refst9 = mrdfits(reffile[ln-1], 9, /silent)
refframe = refst9[n_elements(refst9.frame)-1].frame


file = file_search(backupmidi+indmidi[1]+'*', count=nfiles)
nframes = uintarr(nfiles)
for i=0,nfiles-1 do begin

  st9 = mrdfits(file[i], 9, /silent)
  nframes[i] = n_elements(st9.frame)

  if (i eq 0) then framenum = lindgen(nframes[i])+long(1)+refframe $
    else framenum = lindgen(nframes[i])+long(1)+refframe+long(total(nframes[0:i-1]))
  ;print, framenum[0], framenum[n_elements(framenum)-1]

  ;laufender index
  num = ln+i+1
  if (num lt 10) then num = '0'+strcompress(num, /rem) $
    else num = strcompress(num, /rem)

  ;write MIDI files with updated framenumbers
  for j=0,9 do begin

    data = mrdfits(file[i], j, hdr, /silent)

    if (j ne 9) then begin 

      mwrfits, data, midipath+finalmidi+'.000_'+num+'.fits', hdr, /silent

    endif else begin

      data.frame = framenum
      mwrfits, data, midipath+finalmidi+'.000_'+num+'.fits', hdr, /silent

    endelse

  endfor

endfor

end
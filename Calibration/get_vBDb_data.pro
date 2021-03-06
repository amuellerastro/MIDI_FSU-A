pro get_vBDb_data

cdb = vboekelbase()	;EWS shortcut to the database

id = ''
quest = 'y'

repeat begin

  read, 'ID: ', id

  nfound = cdb->sourceByName(id, sourceData, diamData, photData, specData)
  if (total(nfound) eq 0.) then begin
    print, ''
    print, 'Calibrator is not in vanBoekel Database. Stop.'
    print, ''
    stop
  endif

  print, ''
  print, 'ID     : ', id
  print, 'Diam   : ', diamData.diameter*1.d3
  print, 'Diamerr: ', diamData.diameter_err*1.d3
  print, 'Band   : ', photdata[11].band
  print, 'Flux   : ', photdata[11].flux

  print, ''
  read, 'Another Star? (y/n): ', quest

endrep until (quest ne 'y')

stop
end


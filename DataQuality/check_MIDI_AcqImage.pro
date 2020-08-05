;adopted from
;http://www.strw.leidenuniv.nl/~jaffe/ews/MIA+EWS-Manual/basic_tools.html
;http://www.strw.leidenuniv.nl/~jaffe/ews/MIA+EWS-Manual/photometry.html

pro check_MIDI_AcqImage

commnum = 'P92CHOQUET'

; path = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+commnum+'/MIDIdata/Acquisition_Images/'
path = '/home/amueller/Downloads/fsua/FSUA_2014_04_06/'

;   path = '/home/amueller/Downloads/superbright/'
;path = '/media/disk/MIDI_FSUA/COMM'+commnum+'/MIDIdata/PRIMA_STS_Acquisition_Images_20110829/'

print, ''
print, '***********************************'
print, '*Can only be started using MIA+EWS*'
print, '***********************************'
print, ''

;file selection
  ;cd, path

  print, ''
  print, 'Select ONE acquisition image'
  print, ''

  file = midigui(dir=path)


;unchopped data
  acq = oirgetdata(file)

  data1 = acq.data1
  data2 = acq.data2

  ;average all frames
    sz = size(data1)
    tmp1 = dblarr(sz[1], sz[2])
    tmp2 = dblarr(sz[1], sz[2])
    for xx=0L,sz[1]-1 do begin
      for yy=0L,sz[2]-1 do begin

      tmp1[xx,yy] = mean(data1[xx,yy,*])
      tmp2[xx,yy] = mean(data2[xx,yy,*])

      endfor
    endfor

    ave_data1 = tmp1
    ave_data2 = tmp2

window, 2, xs=1200, ys=600
!p.multi=[0,2,1]

  contour, ave_data1;, xst=1, yst=1
  tvscl, rebin(ave_data1, 3*sz[1], 3*sz[2]), 400, 400
  contour, ave_data2;, xst=1, yst=1
  tvscl, rebin(ave_data1, 3*sz[1], 3*sz[2]), 1000, 400

!p.multi=[0,1,0]


;data = midiChopImage(file)
chop_nod_phot, file






stop
end
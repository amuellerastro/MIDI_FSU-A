; EXECUTE IN MIA!

pro check_MIDI_data_qual

;spawn, '$mia'

; obs = 'P90RATZKA2'
; path = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obs+'/MIDIreduced_EWS_CustomMask/'
; ; star = 'CQTau_03:34:26'
; star = 'CQTau_03:30:16'

obs = 'P92CHOQUET'
path = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_'+obs+'/MIDIreduced_EWS_CustomMask/'
star = 'HD76111_03:34:51'

; star = '24Psc_01:18:17'
; star = 'NGC1068_02:44:20'
; star = 'HD16212_05:18:34'
; star = '24Psc_02:24:10'
; star = 'HD53179_07:11:10'
;star = 'RYTau_04:03:29'
;star = 'RYTau_04:55:48'
;star = 'betaPic_06:30:08'

;COMMP18
; path = '/media/disk/MIDI_FSUA/COMM18/MIDIreduced_EWS/'
;star = 'CD-482256_04:01:15'
; star = 'HD40270_03:05:32'
;star = 'HD32820_02:35:47'
; star = 'HD81797_07:38:04'
; star = 'HD23249_04:43:41'

;COMM16
;path = '/media/disk/MIDI_FSUA/COMM16/MIDIreduced_OPD_ZERO/MIDIreduced_EWS/'
;star = 'HD218594_09:20:46'
;star = 'HD218594_09:43:28'
;star = 'HD218594_09:52:41'
;star = 'HD1014_06:19:00'
;star = '24Psc_06:49:28'
;star = 'Gl86_07:20:16'
;star = 'HD11353_08:19:01'
;star = 'WDSJ05320-0018AaAb_10:07:26'

; path = '/media/disk/MIDI_FSUA/COMMtest/MIDIreduced_EWS/'
; star = 'HD1014_06:19:00'
;star = 'WDSJ05320-0018AaAb_10:07:26'


; path = '/media/disk/MIDI_data/HD144432/MIDIarchive/2006-05-15/results_2006-05-15/'
; star = 'HD144432_1'
; path = '/media/disk/MIDI_data/HD144432/MIDIarchive/2006-07-10/'
; star = 'HD135344.MIDI.2006-07-11T00:54:37'
; path = '/media/disk/MIDI_FSUA/COMM16/MIDIdata/HD25025_MIDIonly/'
; star = 'HD25025_10:10:08'



file1 = path+star+'.groupdelay.fits'
file2 = path+star+'.corr.fits'

; raw correlated amplitude and phase
v = oirgetvis(file2,wave=wave)

window, 3,xs=600,ys=450
!p.multi=[0,1,2]
plot, wave,v.visamp, xr=[6,14], charsize=1.2
oplot,wave,v.visamp*cos(v.visphi/!radeg),col=255
plot,wave,v.visphi, xr=[6,14], charsize=1.2, yr=[-10.,10.], yst=1
!p.multi=[0,1,0]

; GD
c=midigetcomplex(path+star,'groupdelay')
window, 4, xs=1000, ys=512
tv,rebin(abs(c)*0.01,1000,512)/9.d0

g=oirgetdata(file1)
gg=pseudocomplex(g.data1)
window, 0,xs=1024,ys=170
tvsclm,abs(gg[1000+indgen(1024),160:330])


gs=csmooth2(gg,10.)
window, 1,xs=1524,ys=512
tvsclm,abs(gs[1000+indgen(1524),*])
;tvsclm,rebin(abs(gs[1000+indgen(1024),256-127:256+128]),1024,512)


; ;phase delay???
; window, 2,xs=700,ys=512
; pd = float(gs[1000+indgen(1024),*])
; tvsclm, pd

window, 6, xs=400, ys=400
gm=total(gg[1000:*,*],1)
plot,abs(gm), charsize=1.5

stop
end
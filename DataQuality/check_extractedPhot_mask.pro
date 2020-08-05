;midiphotopipe, 'alfBoo_20130412', f, mask='median.srcmask.fits', skymask='median.skymask.fits'

pro check_extractedPhot_mask

;path of A and B Photometry

;================================================================================
;================================================================================

path = 'reduced_phot/'
srcmask = 'ftk/median.srcmask.fits'
skymask = 'ftk/median.skymask.fits'

;================================================================================
;================================================================================


file = file_search(path+'*.Aphotometry.fits', count=nfiles)

tmpID = file
for i=0,nfiles-1 do begin

  pos1 = strpos(file[i], '/')
  tmpID[i] = strmid(file[i], pos1+1, strlen(file[i])-pos1-18)

endfor

print, ''
for i=0,nfiles-1 do print, i+1, '  ', tmpID[i]
print, ''
read, 'Select photometry: ', quest

file = file[quest-1]
fileA = file

pos1 = strpos(file, 'Aphotometry')
strput, file, 'B', pos1
fileB = file

file = [fileA, fileB]


; file = ['HD181109_2013-05-22T06:45:57.Aphotometry.fits', 'HD181109_2013-05-22T06:45:57.Bphotometry.fits']

;================================================================================


;Masks

stm = mrdfits(srcmask,2, /silent)
; stm = mrdfits('HD102608-01.srcmask.fits',2, /silent)
; stm = mrdfits('minrtsMask_FIELD_PRISM_HIGH_SENS_SLIT.fits', 2, /silent)
mta = stm.data1
mtb = stm.data2

sts = mrdfits(skymask,2, /silent)
; sts = mrdfits('HD102608-01.skymask.fits',2, /silent)
msa = sts.data1
msb = sts.data2


;================================================================================


;Photometry

; file = 'MIDI.2013-04-30T09:44:12.000_01.fits'	;pA
; file = 'MIDI.2013-04-30T09:46:30.000_01.fits'	;pB

; file = 'MIDI.2013-04-30T10:00:03.000_01.fits'	;pA
; file = 'MIDI.2013-04-30T10:02:17.000_01.fits'	;pB

st0 = mrdfits(file[0],2, /silent)
data0A = st0[0].data1
data0B = st0[0].data2
st1 = mrdfits(file[1],2, /silent)
data1A = st1[0].data1
data1B = st1[0].data2

dum = mrdfits(file[0], 0, hdr, /silent)
cf = get_eso_keyword(hdr, 'HIERARCH ESO ISS CHOP FREQ')
object = get_eso_keyword(hdr, 'HIERARCH ESO OBS TARG NAME')

;================================================================================


x = dindgen(171)+1
y = dindgen(41)+1

window, 0, xs=1200, ys=900;, title=fileA
!p.multi=[0,2,2]

loadct, 3
  cgimage, data0A, /axes, xrange=[min(x), max(x)], yrange=[min(y), max(y)], $
      axkeywords={CHARSIZE:1.5, XTITLE:'Pixel', ytitle:'Pixel', xticklen:0.05}, position=[0.1,0.14,0.99,0.96], minvalue=min(data0A), maxvalue=max(data0A), background=fsc_color('black');, ctindex=32
  cgcontour, mta, x, y, /overplot, color='green', nlevels=3
  cgcontour, msa, x, y, /overplot, color='orange', nlevels=3
  plots, !x.crange, [7,7], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [11,11], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [23,23], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [27,27], color=fsc_color('gray'), linestyle=1

  legend, ['A OPEN / Beam A'], box=0, margin=0, charsize=1.5, /left
;   legend, [tmpID[quest-1], 'Chopping Frequency: '+sigfig(cf,2)], box=0, margin=0, charsize=2.0, /right
  legend, [tmpID[quest-1]], box=0, margin=1, charsize=2.0, /right

  cgimage, data0B, /axes, xrange=[min(x), max(x)], yrange=[min(y), max(y)], $
      axkeywords={CHARSIZE:1.5, XTITLE:'Pixel', ytitle:'Pixel', xticklen:0.05}, minvalue=min(data0A), maxvalue=max(data0A), position=[0.1,0.14,0.99,0.96]
  cgcontour, mtb, x, y, /overplot, color='green', nlevels=3
  cgcontour, msb, x, y, /overplot, color='orange', nlevels=3
  plots, !x.crange, [7,7], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [11,11], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [23,23], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [27,27], color=fsc_color('gray'), linestyle=1

  legend, ['A OPEN / Beam B'], box=0, margin=0, charsize=1.5, /left

;----------

  cgimage, data1A, /axes, xrange=[min(x), max(x)], yrange=[min(y), max(y)], $
      axkeywords={CHARSIZE:1.5, XTITLE:'Pixel', ytitle:'Pixel', xticklen:0.05}, position=[0.1,0.14,0.99,0.96], minvalue=min(data0A), maxvalue=max(data0A)

  cgcontour, mta, x, y, /overplot, color='green', nlevels=3
  cgcontour, msa, x, y, /overplot, color='orange', nlevels=3
  plots, !x.crange, [7,7], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [11,11], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [23,23], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [27,27], color=fsc_color('gray'), linestyle=1

  legend, ['B OPEN / Beam A'], box=0, margin=0, charsize=1.5, /left


  cgimage, data1B, /axes, xrange=[min(x), max(x)], yrange=[min(y), max(y)], $
      axkeywords={CHARSIZE:1.5, XTITLE:'Pixel', ytitle:'Pixel', xticklen:0.05}, minvalue=min(data0A), maxvalue=max(data0A), position=[0.1,0.14,0.99,0.96]
  cgcontour, mtb, x, y, /overplot, color='green', nlevels=3
  cgcontour, msb, x, y, /overplot, color='orange', nlevels=3
  plots, !x.crange, [7,7], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [11,11], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [23,23], color=fsc_color('gray'), linestyle=1
  plots, !x.crange, [27,27], color=fsc_color('gray'), linestyle=1

  legend, ['B OPEN / Beam B'], box=0, margin=0, charsize=1.5, /left


!p.multi=[0,1,0]



stop
end
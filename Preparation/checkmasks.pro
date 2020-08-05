pro checkmasks

path = 'MIDI_CustomMasks/'

; srcmask = 'HD76111_03:34:51.srcmask.fits'
; skymask = 'HD76111_03:34:51.skymask.fits'

srcmask = file_search(path+'*srcmask.fits', count=n)
skymask = file_search(path+'*skymask.fits', count=n)

data1src = dblarr(n,171,41)
data2src = dblarr(n,171,41)
data1sky = dblarr(n,171,41)
data2sky = dblarr(n,171,41)

for xx=0,n-1 do begin

  st = mrdfits(srcmask[xx],2,/silent)
  data1src[xx,*,*] = st.data1
  data2src[xx,*,*] = st.data2

  st = mrdfits(skymask[xx],2,/silent)
  data1sky[xx,*,*] = st.data1
  data2sky[xx,*,*] = st.data2

endfor

;median of all masks
mdata1src = dblarr(171,41)
mdata2src = dblarr(171,41)
mdata1sky = dblarr(171,41)
mdata2sky = dblarr(171,41)

;mean of all masks
mmdata1src = dblarr(171,41)
mmdata2src = dblarr(171,41)
mmdata1sky = dblarr(171,41)
mmdata2sky = dblarr(171,41)

for i=0,170 do begin

  for j=0,40 do begin

    mdata1src[i,j] = median(data1src[*,i,j], /even)
    mdata2src[i,j] = median(data2src[*,i,j], /even)
    mdata1sky[i,j] = median(data1sky[*,i,j], /even)
    mdata2sky[i,j] = median(data2sky[*,i,j], /even)

    mmdata1src[i,j] = mean(data1src[*,i,j])
    mmdata2src[i,j] = mean(data2src[*,i,j])
    mmdata1sky[i,j] = mean(data1sky[*,i,j])
    mmdata2sky[i,j] = mean(data2sky[*,i,j])

  endfor

endfor

;-------------------------------------------------------------------

;write out median averaged mask

;copy a sky and src mask as reference
refsrc = srcmask[0]	;'HD76111_03:37:09.srcmask.fits'
refsky = skymask[0]	;'HD76111_03:37:09.skymask.fits'
spawn, 'cp '+refsrc+' .'
spawn, 'cp '+refsky+' .'

;src mask
for i=0,2 do begin

  st = mrdfits(refsrc, i, hdr, /silent)

  if (i ne 2) then begin

    mwrfits, st, 'median.srcmask.fits', hdr, /silent

  endif else begin

    st.data1 = mdata1src
    st.data2 = mdata2src

    mwrfits, st, 'median.srcmask.fits', hdr, /silent

  endelse

endfor

;sky mask
for i=0,2 do begin

  st = mrdfits(refsky, i, hdr, /silent)

  if (i ne 2) then begin

    mwrfits, st, 'median.skymask.fits', hdr, /silent

  endif else begin

    st.data1 = mdata1sky
    st.data2 = mdata2sky

    mwrfits, st, 'median.skymask.fits', hdr, /silent

  endelse

endfor

; spawn, 'rm '+refsrc
; spawn, 'rm '+refsky
spawn, 'mv median.skymask.fits '+path+'.'
spawn, 'mv median.srcmask.fits '+path+'.'

;-------------------------------------------------------------------

;plot stuff

msrc1 = dblarr(n,41) & msrc2 = msrc1
msky1 = msrc1 & msky2 = msrc1

;for median
medsrc1 = dblarr(41)
medsrc2 = dblarr(41)
medsky1 = dblarr(41)
medsky2 = dblarr(41)

;for mean
mmedsrc1 = dblarr(41)
mmedsrc2 = dblarr(41)
mmedsky1 = dblarr(41)
mmedsky2 = dblarr(41)

for i=0,40 do begin

  medsrc1[i] = total(mdata1src[*,i])
  medsrc2[i] = total(mdata2src[*,i])
  medsky1[i] = total(mdata1sky[*,i])
  medsky2[i] = total(mdata2sky[*,i])

  mmedsrc1[i] = total(mmdata1src[*,i])
  mmedsrc2[i] = total(mmdata2src[*,i])
  mmedsky1[i] = total(mmdata1sky[*,i])
  mmedsky2[i] = total(mmdata2sky[*,i])

endfor

for xx=0,n-1 do begin

  st = mrdfits(srcmask[xx],2,/silent)
  data1 = st.data1
  data2 = st.data2

  st = mrdfits(skymask[xx],2,/silent)
  data3 = st.data1
  data4 = st.data2

  for i=0,40 do begin

    msrc1[xx,i] = total(data1[*,i])
    msrc2[xx,i] = total(data2[*,i])
    msky1[xx,i] = total(data3[*,i])
    msky2[xx,i] = total(data4[*,i])


  endfor


endfor

window, 0
!p.multi=[0,2,1]
plot, msrc1[0,*], /nodata
for i=0,n-1 do begin
  oplot, msrc1[i,*]
;   if (i eq 1) then oplot, msrc1[i,*], color=cgcolor('red')
endfor
  oplot, medsrc1, color=cgcolor('red'), thick=2
  oplot, mmedsrc1, color=cgcolor('green'), thick=2

plot, msrc2[0,*], /nodata
for i=0,n-1 do begin
  oplot, msrc2[i,*]
;   if (i eq 1) then oplot, msrc2[i,*], color=cgcolor('red')
endfor
  oplot, medsrc2, color=cgcolor('red'), thick=2
  oplot, mmedsrc2, color=cgcolor('green'), thick=2

window, 1
!p.multi=[0,2,1]
plot, msky1[0,*], /nodata
for i=0,n-1 do begin
  oplot, msky1[i,*]
;   if (i eq 1) then oplot, msky1[i,*], color=cgcolor('red')
endfor
  oplot, medsky1, color=cgcolor('red'), thick=2
  oplot, mmedsky1, color=cgcolor('green'), thick=2
  oplot, medsrc1, color=cgcolor('blue'), thick=2

plot, msky2[0,*], /nodata
for i=0,n-1 do begin
  oplot, msky2[i,*]
;   if (i eq 1) then oplot, msky2[i,*], color=cgcolor('red')
endfor
oplot, medsky2, color=cgcolor('red'), thick=2
  oplot, mmedsky2, color=cgcolor('green'), thick=2
  oplot, medsrc2, color=cgcolor('blue'), thick=2

!p.multi=[0,1,0]

window, 2
plot, medsrc1
oplot, medsky1


stop
end
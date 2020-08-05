function normalize, x

  return, (x-min(x))/(max(x)-min(x))

end


pro fringe

; ;lambdas and bandwidth
; wl = dblarr(5,4,6)
; for i=0,4 do begin
; 
;   readcol, 'l'+strcompress(i+1, /rem)+'.dat', dum, format='d'
;   wl[i,*,0] = dum[0:3]
;   wl[i,*,1] = dum[4:7]
;   wl[i,*,2] = dum[8:11]
;   wl[i,*,3] = dum[12:15]
;   wl[i,*,4] = dum[16:19]
;   wl[i,*,5] = dum[20:23]
; 
; endfor
; 
; lam = dblarr(6)
; slam = lam
; 
; lam[0] = median(wl[*,*,0]) & slam[0] = stddev(wl[*,*,0])
; lam[1] = median(wl[*,*,1]) & slam[1] = stddev(wl[*,*,1])
; lam[2] = median(wl[*,*,2]) & slam[2] = stddev(wl[*,*,2])
; lam[3] = median(wl[*,*,3]) & slam[3] = stddev(wl[*,*,3])
; lam[4] = median(wl[*,*,4]) & slam[4] = stddev(wl[*,*,4])
; lam[5] = median(wl[*,*,5]) & slam[5] = stddev(wl[*,*,5])
; 
; diff = ts_diff(lam, 1)
; diff = diff[1:4]
; dl = abs(median(diff))
; lcoh = (lam^2.-(dl^2.)/4.d0)/dl
; print, lcoh
; 
; delta = (dindgen(2001.)-1000.d0)/10.d6
; eta = 1.d0
; I0 = 1.d0
; 
; k = 2.d0*!DPI/lam
; sinc = dblarr(6, n_elements(delta))
; Iint = dblarr(6, n_elements(delta))
; 
; for i=0,5 do begin
; 
;   sinc[i,*] = (sin(!DPI*delta/lcoh[i]))/(!DPI*delta/lcoh[i])
;   Iint[i,*] = (2.d0*I0*dl*eta*sinc[i,*]*sin(k[i]*delta))
; 
;   idx = where(finite(Iint[i,*]) ne 1)
;   Iint[i,idx] = 0.
;   Iint[i,*] = Iint[i,*]/min(Iint[i,*])
; 
; endfor
; 
; ; window, 1
; ; plot, delta*1d6, Iint, charsize=2, xst=1, xr=[-20,20]
; 
; ;-------
; ;lambdas from sky calib for comparison
; lam2 = lam
; 
; lam2[0] = median([2.2447493e-06,2.2380676e-06,2.2604940e-06,2.2640065e-06])
; lam2[1] = median([2.0274317e-06,2.0449362e-06,2.0245543e-06,2.0396177e-06])
; lam2[2] = median([2.1533725e-06,2.1611392e-06,2.1561737e-06,2.1509340e-06])
; lam2[3] = median([2.2779620e-06,2.2858453e-06,2.2835266e-06,2.2815958e-06])
; lam2[4] = median([2.4154613e-06,2.4216423e-06,2.4280512e-06,2.4199824e-06])
; lam2[5] = median([2.5148293e-06,2.5179607e-06,2.5396787e-06,2.5298647e-06])
; 
; diff2 = ts_diff(lam2, 1)
; diff2 = diff2[1:4]
; dl2 = abs(median(diff2))
; lcoh2 = (lam2^2.-(dl2^2.)/4.d0)/dl2
; print, lcoh2
; 
; 
; delta = (dindgen(2001.)-1000.d0)/10.d6
; eta = 1.d0
; I0 = 1.d0
; 
; k2 = 2.d0*!DPI/lam2
; sinc2 = dblarr(6, n_elements(delta))
; Iint2 = dblarr(6, n_elements(delta))
; 
; for i=0,5 do begin
; 
;   sinc2[i,*] = (sin(!DPI*delta/lcoh2[i]))/(!DPI*delta/lcoh2[i])
;   Iint2[i,*] = (2.d0*I0*dl2*eta*sinc2[i,*]*sin(k2[i]*delta))
; 
;   idx = where(finite(Iint2[i,*]) ne 1)
;   Iint2[i,idx] = 0.
;   Iint2[i,*] = Iint2[i,*]/min(Iint2[i,*])
; 
; endfor
; 


;=================

;see PhD Ratzka p126, Tristram p17+

;parameters
;N band
; l0 = 10.5d-6
; dl0 = 5.d-6
; r = 30.d0
;K band
l0 = 2.4d-6
dl0 = 0.145d-6
r = 5.d0
doffset = -30.d-6
Ioffset = 0.5d0

;coherence length
lcoh1 = (l0^2.-(dl0^2.)/4.d0)/dl0
lcoh2 = (l0^2.)/dl0
lcoh3 = r*l0-l0/(4.d0*r)
print, 'in microns', lcoh1*1.d6, lcoh2*1.d6, lcoh3*1.d6


;fringe pattern
lcoh = lcoh1
delta = (dindgen(2001.)-1000.d0)/10.d6
;delta = delta[where(delta ne 0.d0)]
eta0 = 1.d0
I0 = -1000000.d0

k0 = 2.d0*!DPI/l0
sinc = (sin(!DPI*(delta-doffset)/lcoh))/(!DPI*(delta-doffset)/lcoh)
idx = where(finite(sinc) ne 1)
if (idx[0] ne -1) then sinc[idx] = 1.d0

gauss = 1.*exp(-0.5d0*((delta-doffset)/(lcoh*2.d0))^2.)

Iint = (2.d0*I0*dl0*eta0*sinc*sin(k0*(delta-doffset)))+Ioffset
Iintphase = (2.d0*I0*dl0*eta0*sinc*sin(k0*(delta-doffset)+!DPI/2.))+Ioffset

stop
;Iint = Iint/min(Iint)

delta = delta*1.d6
window, 0
plot, delta, Iint, charsize=2, xst=1, xr=[-70,00], yr=[0,1], yst=1
oplot, delta, Iintphase, color=rgb(255,0,0)

;plot, delta-3.9d-3, Iint, charsize=2, xst=1, xr=[-3.94d-3,-3.86d-3]



;for a fit variables are: dl, l0, I0, doffset, Ioffset


stop
end
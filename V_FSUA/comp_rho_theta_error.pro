@boot_mean.pro

pro comp_rho_theta_error

;HD155826
;values
; period = 14.215d0
; T = 1985.98d0
; ec = 0.4912d0
; a = 0.2527d0
; i = 115.2d0
; omega = 190.41d0
; omega2 = 135.2d0
; ; t2 = 2014.d0
; ;errors
; periode = 0.05d0
; Te = 0.17d0
; ece = 0.0048d0
; ae = 0.0043d0
; ie = 1.1d0
; omegae = 0.62d0
; omega2e = 2.5d0
; 
; bina = 0.1
; bind = 1.
; 
; tmp = (3.d0+30.d0/60.d0)/24.d0	;hh,mm
; t2 = 2013.d0+7.d0/12.d0+(7.d0+tmp)/365.25d0
; CDF_EPOCH, MergeDate, 2013.d0, 7.d0, 7.d0, 3.d0, 30.d0, /COMPUTE_EPOCH 

;===================================================

;24Psc
; 
period = 22.81d0
T = 1988.72d0
ec = 0.422d0
a = 0.0832d0
i = 133.7d0
omega = 209.5d0	;node
omega2 = 298.3d0	;argument of periastron

periode = 0.15d0
Te = 0.14d0
ece = 0.013d0
ae = 0.0014d0
ie = 1.8d0
omegae = 2.7d0
omega2e = 2.4d0

bina = 0.1
bind = 0.1

tmp = (4.d0+30.d0/60.d0)/24.d0	;hh,mm
t2 = 2013.d0+10.d0/12.d0+(29.d0+tmp)/365.25d0

;===================================================
;===================================================

;reference
result = rhotheta(period, T, ec, a, i, omega, omega2, t2)
refrho = result.rho
reftheta = result.theta

;create N values based on error bars
;===================================================
n = 1.d4
;===================================================
nperiod = (period-periode)+(2.*periode)*(randomu(seed,n))
nT = (T-Te)+(2.*Te)*(randomu(seed,n))
nec = (ec-ece)+(2.*ece)*(randomu(seed,n))
na = (a-ae)+(2.*ae)*(randomu(seed,n))
ni = (i-ie)+(2.*ie)*(randomu(seed,n))
nomega = (omega-omegae)+(2.*omegae)*(randomu(seed,n))
nomega2 = (omega2-omega2e)+(2.*omega2e)*(randomu(seed,n))

rho = dblarr(n)
theta = rho

for i=0L,n-1 do begin

  tmp = rhotheta(nperiod[i], nT[i], nec[i], na[i], ni[i], nomega[i], nomega2[i], t2)
  rho[i] = tmp.rho	;arcsec
  theta[i] = tmp.theta	;deg

endfor

alpha = rho*sin(theta*!dtor)*1.d3	;mas
delta = rho*cos(theta*!dtor)*1.d3	;mas

;non-gaussian distribution
;===================================================
nboot = 100000.d0
;===================================================
boot_mean, alpha, nboot, tmp
meanalpha = mean(tmp)
boot_mean, delta, nboot, tmp
meandelta = mean(tmp)

plothist, alpha, xhista, yhista, bin=0.1, /noplot
plothist, delta, xhistd, yhistd, bin=1., /noplot

areaa = tsum(xhista, yhista)
aread = tsum(xhistd, yhistd)

i=1
repeat begin

  tmp = tsum(xhista[0:i], yhista[0:i])
  i = i+1

endrep until (tmp ge 0.16*areaa)
minusalpha = meanalpha-xhista[i]

i=1
repeat begin

  tmp = tsum(xhista[0:i], yhista[0:i])
  i = i+1

endrep until (tmp ge 0.84*areaa)
plusalpha = xhista[i]-meanalpha

i=1
repeat begin

  tmp = tsum(xhistd[0:i], yhistd[0:i])
  i = i+1

endrep until (tmp ge 0.16*aread)
minusdelta = meandelta-xhistd[i]

i=1
repeat begin

  tmp = tsum(xhistd[0:i], yhistd[0:i])
  i = i+1

endrep until (tmp ge 0.84*aread)
plusdelta = xhistd[i]-meandelta

window, 0, xs=1000, ys=1000
!p.multi=[0,1,2]

plothist, alpha, xhista, yhista, bin=0.1
plots, (meanalpha-minusalpha)*[1,1], !y.crange, color=fsc_color('green'), /data
plots, meanalpha*[1,1], !y.crange, color=fsc_color('red'), /data
plots, (plusalpha+meanalpha)*[1,1], !y.crange, color=fsc_color('green'), /data

plothist, delta, xhistd, yhistd, bin=1.
plots, (meandelta-minusdelta)*[1,1], !y.crange, color=fsc_color('green'), /data
plots, meandelta*[1,1], !y.crange, color=fsc_color('red'), /data
plots, (plusdelta+meandelta)*[1,1], !y.crange, color=fsc_color('green'), /data


!p.multi=[0,1,0]

print, ''
print, 'predicted values'
print, 'alpha minusE plusE'
print, meanalpha, minusalpha, plusalpha
print, 'delta minusE plusE'
print, meandelta, minusdelta, plusdelta


stop
end
function fit_ccf_gauss_3terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)
  return, fit
end
function fit_ccf_gauss_4terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]
  return, fit
end
function fit_ccf_gauss_5terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)
  return, fit
end
function fit_ccf_gauss_6terms, x, p
  fit = p[0]*exp(-0.5d0*((x-p[1])/p[2])^2.)+p[3]+(p[4]*x)+(p[5]*x^2.)
  return, fit
end

function shah, x, p
  fit = p[0] + p[1]*x + p[2]*p[3]^x
  return, fit
end

function threeparamexp, x, p
  fit = exp(p[0] + p[1]*x + p[2]*x^2.)
  return, fit
end

function explinear, x, p
  fit = p[0]*exp(-1.d0*x/p[1])+p[2]+p[3]*x
  return, fit
end

function parabel, x, p
  fit = p[0]*x^2. + p[1]*x + p[2]
  return, fit
end

function binsizeFSU, x

;http://176.32.89.45/~hideaki/res/histogram.html
;Shimazaki and Shinomoto, Neural Comput 19 1503-1527, 2007

; I. Divide the data range into $N$ bins of width $\Delta $. Count the number of events $k_i$ that enter the i'th bin.
; II. Calculate the mean and variance of the number of events as <formula>
; III. Compute a <formula>,
; IV. Repeat i-iii while changing $\Delta $. Find $\Delta^{*}$ that minimizes $C_n (\Delta )$.
; *VERY IMPORTANT: Do NOT use a variance that uses N-1 to divide the sum of squared errors. Use the biased variance in the method.
; *To obtain a smooth cost function, it is recommended to use an average of cost functions calculated from multiple initial partitioning positions. 

; bins = dindgen(2500.)+2.d0	;max. number of bins to be considered, bins[0]=2 but there is only 1 bin
bins = dindgen(1000.)+2.d0	;max. number of bins to be considered, bins[0]=2 but there is only 1 bin

binwidth = dblarr(n_elements(bins))
nbins = dblarr(n_elements(bins))
k = dblarr(n_elements(bins)) & v = dblarr(n_elements(bins))
c = dblarr(n_elements(bins))

for i=0L,n_elements(bins)-1 do begin

  binarr = linspace(min(x), max(x), bins[i], step=stepsize)	;dividing data into n-bins, i.e. i+1 bins

  binwidth[i] = stepsize
  nbins[i] = i+1

  events = dblarr(nbins[i])
  aveevents = dblarr(nbins[i])
  stddevevents = dblarr(nbins[i])
  for j=0,nbins[i]-1 do begin

    idx = where(x ge binarr[j] and x le binarr[j+1])
    if (idx[0] ne -1) then events[j] = n_elements(idx) else events[j] = 0.d0

  endfor

  k[i] = (total(events))/nbins[i]
  v[i] = total((events-k[i])^2.d0)/nbins[i]
  c[i] = (2.d0*k[i]-v[i])/stepsize^2.d0

;   proceeding_text, loop=n_elements(bins), i=i, prompt='> Processing                     '+string(i+1,form='(I4)')

endfor


; window, 0, xs=1500, ys=500
; plot, binwidth, c, /xlog, xtitle='Bin Size', ytitle='C!Dn!N', xst=1, charsize=1.5;, xr=[1, max(binwidth)]
; 
c = smooth(c,11)
; oplot, binwidth, c, color=fsc_color('green')
; 
; tmp = min(c, idxmin)
; binopt = binwidth[idxmin]
; plots, [binopt, binopt], !y.crange, linestyle=1

;fit a polynomial nth degree to c
params = poly_fit(alog10(binwidth), c, 6, yfit=yfit)
; oplot, binwidth, yfit, color=fsc_color('red')
tmp = min(yfit, idxmin)
binopt = binwidth[idxmin]
; plots, [binopt, binopt], !y.crange, linestyle=1, color=fsc_color('red')






;uncomment if you want to fit the data

; ; ; 
; ; ; print, 'select region to fit'
; ; ; cursor, x1, y1, /data, /down
; ; ; cursor, x2, y2, /data, /down
; ; ; ; print, 'select minimum'
; ; ; ; cursor, x3, y3, /data, /down
; ; ; 
; ; ; idx = where(binwidth ge x1 and binwidth le x2)
; ; ; binwidth_cut = binwidth[idx]
; ; ; c_cut = c[idx]
; ; ; 
; ; ; ; dummy_err = dblarr(n_elements(c_cut))
; ; ; ; dummy_err[*] = 1.d0
; ; ; ; start_val = [y3, 1., 1., 1.]
; ; ; ; fitparams = mpfitfun('shah', binwidth_cut, c_cut, dummy_err, start_val, weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
; ; ; ; oplot, binwidth_cut, yfit, color=fsc_color('yellow')
; ; ; ; print, 'optimum bin size: ', alog(-fitparams[1]/(fitparams[2]*alog(fitparams[3])))/alog(fitparams[3])
; ; ; 
; ; ; ; dummy_err = dblarr(n_elements(c_cut))
; ; ; ; dummy_err[*] = 1.d0
; ; ; ; start_val = [1., 1., 0.5]
; ; ; ; fitparams = mpfitfun('threeparamexp', binwidth_cut, c_cut, dummy_err, start_val, weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof)
; ; ; ; oplot, binwidth_cut, yfit, color=fsc_color('red')
; ; ; 
; ; ; ;GAUSS
; ; ; gauss_nterms = '4'
; ; ; dummy_err = dblarr(n_elements(c_cut))
; ; ; dummy_err[*] = 1.d0
; ; ; start_val = [-0.02d0, 30.d0, 100.d0, 0.d0]
; ; ; fitparams = mpfitfun('fit_ccf_gauss_'+gauss_nterms+'terms', alog10(binwidth_cut), c_cut, dummy_err, start_val, weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
; ; ; oplot, binwidth_cut, yfit, color=fsc_color('red')
; ; ; 
; ; ; binopt = 10.d0^fitparams[1]
; ; ; print, 'optimum bin size: ', binopt
; ; ; ;print, fitparams[1]*fitparams[3]+fitparams[2]
; ; ; plots, [binopt, binopt], !y.crange, linestyle=1, color=fsc_color('red')
; ; ; 
; ; ; stop
; ; ; 
; ; ; ;PARABEL
; ; ; dummy_err = dblarr(n_elements(c_cut))
; ; ; dummy_err[*] = 1.d0
; ; ; start_val = [10.d0, 6.d0, -0.018d0]
; ; ; fitparams = mpfitfun('parabel', alog10(binwidth_cut), c_cut, dummy_err, start_val, weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
; ; ; oplot, binwidth_cut, yfit, color=fsc_color('red')
; ; ; 
; ; ; binopt = 10.d0^(-1.d0*fitparams[1]/(2.d0*fitparams[0]))
; ; ; print, 'optimum bin size: ', binopt
; ; ; ;print, fitparams[1]*fitparams[3]+fitparams[2]
; ; ; plots, [binopt, binopt], !y.crange, linestyle=1, color=fsc_color('red')
; ; ; 
; ; ; ;EXPLINEAR
; ; ; dummy_err = dblarr(n_elements(c_cut))
; ; ; dummy_err[*] = 1.d0
; ; ; start_val = [0.02d0, 6.d0, -0.02d0, 0.d0]
; ; ; fitparams = mpfitfun('explinear', binwidth_cut, c_cut, dummy_err, start_val, weights=dummy_err, maxiter=2000, niter=niter, status=status, bestnorm=bestnorm, yfit=yfit, perror=perror, dof=dof, /quiet)
; ; ; oplot, binwidth_cut, yfit, color=fsc_color('yellow')
; ; ; 
; ; ; binopt = -1.d0*fitparams[1]*alog(fitparams[1]*fitparams[3]/fitparams[0])
; ; ; print, 'optimum bin size: ', binopt
; ; ; ;print, fitparams[1]*fitparams[3]+fitparams[2]
; ; ; plots, [binopt, binopt], !y.crange, linestyle=1, color=fsc_color('yellow')
; ; ; 
; ; ; stop

return, binopt

end
@showsym.pro

pro plot_Vis_P94KOEHLER

file1 = 'summary_FSUAscans_P94KOEHLERD0I1.txt'
readcol, file1, jd1, bl1, pa1, v1, ve1, format='d,d,d,d,d'
file2 = 'summary_FSUAscans_P94KOEHLERH0I1.txt'
readcol, file2, jd2, bl2, pa2, v2, ve2, format='d,d,d,d,d'


; jd = jd-2450000.5d0
; jd = (jd-fix(jd[0]))*24.

window, 0, xs=500, ys=500

; ploterror, bl, v, ve, psym=sym(1), yr=[0,0.5], symsize=2, $
;   xtitle='test', ytitle='V^2'

ploterror, [bl1,bl2], [v1,v2], [ve1,ve2], psym=sym(1), yr=[0,0.5], symsize=2, $
  xtitle='Projected Baseline [m]', ytitle='V^2', charsize=1.5



stop
end

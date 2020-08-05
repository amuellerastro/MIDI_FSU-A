;*******************************************************************************
;+
;*NAME:
;
;      CRSPROD (IUE Production Library) (01 April 1988)
;
;*CLASS:
;
;     cross-correlation
;
;*CATEGORY:
;
;*PURPOSE:
;
;   To compute a normalized cross correlation for two spectral segments
;   which are sampled on the same linear or log(lambda) scale. 
;
;*CALLING SEQUENCE:
;
;   CRSPROD,FFIR,FSEC,NSPR,CROSS,CRMIN,CRMAX
;
;*PARAMETERS:
;
;   FFIR    (REQ) (I)  (1) (F)
;           Required input vector giving the flux data for the first
;           spectrum.
;
;   FSEC    (REQ) (I)  (1) (F)
;           Required input vector giving the flux data for the second 
;           spectrum.
;
;   NSPR    (REQ) (I)  (0) (F)
;           Required input parameter specifying the spectral range to
;           be considered in the cross-correlation function.
;
;   CROSS   (REQ) (O)  (1) (F)
;           Required output vector containing the cross-correlation 
;           function.
;
;   CRMIN   (REQ) (O)  (0) (F)
;           Required output vector containing the minimum of the 
;           cross-correlation function. 
;
;   CRMAX   (REQ) (O)
;           Required output vector containing the maximum of the 
;           cross-correlation function. 
;
;*EXAMPLES:
;
;    To compute the cross-correlation function for two spectra, FIRST
;    and SECOND, using the recommended initial spectral range from CRSCOR,
;
;    CRSPROD,FIRST,SECOND,15,CROSS,CRMIN,CRMAX
;
;*SYSTEM VARIABLES USED:
;
;*INTERACTIVE INPUT:
;
;*SUBROUTINES CALLED:
;
;    PARCHECK
;
;*FILES USED:
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;*NOTES:
;       Assumes same number of elements in both spectra. (Both fluxes are
;       divided by the number of elements in the first spectrum.)
;
;*PROCEDURE:
;
;     CROSS is determined for (2*nspr + 1) tags or shifts going from -15
;     to +15 shifts from the starting locations. 
;     After subtracting the average flux from each spectrum, the cross
;     correlation function is computed as follows for each point in 
;     the spectra, 
;      TEMP = (second spectrum) * SHIFT(first spectrum,ns)
;      CROSS(L) = TOTAL(TEMP(ls:us)/nele) 
;
;
;*MODIFICATION HISTORY:
;
;	25 Jun 1991  PJL cleaned up; added PARCHECK and parameter eq 0
;			 print; tested on SUN and VAX; updated prolog
;
;-
;*******************************************************************************
 pro crsprod_MOD,ffir,fsec,nspr,cross,crmin,crmax
;
;  npar = n_params(0)
;  if npar eq 0 then begin
;     print,' CRSPROD,FFIR,FSEC,NSPR,CROSS,CRMIN,CRMAX'
;     retall
;  endif  ; npar
;  parcheck,npar,6,'CRSPROD'
;
;  compute a normalized cross correlation function
;
 mi   = n_elements(ffir) 
  avg1 = total(ffir)/mi & ff= ffir - avg1    ; subtract average fluxes 
  avg2 = total(fsec)/mi & fs= fsec - avg2
; ff = ffir
; fs = fsec
 ntot = nspr+nspr+1
 cross= dblarr(ntot)
 temp = fs
 for l=0,ntot-1 do begin
    ns = nspr - l
    temp = fs*shift(ff,ns)
    ls = ns > 0
    us = mi  - 1 + (ns < 0)
    nele = us - ls + 1
    cross(l) = total(temp(ls:us)) / nele
 endfor  ; l

crmax= max(cross) & crmin= min(cross)
cross= (cross -crmin)/(crmax-crmin)       ; normalize function

 return     
 end  ; crsprod


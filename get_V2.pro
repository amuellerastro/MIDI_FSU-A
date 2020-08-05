;***************
;ONLY FSU A!!!!*
;***************

; compte V^2(t) from FSU data

pro get_V2, file, V2

;file='PACMAN_OBS_GENERIC207_0013.fits'
;target = 'test'

;read main header
dummy = readfits(file,head)

darkd = get_eso_keyword(head,'HIERARCH ESO OCS FSUA KWDARK')
flatd = get_eso_keyword(head,'HIERARCH ESO OCS FSUA KWFLAT')
phased = get_eso_keyword(head,'HIERARCH ESO OCS FSUA KWPHAS')


phase = dblarr(4) & dark = phase & flat = phase
for i=0,3 do begin

  pos = strpos(phased,',')
  phase(i) = strmid(phased,0,pos)
  phased = strmid(phased,pos+1,strlen(phased))
  if (pos eq -1) then phase(i) = phased


  pos = strpos(darkd,',')
  dark(i) = strmid(darkd,0,pos)
  darkd = strmid(darkd,pos+1,strlen(darkd))
  if (pos eq -1) then dark(i) = darkd


  pos = strpos(flatd,',')
  flat(i) = strmid(flatd,0,pos)
  flatd = strmid(flatd,pos+1,strlen(flatd))
  if (pos eq -1) then flat(i) = flatd

endfor

alpha = phase(0)
beta = phase(1)
gamma = phase(2)
delta = phase(3)


csc = 2./(beta*gamma - alpha*delta)

;read table (imaging_data_fsuX)
;if (fsu eq 'FSUA') then begin

st = mrdfits(file,9)

flux1 = (st.data1[0] - dark[0]) / (flat[0] - 2.*dark[0])
flux2 = (st.data2[0] - dark[1]) / (flat[1] - 2.*dark[1])
flux3 = (st.data3[0] - dark[2]) / (flat[2] - 2.*dark[2])
flux4 = (st.data4[0] - dark[3]) / (flat[3] - 2.*dark[3])


;4 data arrays (A,B,C,D channel) each with 6 columns where 1st columns is the UNNORMALIZED white light pixel and the last are the dispersed pixel

;X=[(A-C)*gamma - (B-D)*alpha]*csc
;Y=[(B-D)*beta - (A-C)*delta]*csc
X = ((flux1-flux3)*gamma - (flux2-flux4)*alpha)*csc
Y = ((flux2-flux4)*beta - (flux1-flux3)*delta)*csc
;stop
V2 = (0.5*sqrt(X*X + Y*Y))^2


end
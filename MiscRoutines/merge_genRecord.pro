;in case the fits files are not merged
;merged tha they meet the structure of a merged file even some extensions might be missing


pro merge_genRecord

;============================================================
;directory with files to be merged
path = '/media/disk_MIDIFSU/MIDI_FSUA/MIDI_FSUA_P90POTT/FSUAdata/'

;file ID to be merged
id = '310_0007'

;reference file for missing extensions
ref_id = '310_0006'


filebase = 'PACMAN_GEN_RECORD_'
hdr = ['_hdr.fits', '_rmn_hdr.fits']

;extension to be merged
ext = ['', '', '', '', '', '', '_DOPDC.fits', '_IMAGING_DATA_FSUA.fits', '_IMAGING_DATA_FSUB.fits', '_METROLOGY_DATA.fits', '_METROLOGY_DATA_FSUB.fits', '_OPDC.fits']

extnum = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]	;extension number, missing extensions will be filled with dummy

;============================================================
;next = n_elements(ext)

result = filebase+id+'.fits'


;new header, extension 0
struc = mrdfits(path+filebase+ref_id+'.fits', 0, oldhdr, /silent)
readcol, path+filebase+id+hdr[0], t1, t2, format='a,a', delimiter='='
readcol, path+filebase+id+hdr[1], t3, t4, format='a,a', delimiter='='
; h1 = t1+'='+t2
; h2 = t3+'='+t4
; prime = "'"
; tmp = 'SIMPLE   = '+prime+'T'+prime+'    /Standard FITS (NOST-100-2.0)'
; newhdr = [tmp, h1, h2]
; modfits, path+result, 0, newhdr

keyword = [t1,t3]
tmp = [t2,t4]

pos = strpos(tmp, '/')
value = strarr(n_elements(keyword))
comment = strarr(n_elements(keyword))

for i=0,n_elements(keyword)-1 do begin

  value[i] = strmid(tmp[i], 0, pos[i])
  comment[i] = strmid(tmp[i], pos[i]+1, strlen(tmp[i])-pos[i])
  sxaddpar, oldhdr, keyword[i], value[i], comment[i]

endfor

mwrfits, struc, path+result, oldhdr, /silent


for i=1,11 do begin

;   if (i eq 0) then begin
; 
;     struc = mrdfits(path+filebase+ref_id+'.fits', i, refhead, /silent)
;     mwrfits, struc, path+result, refhead, /silent
; 
;   endif


  if (i lt 6 and i ge 1) then begin

    struc = mrdfits(path+filebase+ref_id+'.fits', i, head, /silent)
    mwrfits, struc, path+result, head, /silent

  endif

  if (i ge 6) then begin

    struc = mrdfits(path+filebase+id+ext[i], 1, head, /silent)
    mwrfits, struc, path+result, head, /silent

  endif

endfor






stop
end
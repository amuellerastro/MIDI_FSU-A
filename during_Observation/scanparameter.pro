pro scanparameter

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; % small tool to compute the (D)OPDC scanning parameter for the simple scan
; % to find the fringe for the first time 
; % (one single scan and recording of file)
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; % csc 2010-09-15
; % -------------------------------------------------------------------------

scanrange = 0.04d0	;% range that should be scanned in meter
offset = 0.0d0	;% position around which should be scanned in meter
dlsign = -1	;% (D)DL sign
dit = 0.001d0	;% Detector Integration Time in sec
Frate = 1000.0d0;         % in Hz

;%% ========================================================================

; %speed = 0.01/50;        % a speed of 1 cm in 50 s has turned
;                         % out to be reasonable for FSU @ 1KHz 

nsample = 6.0d0	;6	; % number of samples per fringes
speed = 2.2d-6/nsample*Frate; % nsample per 2.2um at Frate 
; % =========================================================================
; % ------------------ Do NOT edit below!!! ---------------------------------
; % =========================================================================

ZPDoffset = ((offset - (scanrange/2.d0*dlsign)))*dlsign
scanperiod = 1.d0/speed * scanrange * 4.d0 * dit/0.001d0
filel = scanrange/speed * dit/0.001 + 10
print, ''
print, '-- in Gui -----------'
print, 'ZPD offset     :', ZPDoffset
print, 'Scan amplitude :', scanrange
print, 'Scan period    :', scanperiod
print, '-- in BOB -----------'
print, 'File length    :', filel
print, ''
;% -------------------------------------------------------------------------



stop
end
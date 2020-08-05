; ;wrap [rad] angle to [-PI..PI)
; function WrapPosNegPI, data
; 
;   return, (((data mod (ddp)) + ddp) mod ddp) - ddp/2.
; 
; end
; 
; 
; ;wrap [deg] angle to [-180..180)
; function WrapPosNeg180, data
; 
;   return, (((data mod 360.d0) + 360.d0) mod 360.d0) - 180.d0
; 
; end


;compared to the previous function the following functions have a 2*PI offset?...

;http://www.22answers.com/posts/answers/en/4633177

function specialmod, x, y

  return, x-y*floor(x/y)

end

;wrap [rad] angle to [-PI..PI)
function WrapPosNegPI, data

  return, specialmod(data+!DPI, 2.d0*!DPI)-!DPI

end


;wrap [rad] angle to [0..TWO_PI)
function WrapTwoPI, data

  return, specialmod(data, 2.d0*!DPI)

end

;wrap [deg] angle to [-180..180)
function WrapPosNeg180, data

  return, specialmod(data + 180.d0, 360.d0) - 180.d0

end

;wrap [deg] angle to [0..360)
function Wrap360, data

  return, specialmod(data ,360.d0)

end

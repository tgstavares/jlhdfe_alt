clear all
cd "/home/tgst/Desktop/Testes_stata"

// net install jlhdfe_alt, replace from("https://raw.githubusercontent.com/tgstavares/jlhdfe_alt/main/")
// net from "https://raw.githubusercontent.com/tgstavares/jlhdfe_alt/main/"
// net describe jlhdfe_alt
// net install jlhdfe_alt, replace

* Create some made-up data
local ebase=1e4
local nobs=10000*`ebase'
set obs `nobs'
gen id1 = runiformint(1, 1000*`ebase')
gen id2 = runiformint(1, 100*`ebase')
gen id3 = runiformint(1, 10*`ebase')
gen id4 = runiformint(1, 1*`ebase')
gen aux = cos(id1) + sin(id2) + sin(id3)^2 + cos(id4)^2

local nvars=20
local X
local erro_X
forvalues i=1/`nvars'{
	local X `X' x`i'
	local erro_X `erro_X' erro_x`i'
}
drawnorm `erro_X'
forvalues i=1/`nvars'{
	gen x`i' = aux*exp(erro_x`i')
}
gen y0 = 0.0
local XX
local nvars2=`nvars'-3
forvalues i=1/`nvars2'{
	replace y0 = y0 + `i'*x`i' + x`i'*x20 + x`i'*x19 + x`i'*x18
	local XX `XX' c.x`i'#c.x20 c.x`i'#c.x19 c.x`i'#c.x18
}
replace y0 = y0 + aux + rnormal() * `ebase'

compress _all

capture quietly log using "Performance_`ebase'.txt", text replace
capture set rmsg on
reg        y0 `X' `XX'
reghdfejl  y0 `X' `XX', absorb(id1 id2 id3 id4) vce(cluster id1) keepsingletons gpu
jlhdfe_alt y0 `X' `XX', absorb(id1 id2 id3 id4) vce(cluster id1) method(cuda) threads(48) outprefix(fem_out) display
capture set rmsg off
capture quietly log close

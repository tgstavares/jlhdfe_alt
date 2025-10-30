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
gen aux = cos(id1) + sin(id2) + sin(id3)^4

local nvars=10
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
forvalues i=1/`nvars'{
	replace y0 = y0 + `i'*x`i' + x`i'*x10 + x`i'*x9 + x`i'*x8
	local XX `XX' c.x`i'#c.x10 c.x`i'#c.x9 c.x`i'#c.x8
}
replace y0 = y0 + aux + rnormal() * `ebase'

compress _all

capture quietly log using "Performance_`ebase'.txt", text replace
capture set rmsg on
reg        y0 `X' `XX'
reghdfejl  y0 `X' `XX', absorb(id1 id2 id3) vce(cluster id1 id2 id3) keepsingletons gpu
jlhdfe_alt y0 `X' `XX', absorb(id1 id2 id3) vce(cluster id1 id2 id3) method(cuda) threads(48) outprefix(fem_out) display
capture set rmsg off
capture quietly log close

/// jlhdfe_alt y `X'  x_end         , absorb(id1 id2 id3 id4) vce(cluster id1 id2 id3 id4) method(cuda) threads(48) outprefix(fem_out) display
// jlhdfe_alt y `X' (x_end = z_end), absorb(id1 id2 id3 id4) vce(cluster id1 id2 id3 id4) method(cuda) threads(48) outprefix(fem_out) display

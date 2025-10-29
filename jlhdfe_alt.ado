program define jlhdfe_alt, eclass
	version 15.1
	// Stata-style front-end for fem_cli.jl with [if] [in]

	// --- PARSE OPTIONS ---
	syntax anything(name=spec) [if] [in] , ///
		[ absorb(varlist) vce(string) ///
		method(name) threads(integer 0) outprefix(name) pqfile(string) ///
		julia(string) femcli(string) display ]

	di as txt ">> jlhdfe_alt: options parsed OK"

	// defaults
	if "`method'"    == "" local method    "CUDA"
	if "`outprefix'" == "" local outprefix "fem_out"
	if "`pqfile'"    == "" local pqfile    "temp.parquet"
	if "`julia'"     == "" local julia     "julia"
	if "`femcli'"    == "" local femcli    "fem_cli.jl"

	* --- Resolve fem_cli.jl path (option > adopath > next to ado) ---
	if "`femcli'" != "" {
		capture confirm file "`femcli'"
		if _rc {
			di as err "femcli(`femcli') not found; provide a full path."
			exit 601
		}
	}
	else {
		quietly capture findfile fem_cli.jl
		if !_rc {
			local femcli "`r(fn)'"
		}
		else {
			quietly capture findfile jlhdfe_alt.ado
			if !_rc {
				local __pkgdir = subinstr("`r(fn)'","jlhdfe_alt.ado","",.)
				capture confirm file "`__pkgdir'fem_cli.jl"
				if !_rc local femcli "`__pkgdir'fem_cli.jl"
			}
		}
		if "`femcli'" == "" {
			di as err "Could not locate fem_cli.jl on adopath."
			di as err "Either net install jlhdfe_alt (so fem_cli.jl lands in PLUS) or pass femcli(/full/path/to/fem_cli.jl)."
			exit 601
		}
	}

	di as txt ">> fem_cli.jl: " as res "`femcli'"

	// 1) Split spec into depvar and RHS (RHS may include (endog = instr))
	local lhsrhs `spec'
	gettoken y lhsrhs : lhsrhs
	confirm variable `y'

	// Extract IV block "( ... = ... )" if present
	local ivblock ""
	if regexm(`"`lhsrhs'"', "\(([^)]*)\)") {
		local ivblock `"`=regexs(1)'"'
		local lhsrhs : subinstr local lhsrhs "(`ivblock')" "", all
	}

	// Separate plain vars vs factor terms
	local W ""
	local FTERMS ""
	local nt : word count `lhsrhs'
	forvalues i = 1/`nt' {
		local tok : word `i' of `lhsrhs'
		if strpos("`tok'", ".") | strpos("`tok'", "#") {
			if "`FTERMS'" == "" local FTERMS "`tok'"
			else               local FTERMS "`FTERMS' `tok'"
		}
		else if "`tok'" != "" {
			if "`W'" == "" local W "`tok'"
			else           local W "`W' `tok'"
		}
	}

	// vce(): accept vce(robust), vce(cluster ntrab), vce(cluster(ntrab estab_id))
	local cluster ""
	if "`vce'" != "" {
		// lowercase safely, then normalize whitespace
		local vcel = lower("`vce'")
		local vcel : list retokenize vcel     // ← name, not `vcel'

		if strpos("`vcel'","robust") {
			local cluster ""
		}
		else if regexm("`vcel'","cluster[[:space:]]*\(([^)]*)\)") {
			local cluster = trim(regexs(1))
		}
		else if regexm("`vcel'","cluster[[:space:]]+(.+)$") {
			local cluster = trim(regexs(1))
		}

		if "`cluster'" != "" {
			local cluster : list retokenize cluster   // ← name, not `cluster'
			local cluster : subinstr local cluster " " ",", all
			local cluster : subinstr local cluster ",," ",", all
		}
	}

	// FE from absorb()
	local FEplus ""
	if "`absorb'" != "" {
		local _tmp : subinstr local absorb " " ") + fe(", all
		local FEplus "fe(`_tmp')"
	}

	// IV: accept space- or plus-separated on each side, emit plus-separated for Julia
	local IVSTR ""
	if "`ivblock'" != "" {
		// split raw block at '='
		local left  ""
		local right ""
		gettoken left right : ivblock, parse("=")
		local right : subinstr local right "=" "", all

		// trim using function form (not an extended macro function)
		local left  = trim("`left'")
		local right = trim("`right'")

		// If a side has no '+', treat whitespace as separators and insert ' + '
		if strpos("`left'","+")==0 {
			local left : list retokenize left
			local left : subinstr local left " " " + ", all
		}
		if strpos("`right'","+")==0 {
			local right : list retokenize right
			local right : subinstr local right " " " + ", all
		}

		// final IV chunk for the formula
		local IVSTR `"`left' ~ `right'"'
	}

	// Build RHS parts cleanly (no leading " + ")
	local Wplus ""
	if "`W'" != "" {
		local Wplus : subinstr local W " " " + ", all
	}
	local Fplus ""
	if "`FTERMS'" != "" {
		local Fplus : subinstr local FTERMS " " " + ", all
	}

	local RHS ""
	if "`IVSTR'" != "" {
		local part "(`IVSTR')"
		if "`RHS'" != "" local RHS "`RHS' + `part'"
		else              local RHS "`part'"
	}
	foreach part in Wplus Fplus FEplus {
		if "``part''" != "" {
			if "`RHS'" != "" local RHS "`RHS' + ``part''"
			else              local RHS "``part''"
		}
	}
	local FORMULA "`y' ~ `RHS'"

	// Derive xvars / zvars from ivblock (for keeping)
	local XVARS ""
	local ZVARS ""
	if "`ivblock'" != "" {
		local left ""
		local right ""
		gettoken left right : ivblock, parse("=")
		local right : subinstr local right "=" "", all

		foreach side in left right {
			local __s ``side''
			local __s : subinstr local __s "+" " ", all
			local __s : subinstr local __s "#" " ", all
			local __s : subinstr local __s "(" " ", all
			local __s : subinstr local __s ")" " ", all
			local tmp ""
			local m : word count `__s'
			forvalues j = 1/`m' {
				local t : word `j' of `__s'
				if "`t'" == "" continue
				local base "`t'"
				local p = strpos("`t'", ".")
				if `p' > 0 {
					local L = strlen("`t'")
					local base = substr("`t'", `p'+1, `L' - `p')
				}
				capture confirm variable `base'
				if !_rc local tmp "`tmp' `base'"
			}
			local tmp : list uniq tmp
			if "`side'" == "left"  local XVARS "`tmp'"
			if "`side'" == "right" local ZVARS "`tmp'"
		}
	}

	// --- collect base variables from factor/interaction terms in FTERMS ---
	local FBASES ""
	if "`FTERMS'" != "" {
		// Flatten operators to spaces and strip parentheses
		local _flat "`FTERMS'"
		local _flat : subinstr local _flat "##" " ", all
		local _flat : subinstr local _flat "#"  " ", all
		local _flat : subinstr local _flat "+"  " ", all
		local _flat : subinstr local _flat "("  " ", all
		local _flat : subinstr local _flat ")"  " ", all

		// Now tokens look like: c.Z_mt  ib3.k_new  i.female  c.age
		local __n : word count `_flat'
		forvalues __i = 1/`__n' {
			local __tok : word `__i' of `_flat'
			if "`__tok'" == "" continue

			// peel prefix up to last dot if present (c., i., ib3., ib(3)., etc.)
			local __base "`__tok'"
			local __p = strpos("`__tok'", ".")
			if `__p' > 0 {
				local __base = substr("`__tok'", `__p'+1, .)
			}

			// keep only if this is a real variable in the data
			capture confirm variable `__base'
			if !_rc local FBASES `FBASES' `__base'
		}
		// de-dup
		local FBASES : list uniq FBASES
	}


	// Restrict sample [if][in] and export parquet
	local KEEPLIST `y' `W' `FBASES' `absorb' `XVARS' `ZVARS'
	tokenize "`cluster'", parse(",")
	if "`1'" != "" & "`2'" == "" local KEEPLIST `KEEPLIST' `1'

	// prepare cluster positional arg (literal "" if empty)
	local CLARG = trim("`cluster'")
	local CLARGQ `""""'
	if `"`CLARG'"' != "" local CLARGQ `"`CLARG'"'

	preserve
    // Apply [in] then [if]
	// Apply [in] then [if]
        if "`in'" != "" keep `in'
	if "`if'" != "" keep `if'

        // Trim to variables we'll send to Julia
        if "`KEEPLIST'" != "" keep `KEEPLIST'
        compress

        // parquet export (ssc install parquet)
        capture noisily {
		set rmsg on
		pq save using "`pqfile'", replace
		set rmsg off
        }
        if _rc {
		di as err "Parquet write failed. (ssc install parquet)"
		error _rc
        }

        // write formula
        tempname fh
        file open `fh' using "formula.txt", write replace text
        file write `fh' `"`FORMULA'"' _n
        file close `fh'
        confirm file "formula.txt"

        // build + run Julia cmd
        local CMD ""
        if `threads' > 0 local CMD `"JULIA_NUM_THREADS=`threads' "'
        local CMD `"`CMD'`julia' "`femcli'" "`pqfile'" "formula.txt" `CLARGQ' "`method'" "`outprefix'" 0"'
        di as txt ">> " as res `"`CMD'"'
        ! `CMD'

	// === collect results into matrices (inside preserve) ===
	quietly count
	local N = r(N)

	// read coefficients CSV produced by fem_cli
	capture noisily import delimited using "`outprefix'_coef.csv", clear varnames(1)
	if _rc {
		di as err "Could not read `outprefix'_coef.csv'."
		restore
		error _rc
	}

	// remove rows with missing coeff/SE (collinear/dropped terms)
	confirm variable term
	confirm variable estimate
	confirm variable stderr

	local Nraw = _N
	quietly levelsof term if missing(estimate) | missing(stderr), local(_dropped) clean
	if "`_dropped'" != "" {
		local ndrop : word count `_dropped'
		di as txt "note: dropping `ndrop' collinear term(s): " as res "`_dropped'"
	}
	quietly drop if missing(estimate) | missing(stderr)

	count
	if r(N) == 0 {
		restore
		di as err "All coefficients are missing (likely all RHS collinear with FE)."
		error 504
	}

	ereturn local dropped "`_dropped'"

	// Build b and V (diagonal from stderr)
	order term estimate stderr
	tempname _B _SE
	quietly mkmat estimate, matrix(`_B') rowname(term)
	quietly mkmat stderr,   matrix(`_SE') rowname(term)

	// Build clean coefficient names (unique, valid) & special-case intercept
	local rawnames : rownames `_B'
	local clean ""
	foreach nm of local rawnames {
		local s "`nm'"
		if "`nm'" == "(Intercept)" local s "_cons"
		else {
			local s : subinstr local s " " "_" , all
			local s : subinstr local s ":" "_" , all
			local s : subinstr local s "(" "_" , all
			local s : subinstr local s ")" ""  , all
			local s : subinstr local s "." "_" , all
			// if first char not a letter, prefix
			if !regexm("`s'","^[A-Za-z]") local s = "b_`s'"
		}
		local clean `clean' `s'
	}
	// ensure uniqueness just in case
	local clean : list uniq clean

	// b: 1×k row vector with colnames + eq stripe
	tempname Bvec Vmat
	matrix `Bvec' = `_B''
	matrix colnames `Bvec' = `clean'
	matrix coleq    `Bvec' = eq1

	// V: k×k diagonal with matching row/col names + same eq stripe
	local k = rowsof(`_SE')
	matrix `Vmat' = J(`k', `k', 0)
	forvalues i = 1/`k' {
		matrix `Vmat'[`i',`i'] = (`_SE'[`i',1])^2
	}
	matrix rownames `Vmat' = `clean'
	matrix colnames `Vmat' = `clean'
	matrix roweq    `Vmat' = eq1
	matrix coleq    `Vmat' = eq1

	restore

	// === post to e() ===
	ereturn clear
	capture noisily ereturn post `Bvec' `Vmat', obs(`N')
	if _rc {
		// last-resort fallback: anonymize names to guarantee posting
		tempname Bbare Vbare
		matrix `Bbare' = `Bvec'
		matrix `Vbare' = `Vmat'
		matrix colnames `Bbare' =
		matrix rownames `Vbare' =
		matrix colnames `Vbare' =
		matrix coleq    `Bbare' = eq1
		matrix roweq    `Vbare' = eq1
		matrix coleq    `Vbare' = eq1
		ereturn clear
		ereturn post `Bbare' `Vbare', obs(`N')
	}
	// tag metadata
	local vcedesc = cond("`cluster'"=="", "robust", "cluster(`cluster')")
	ereturn local cmd      "jlhdfe_alt"
	ereturn local depvar   "`y'"
	ereturn local title    "FixedEffectModels (Julia) via fem_cli"
	ereturn local method   "`method'"
	ereturn local absorb   "`absorb'"
	ereturn local vce      "`vcedesc'"
	ereturn local cluster  "`cluster'"
	ereturn local outfix   "`outprefix'"
	ereturn local formula  `"`FORMULA'"'

	// optional display of Julia outputs (summary + a peek at coefficients)
	if "`display'" != "" {
		di as txt "----- fem_cli summary (first 80 lines) -----"
		capture confirm file "`outprefix'_summary.txt"
		if _rc {
			di as err "(no summary file `outprefix'_summary.txt' found)"
		}
		else {
			tempname fh
			file open `fh' using "`outprefix'_summary.txt", read text
			file read `fh' __line
			local shown = 0
			while (r(eof)==0 & `shown' < 80) {
				di as txt `"`__line'"'
				local ++shown
				file read `fh' __line
			}
			file close `fh'
			if `shown' == 80 di as txt "... [truncated]"
		}
	}
end

*! version 0.8.0 28mar2024
program eventbaseline, eclass
    version 18
    syntax [, pre(integer 1) post(integer 3) baseline(string) generate(string) level(real 95) weight(varname)] [graph]
	if ("`level'" == "") {
		local level 95
	}    
    if ("`baseline'" == "") {
        local baseline "-1"
    }   
    if (e(cmd) != "xthdidregress") {
        display in red "eventbaseline can only be used after xthdidregress"
        error 198
    }
    local T1 = `pre'-1
    local K = `pre'+`post'+1

    local depvar = e(depvar)
    local cohortvar = e(cohortvar)
    local timevar = e(timevar)

    * reweight cohorts if needed
    levelsof `cohortvar' if `cohortvar' != 0, local(gs)
    levelsof `timevar', local(ts)
    local G : word count `gs'
    local T : word count `ts'

    tempname w bad_coef bad_Var Wcum Wsum D W0 b V Nobs coefplot tlabels

    matrix `w' = J(`G', 1, .)
    forvalues g = 1/`G' {
        local i = `g' + 1
        local cohort : word `g' of `gs'
        if "`weight'" == "" {
            matrix `w'[`g', 1] = e(cohort)[`i', 1]
        }
        else {
            summarize `weight' if `cohortvar' == `cohort', meanonly
            matrix `w'[`g', 1] = r(mean)
        }
    }

    matrix `bad_coef' = e(b)
    matrix `bad_Var' = e(V)
    local GT = colsof(`bad_coef')

    assert `GT' == `G' * `=`T'-1'
    assert colsof(`bad_Var') == `GT'

    matrix `Wcum' = J(`GT', `K', 0)
    local i = 1
    forvalues g = 1/`G' {
        forvalues t = 2/`T' {
            local time : word `t' of `ts'
            local start : word `g' of `gs'
            local e = `time' - `start'
            if inrange(`e', -`pre', `post') {
                matrix `Wcum'[`i', `e' + `pre' + 1] = `w'[`g', 1]
            }
            local i = `i' + 1
        }
    }
    matrix `Wsum' = J(1, `GT', 1) * `Wcum' 
    matrix `D' = diag(`Wsum')
    matrix `Wcum' = `Wcum' * inv(`D')

    tempvar exclude esample
    * exclude observations outside of the event window
    quietly generate `exclude' = cond(`cohortvar' == 0, 0, !inrange(`timevar' - `cohortvar', -`pre', `post'))
    quietly generate `esample' = e(sample) & (`exclude' == 0)
    quietly count if `esample'
    local Nobs = r(N)

    if ("`baseline'" == "average") {
        matrix `W0' = I(`K') - (J(`K', `pre', 1/`pre'), J(`K', `post'+1, 0))
    }
    else if ("`baseline'" == "atet") {
        matrix `W0' = (J(1, `pre', -1/`pre'), J(1, `post'+1, 1/(`post'+1)))
    }
    else {
        if (!inrange(`baseline', -`pre', -1)) {
            display in red "Baseline must be between -`pre' and -1"
            error 198
        }
        matrix `W0' = I(`K')
        local bl = `pre' + `baseline' + 1
        forvalues i = 1/`K' {
            matrix `W0'[`i', `bl'] = `W0'[`i', `bl'] - 1.0
        }
    }
    matrix `b' = `bad_coef' * `Wcum' * `W0''
    matrix `V' = `W0' * `Wcum'' * `bad_Var' * `Wcum' * `W0''

    if ("`baseline'" == "atet") {
        local colnames "ATET"
    }
    else {
        * label coefficients
        forvalues t = -`pre'/`post' {
            local colnames `colnames' `t'
        }
    }
    matrix colname `b' = `colnames'
    matrix colname `V' = `colnames'
    matrix rowname `V' = `colnames'

    matrix `coefplot' = J(`K', 4, .)
    matrix colname `coefplot' = xvar b ll ul
    local tlabels ""
    forvalues t = -`pre'/`post' {
        local tlabels `tlabels' `t'
        local i = `t' + `pre' + 1
        matrix `coefplot'[`i', 1] = `t''
        matrix `coefplot'[`i', 2] = `b'[1, `i']
        matrix `coefplot'[`i', 3] = `b'[1, `i'] + invnormal((100-`level')/200) * sqrt(`V'[`i', `i'])
        matrix `coefplot'[`i', 4] = `b'[1, `i'] - invnormal((100-`level')/200) * sqrt(`V'[`i', `i'])
    }

    tempname coef lower upper
    if ("`generate'" != "") {
        capture frame drop `generate'
        frame create `generate' time coef lower upper
        forvalues t = -`pre'/`post' {
            local i = `t' + `pre' + 1
            scalar `coef' = `b'[1, `i']
            scalar `lower' = `b'[1, `i'] + invnormal((100-`level')/200) * sqrt(`V'[`i', `i'])
            scalar `upper' = `b'[1, `i'] - invnormal((100-`level')/200) * sqrt(`V'[`i', `i'])
            frame post `generate' (`t') (`coef') (`lower') (`upper')
        }
        frame `generate': tsset time
        frame `generate': format coef lower upper %9.3f
    }

	ereturn post `b' `V', obs(`Nobs') esample(`esample')
    ereturn local depvar `depvar'
	ereturn local cmd eventstudy
	ereturn local cmdline eventstudy `0'

    _coef_table_header, title(Event study relative to `baseline') width(62)
	display
	_coef_table, bmat(e(b)) vmat(e(V)) level(`level') 	///
		depname(`depvar') coeftitle(ATET)

    if ("`graph'" == "graph") {
        hetdid_coefplot, mat(`coefplot') title(Event study relative to `baseline') ///
            ylb(`depvar') xlb("Length of exposure to the treatment") ///
            yline(0) legend(off) level(`level') yline(0,  extend) ytick(0, add) ylabel(0, add) xlabel(`tlabels')
    }
end


{smcl}
{title:jlhdfe_alt} — Stata front-end for Julia FixedEffectModels

{p 4 4 2}
Syntax:
{p_end}
{p 8 8 2}
{cmd:jlhdfe_alt} {it:depvar} [(endog = instr)] [controls and factor terms] {cmd:,} {cmd:absorb(}{it:fevars}{cmd:)} {cmd:vce(}{it:robust|cluster ...}{cmd:)} [method(cpu|cuda) threads(#) outprefix(name) femcli(path) display]
{p_end}

{title:Examples}
{p 8 8 2}{cmd:jlhdfe_alt lw age age2, absorb(ano) vce(cluster ntrab) method(cuda) display}{p_end}

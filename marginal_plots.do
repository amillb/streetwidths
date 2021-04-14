// Regressions and marginal effects plots for street widths
// Adam Millard-Ball, December 2020

// File created in Python
use "/Users/adammb/Documents/GDrive/Research/Street widths/Working/streets_residential.dta", clear
set more off
set scheme s1mono
gen fbuy_sq = fbuy*fbuy
encode fips, gen(county)

log using "/Users/adammb/Documents/GDrive/Research/Street widths/Working/stata_log.smcl", replace
// continuous (quadratic) year
regress median_width pop_sqkm length c.fbuy##c.fbuy  i.county, robust beta
stdBeta, se

// indicator variable for year, without landvalue
regress median_width pop_sqkm length i.fbuy i.county, robust beta
stdBeta, se

// with landvalue: preferred specification
regress median_width pop_sqkm length landvalue_acre_std i.fbuy i.county, robust beta
stdBeta, se
log close

// marginal effects plots
margins, dydx(fbuy) atmeans
marginsplot,  xlabel(4 "1825" 9 "1850" 14 "1875" 19 "1900" 24 "1925" 29 "1950" 34 "1975" 39 "2000", grid) title("") xti("") yti("Marginal effect on street width (meters)")
graph export "/Users/adammb/Documents/GDrive/Research/Street widths/Figures/marginaleffects_byyear.pdf", replace

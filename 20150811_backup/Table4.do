/*

This will estimate the income elasticity of migration, as reported in Table 4

*/

/* Load data on migration shares, real income, and geographic bilateral distance */
use migration_data, clear

/* Define necessary variables for the regressions */
gen lvjvi=ln(Vj/Vi)
gen ld=ln(distance)
gen logmijmii=ln(mij_mii2000)
gen logmijmii2005=ln(mij_mii2005)
gen logmij2000=ln(mij2000)
gen logmij2005=ln(mij2005)
egen FE=group(importer exporter)
egen origin=group(i)
gen lrYLj=ln(Vj)
gen IVj=ln(rYLj_bartik_sector) /* the bartik-style instrument of expected income, 2005 Census Income */
gen IVj_uhs2000=ln(rYLj_bartik_sector_uhs2000) /* the bartik-style instrument of expected income, 2000 Urban Household Survey */
gen IVj_neighbour=ln(rYLj_IVneighbour)

/* OLS: Columns 1 and 2 */
est clear
qui regress logmijmii lrYLj ld i.origin 
est store r1
qui regress logmijmii lrYLj i.FE i.origin 
est store r2
esttab r*, se lab r2 nogap bracket b(2) star(* .1 ** .05 *** .01) keep(lrYLj ld _cons)

/* IV: Columns 3 and 4 */
est clear
qui ivregress 2sls logmij2000 (lrYLj = IVj_neighbour) ld i.origin
est store r1
qui ivregress 2sls logmij2000 (lrYLj = IVj_neighbour) i.FE i.origin
est store r2
esttab r*, se lab r2 nogap bracket b(2) star(* .1 ** .05 *** .01) keep(lrYLj ld _cons)
ivregress 2sls logmij2000 (lrYLj = IVj_neighbour) ld i.origin, first
ivregress 2sls logmij2000 (lrYLj = IVj_neighbour) i.FE i.origin, first

/* IV: Columns 5 and 6 */
est clear
qui ivregress 2sls logmijmii2005 (lrYLj = IVj) ld i.origin
est store r1
qui ivregress 2sls logmijmii2005 (lrYLj = IVj) i.FE i.origin
est store r2
esttab r*, se lab r2 nogap bracket b(2) star(* .1 ** .05 *** .01) keep(lrYLj ld _cons)
ivregress 2sls logmijmii2005 (lrYLj = IVj) ld i.origin, first
ivregress 2sls logmijmii2005 (lrYLj = IVj) i.FE i.origin, first

/* IV: Columns 7 and 8 - Subset of provinces since this is UHS instead of Census */
est clear
qui ivregress 2sls logmijmii (lrYLj = IVj_uhs2000) ld i.origin
est store r1
qui ivregress 2sls logmijmii (lrYLj = IVj_uhs2000) i.FE i.origin
est store r2
esttab r*, se lab r2 nogap bracket b(2) star(* .1 ** .05 *** .01) keep(lrYLj ld _cons)
ivregress 2sls logmijmii (lrYLj = IVj_uhs2000) ld i.origin, first
ivregress 2sls logmijmii (lrYLj = IVj_uhs2000) i.FE i.origin, first

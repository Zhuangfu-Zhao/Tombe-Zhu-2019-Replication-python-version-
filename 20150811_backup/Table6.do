/*

This estimates and reports the regional trade costs for Table 6
It also estimates the changes in those costs between 2002 and 2007

*/

/* Load the Trade Share and Trade Cost Data */
use trade_data, clear

/* Construct Asymmetric Trade Costs, Given Symmetric HRI Costs tij */
gen tij_asym_ag=tij_ag*(ti_ag/tj_ag)^0.5
gen tij_asym_na=tij_na*(ti_na/tj_na)^0.5

/* Construct the relative change in trade costs from 2002 to 2007 */
keep importer exporter tij* year i j Eij_ag Eij_na Ei_ag Ei_na
reshape wide tij_asym_ag tij_asym_na tij_ag tij_na Eij_ag Eij_na Ei_ag Ei_na, i(importer exporter) j(year)
gen dij_ag=tij_ag2007/tij_ag2002
gen dij_na=tij_na2007/tij_na2002
gen dij_asym_ag=tij_asym_ag2007/tij_asym_ag2002
gen dij_asym_na=tij_asym_na2007/tij_asym_na2002

/* Table 6: expenditure wegithed average of sectoral trade costs and their changes */
gen tij_asym_avg2002=(tij_asym_ag2002*Eij_ag2002+tij_asym_na2002*Eij_na2002)/(Eij_ag2002+Eij_na2002)
table i j, c(mean tij_asym_avg2002)
gen dij_asym_avg=(dij_asym_ag*Eij_ag2002+dij_asym_na*Eij_na2002)/(Eij_ag2002+Eij_na2002)
table i j, c(mean dij_asym_avg)


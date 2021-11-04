# Endothelial-cell-sprouting-model

1. Cell proliferation module 
Codes can be found under ‘CellProliferation’ folder.
bestparams.m - Select the best fits from the total fits
coreFile_moleculardetail.m - Odes for the molecularly detailed ERK-Akt model  
findParamValue.m - Generate random numbers that follow a defined normal distribution
plotfits.m - Plot the fitting results
RunModelrcp.m - Calculate the cell proliferation rates
shadedErrorBar.m - Plot shaded error bars

initvalue.mat - Initial concentrations of the ERK-Akt model 
params.mat - Parameter values of the ERK-Akt model 
results.mat - All fitting results
bestfits.mat - 21 sets of the best fitted parameter values
besterrors.mat - 21 sets of the best fitted parameter errors

r_g_V.mat - Cell proliferation rates stimulated by VEGF for all fits 
r_g.mat - Cell proliferation rates stimulated by FGF for all fits 
r_p_V.mat - Cell proliferation rates contributed by ERK phosphorylation stimulated by VEGF for all fits 
r_p.mat - Cell proliferation rates contributed by ERK phosphorylation stimulated by FGF for all fits 
r_s_V.mat - Cell proliferation rates contributed by Akt phosphorylation stimulated by VEGF for all fits 
r_s.mat - Cell proliferation rates contributed by Akt phosphorylation stimulated by FGF for all fits 

best_r_g_V.mat - Cell proliferation rates stimulated by VEGF for the best fits 
best_r_g.mat - Cell proliferation rates stimulated by FGF for the best fits 
best_r_p_V.mat - Cell proliferation rates contributed by ERK phosphorylation stimulated by VEGF for the best fits 
best_r_p.mat - Cell proliferation rates contributed by ERK phosphorylation stimulated by FGF for the best fits 
best_r_s_V.mat - Cell proliferation rates contributed by Akt phosphorylation stimulated by VEGF for the best fits 
best_r_s.mat - Cell proliferation rates contributed by Akt phosphorylation stimulated by FGF for the best fits 

2. Cell sprouting model 
Codes can be found under ‘SproutingModel’ folder.
Cellresponses.m - Calculate TL, NS, and AL formed from 500-cell spheroids stimulated by VEGF 
coreFile_moleculardetail.m - Odes for the molecularly detailed ERK-Akt model  
findParamValue.m - Generate random numbers that follow a defined normal distribution
RunRates.m - Calculate the cell proliferation rates, sprout growth rates, and the probability of forming a new sprout stimulated by VEGF
shadedErrorBar.m - Plot shaded error bars

initvalue.mat - Initial concentrations of the ERK-Akt model 
params.mat - Parameter values of the ERK-Akt model 
bestfits.mat - 18 sets of best fitted parameter values
besterrors.mat - 18 sets of best fitted parameter errors

P_V.mat - The probabilities of forming a new sprout stimulated by VEGF for the best fits 
rg_V.mat - The cell proliferation rates stimulated by VEGF for the best fits 
rm_V.mat - The sprout growth rates stimulated by VEGF for the best fits 

All_AvgLength_VEGF.mat - AL formed from 500-cell spheroids stimulated by VEGF
All_TotalLength_VEGF.mat - TL formed from 500-cell spheroids stimulated by VEGF
All_TotalTipCell_VEGF.mat - NS formed from 500-cell spheroids stimulated by VEGF


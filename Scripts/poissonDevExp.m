function devExp = poissonDevExp(Ydata,Yhat)
%% Calculate deviance explained for poisson GLM

Dmodel = 2*nansum(Ydata.*log(Ydata./Yhat)-(Ydata-Yhat));

Dnull = 2*nansum(Ydata.*log(Ydata/mean(Ydata))-(Ydata-Yhat));

devExp = 1 - Dmodel/Dnull;
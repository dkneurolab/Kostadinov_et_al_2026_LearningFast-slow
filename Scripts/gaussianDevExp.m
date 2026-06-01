function devExp = gaussianDevExp(Ydata,Yhat)
%% Calculate deviance explained for gaussian GLM

SSR = sum((Ydata-Yhat).^2,'omitnan');

SST = sum((Ydata-mean(Ydata)).^2,'omitnan');

devExp = 1-SSR/SST;

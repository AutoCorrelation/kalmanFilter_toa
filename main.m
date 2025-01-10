% https://github.com/AutoCorrelation/Sms2024
if ~exist("obj","var")
    obj = simulate(1e4, 11);
    obj = setsquareAnchorPos(obj, 10);
    obj = generateRangingInfo(obj);
    obj = toa(obj);
    obj = getQ(obj);
    obj = getR(obj);
    obj = ToaKalmanFilter(obj);
end
obj = ToaUKF(obj);
obj = resultPlot(obj);
% obj = optimizeParam1(obj);
% obj = optimizeParam2(obj);
% obj = optimizeParam3(obj);
obj = optimizeParam4(obj);



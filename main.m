% https://github.com/AutoCorrelation/Sms2024
obj = [];
if isempty(obj)
    obj = simulate(1e3, 11);
    obj = setsquareAnchorPos(obj, 10);
    obj = generateRangingInfo(obj);
    obj = toa(obj);
    obj = getQ(obj);
    obj = getR(obj);
end
obj = ToaKalmanFilter(obj);
obj = resultPlot(obj);
obj = optimizeParam(obj);

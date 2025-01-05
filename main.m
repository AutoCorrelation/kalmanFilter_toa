% https://github.com/AutoCorrelation/Sms2024

if isempty(obj)
    obj = simulate(1e4, 11);
    obj = setsquareAnchorPos(obj, 10);
    obj = generateRangingInfo(obj);
    obj = toa(obj);
    obj = getQ(obj);
    obj = getR(obj);
end
obj = ToaKalmanFilter(obj);
obj = resultPlot(obj);

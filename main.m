% https://github.com/AutoCorrelation/Sms2024

obj = simulate(1e3, 11);
obj = setsquareAnchorPos(obj, 10);
obj = generateRangingInfo(obj);
obj = toa(obj);
obj = getQ(obj);
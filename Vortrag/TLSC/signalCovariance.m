function[c] = signalCovariance(a,c0,t1,t2)
    c = c0 * exp(-a^2 * (t1-t2)^2);
end
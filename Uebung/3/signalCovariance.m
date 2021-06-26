function[K] = signalCovariance(d)
% signalCovariance function
a = 150;
K0 = 1.36;
K = K0 / (1 + (d/a)^2);
end
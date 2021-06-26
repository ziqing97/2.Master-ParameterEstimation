a = 10;
c0 = 10;
r = 1:100;
c = c0 * exp(-a^2 * r.^2);

figure
plot(r,c);

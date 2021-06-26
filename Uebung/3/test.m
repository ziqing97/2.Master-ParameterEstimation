a = 150;
K0 = 1.36;

d = 0:20;

K = K0 ./ (1 + (d./a).^2);
figure
plot(d,K)
xlabel('distance (degree)')
ylabel('signal covariance')
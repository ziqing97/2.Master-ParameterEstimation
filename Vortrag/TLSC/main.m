% clc
% clearvars
% close all

% data
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';

t0 = [0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10, 14, 20];

n = length(t);
p = length(t0);
 
e = rand(n,1) / 10;

% Parameter
% compare case 1,2,3 and case 4,2,5
input = 2;
switch input
    case 1
        a = 1;
        c0 = 0.01;
        string = 'a=1,c0=0.01';
    case 2
        a = 1;
        c0 = 1;
        string = 'a=1,c0=1';
    case 3
        a = 1;
        c0 = 100;
        string = 'a=1,c0=100';
    case 4
        a = 0.1;
        c0 = 1;
        string = 'a=0.1,c0=1';
    case 5
        a = 100;
        c0 = 1;
        string = 'a=100,c0=1';
end

% Begin 
A = [t,ones(7,1)]; 

% covariance 
Qs1s1 = zeros(n);
Qs1s2 = zeros(n,p);

for i = 1:n
    for j = 1:n
        Qs1s1(i,j) = signalCovariance(a,c0,t(i),t(j));
    end
end

for i = 1:n
    for j = 1:p
        Qs1s2(i,j) = signalCovariance(a,c0,t(i),t0(j));
    end
end

Qyy = Qs1s1 + diag(e);
Qsy = Qs1s2';

% results
x_dach = inv(A' * inv(Qyy) * A) * A' * inv(Qyy) * h;
s_dach = Qsy * inv(Qyy) * (h - A * x_dach);

h0_dach = [t0',ones(p,1)] * x_dach + s_dach;

h_dach = A * x_dach;

% plot
figure
plot(t,h,'ok');
hold on
plot([0:20],x_dach(1)*[0:20]+x_dach(2))
plot(t0,h0_dach,'*r')
title(string)
legend('origin points','adjusted line','interpolated point')
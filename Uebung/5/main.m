clc
close all
clear

load('sim_TWS_Greenland.mat')
t = data_real(:,1); 
s = data_real(:,2);
len = length(t);

% c
A = [t, ones(len,1), cos(2 * pi * t / 12), sin(2 * pi * t / 12)];
x_dach = (A'* A) \ A'* s;
s_dach = A * x_dach;

e_dach = s - A * x_dach;
sigma = e_dach' * e_dach / (len - length(x_dach));
Sigma_xx = sigma * inv(A'* A);
Sigma_yy = A * Sigma_xx * A';

t_new = [46;47;115;116;117;118];
len2 = length(t_new);
A_new = [t_new, ones(len2,1), cos(2 * pi * t_new / 12), sin(2 * pi * t_new / 12)];
s_new = A_new * x_dach;

figure
hold on
scatter(t,s,'o')
scatter(t,s_dach,'+')
scatter(t_new,s_new,'*')
xlabel('time')
ylabel('TWSC')
legend('Observation','Adjusted Observation','Filled Gap')

%% Prior
std_xx = diag(sqrt(Sigma_xx));
x_dist = mvnrnd(x_dach, std_xx',10^6);
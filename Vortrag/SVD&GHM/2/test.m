clear all
close all
clc

%%
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';

A = [t];
H = [A,h];

[U, S, V] = svd(H);
E = U(:,end) * S(end,end) * V(:,end)';

x_svd = - V(1:end-1,end) / V(end,end);
H_tilde = - [A,h] * V(:,end) * V(:,end)';
A_tilde = H_tilde(:,1:end-1);


h_svd = (A + A_tilde) * x_svd;
t_svd = H_tilde(:,1) + t;


hold on
%%
A = [t,ones(7,1)];
H = [A,h];

[U, S, V] = svd(H);
E = U(:,end) * S(end,end) * V(:,end)';

x_svd2 = - V(1:end-1,end) / V(end,end);
H_tilde = - [A,h] * V(:,end) * V(:,end)';
A_tilde = H_tilde(:,1:end-1);


h_svd2 = (A + A_tilde) * x_svd2;
t_svd2 = H_tilde(:,1) + t;

%%
figure
plot(t,h,'o')
plot(t_svd,h_svd)
hold on
plot(t_svd2,h_svd2)
ylim([0,20])
xlim([0,10])
legend('origin point','A=[t]', 'A =[t,ones(7,1)]')
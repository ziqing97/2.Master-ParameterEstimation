close all
clc
clearvars

%% Data
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';

%% SVD  (gleich wie LS)
A = [t,ones(7,1)];
[U, S, V] = svd(A,'econ');
x_svd = V * inv(S) * U' * h;
h_svd = A * x_svd;

%% SVD 2 (soll wie TLS gleich sein)
A = [t,ones(7,1)];
H = [A,h];

[U, S, V] = svd(H);
E = U(:,end) * S(end,end) * V(:,end)';

x_svd2 = - V(1:end-1,end) / V(end,end);
H_tilde = - [A,h] * V(:,3) * V(:,3)';
A_tilde = H_tilde(:,1:2);

h_svd12 = h + H_tilde(:,3);
h_svd2 = (A + A_tilde) * x_svd2;
t_svd2 = H_tilde(:,1) + t;

plot(t_svd2,h_svd2)

%% Normal LS
A = [t,ones(7,1)];
x_ls = (A' * A) \ A' * h; % a0,b0
h_ls = A * x_ls;

%% GHM ohne
P = diag(ones(14,1));

A = [t,ones(7,1)];
x0 = (A' * A) \ A' * h; % a0,b0
y = [h;t];
y0 = [h;t];
dx_dach = 1;
j = 1;
x_dach = x0;

xohne = zeros(2,8);
xohne(:,1) = x0;
while norm(dx_dach) > 10e-8
    dy = -(y - y0);
    w0 = y(1:7) - x_dach(1) * y(8:14) - x_dach(2);
    A = [-y(1:7), -ones(7,1)];    
    B1 = eye(7);
    B2 = eye(7) * (-x_dach(1));
    Bt = [B1,B2];
    N = [Bt * inv(P) * Bt' ,-A;
    -A', zeros(2)];

    w = w0 + Bt * dy;
    k = [w;zeros(2,1)];
    L = N \ k;
    lambda = L(1:7);
    dx_dach = L(8:9);
    x_dach = x_dach + dx_dach;
    
    xohne(:,j+1) = x_dach;
    
    e = -inv(P) * Bt' * lambda;
    y0 = y + e;
    j = j+1;
end

t_dach_gh1 = y0(8:14);
h_dach_gh1 = y0(1:7);
x_GHMohne = x_dach;
% 
% 
%% GHM mit
% Naehrungswerte
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
P = diag([1,2,3,4,5,6,7,7,6,5,4,3,2,1]);
P = diag([1,2,3,4,5,6,7,1,2,3,4,5,6,7]);

A = [t,ones(7,1)];
x0 = (A' * A) \ A' * h; % a0,b0
y = [h;t];
y0 = [h;t];
x_dach = x0;

xmit = zeros(2,8);
xmit(:,1) = x0;

dx_dach = [1;1];
j = 1;
while norm(dx_dach) > 10e-7
    dy = -(y - y0);
    w0 = y(1:7) - x_dach(1) * y(8:14) - x_dach(2);
    A = [-y(1:7), -ones(7,1)];    
    B1 = eye(7);
    B2 = eye(7) * (-x_dach(1));
    Bt = [B1,B2];
    N = [Bt * inv(P) * Bt' ,-A;
    -A', zeros(2)];

    w = w0 + Bt * dy;
    k = [w;zeros(2,1)];
    L = N \ k;
    lambda = L(1:7);
    dx_dach = L(8:9);
    x_dach = x_dach + dx_dach;
    xmit(:,j+1) = x_dach;
    
    e = -inv(P) * Bt' * lambda;
    y0 = y + e;
    j = j+1;
end

t_dach_gh2 = y0(8:14);
h_dach_gh2 = y0(1:7);
x_GHMmit = x_dach;

%% Plot
figure
hold on
title("Vergleich")
plot(t,h,"og");
plot(t,h_ls,"*-","color",'m')
plot(t_dach_gh2,h_dach_gh2,'*-',"color",'r');
plot(t_dach_gh1,h_dach_gh1,'*-',"color",'k');
plot(t,h_svd,'x-',"color",'b')
plot(t_svd2,h_svd2,"*-","color",'c')
legend("punkt","normal LS","ghm mit","ghm ohne","svd wie LS",'svd wie GHM')


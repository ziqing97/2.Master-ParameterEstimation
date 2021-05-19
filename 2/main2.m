
clc
clearvars

%% Data
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
figure
scatter(t,h)

%% SVD
t_mean = mean(t);
h_mean = mean(h);
M = [t, h];
Q = M' * M;
[U, S, V] = svd(Q);
% U = V
u1 = U(1,1);
u2 = U(2,1);
phi = atan2(u2 ,u1);

a_svd1 = tan(phi);

hold on
plot([0,12],[0,a_svd1 * 12])

x_svd = a_svd1;

%% SVD 2
A = t;
[U, S, V] = svd(A,'econ');
x_svd2 = V * inv(S) * U' * h;

%% SVD3
A = t;
H = [A, h];

[U, S, V] = svd(H,'econ');
E = U(:,end) * S(end,end) * V(:,end)';
e_dach = E(:,end);
h_dach_svd = h - e_dach;
A_dach = A - E(:,1);
t_dach_svd = A_dach;

% figure
% hold on
% title('svd')
% scatter(t,h)
% plot(t_dach_svd,h_dach_svd,'*-')

%%
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
A = t;
H = [A, h];

[U, S, V] = svd(H);
x = - V(1,2) / V(2,2);
%% GHM ohne
% Naehrungswerte
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
P = diag(ones(14,1));

A = t;
x0 = (A' * A) \ A' * h; % a0,b0
y = [h;t];
y0 = [h;t];
x_dach = x0;

dx_dach = 1;
j = 1;

xohne = zeros(2,8);
xohne(:,1) = x0;
while norm(dx_dach) > 10e-8
    dy = -(y - y0);
    w0 = y(1:7) - x_dach * y(8:14);
    A = -y(8:14);    
    B1 = eye(7);
    B2 = eye(7) * (-x_dach);
    Bt = [B1,B2];
    N = [Bt * inv(P) * Bt' ,-A;
    -A', 0];

    w = w0 + Bt * dy;
    k = [w;0];
    L = N \ k;
    lambda = L(1:7);
    dx_dach = L(8);
    x_dach = x_dach + dx_dach;
    
    xohne(:,j+1) = x_dach;
    
    e = -inv(P) * Bt' * lambda;
    y0 = y + e;
    j = j+1;
end

t_dach_gh1 = y0(8:14);
h_dach_gh1 = y0(1:7);
figure
hold on
title("GHM ohne stochastiche Fehler")
scatter(t,h);
plot(t_dach_gh1,h_dach_gh1,'*-');
plot(t_dach_svd,h_dach_svd,'o-')
legend("origin data","gh modell","svd")
x_GHMohne = x_dach;
% 
% 
%% GHM mit
% Naehrungswerte
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
P = diag([1,2,3,4,5,6,7,7,6,5,4,3,2,1]);
P = diag([1,1,1,1,1,1,1,2,2,2,2,2,2,2]);
P = diag([1,2,1,2,1,2,1,2,1,2,1,2,1,2]);

A = t;
x0 = (A' * A) \ A' * h; % a0,b0
y = [h;t];
y0 = [h;t];
x_dach = x0;

dx_dach = 1;
j = 1;

xohne = zeros(2,8);
xohne(:,1) = x0;
while norm(dx_dach) > 10e-8
    dy = -(y - y0);
    w0 = y(1:7) - x_dach * y(8:14);
    A = -y(8:14);    
    B1 = eye(7);
    B2 = eye(7) * (-x_dach);
    Bt = [B1,B2];
    N = [Bt * inv(P) * Bt' ,-A;
    -A', 0];

    w = w0 + Bt * dy;
    k = [w;0];
    L = N \ k;
    lambda = L(1:7);
    dx_dach = L(8);
    x_dach = x_dach + dx_dach;
    
    xohne(:,j+1) = x_dach;
    
    e = -inv(P) * Bt' * lambda;
    y0 = y + e;
    j = j+1;
end


t_dach_gh2 = y0(8:14);
h_dach_gh2 = y0(1:7);
% figure
% hold on
% title("GHM mit stochastiche Fehler")
% scatter(t,h);
% plot(t_dach_gh2,h_dach_gh2,'*-');
x_GHMmit = x_dach;
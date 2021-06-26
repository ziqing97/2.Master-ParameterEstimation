%% Initialisierung
clc
clear all
close all

%% Data
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
n = length(t);

A2 = t;
A1 = ones(7,1);
B = h;

e = [1,2,3,4,5,6,7,7,6,5,4,3,2,1]';
P = diag(e);
delta = [e(1:n), e(n+1:end)];

%% Step 1
% n1 = n2 = d = 1
H = [A1, A2, B];

[Q,R] = qr(H);
R11 = R(1,1);
R12 = R(1,2:3);
R22 = R(2:3,2:3);

C = delta'* delta;
Rc = chol(C);

%% Step 2
[~,~,V,~,~] = gsvd(R22,Rc);
[Qz,Rz] = qr(V);
V = Qz;
Z2 = V(:,2);
Z1 = R11 \ (-R12) * Z2;

%% Step 3
Rz = [Z1;Z2] / norm([Z1;Z2]);
% Rz = [Z1;Z2] * 400;
Y = Rz(1:2);
Tao = Rz(3);

X_dach_GSVD = -Y * inv(Tao);

R22_tilde = -R22 * Z2 * Z2';
R22_neu = R22 + R22_tilde;
R_neu = [R11,R12];
R_neu = [R_neu;zeros(2,1),R22_neu];
R_neu = [R_neu;zeros(4,3)];
H_neu = Q * R_neu;

t_svd = H_neu(:,2);
h_svd = H_neu(:,3);


%% GHM mit
% Naehrungswerte
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';


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

figure
hold on
title("Vergleich")
plot(t,h,"ok");
plot(t_svd,h_svd,"*-","color",'b')
plot(t_dach_gh2,h_dach_gh2,'x-',"color",'r');
legend('points','svd','gauss helmert')
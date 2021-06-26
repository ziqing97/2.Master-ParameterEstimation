%% Initialisierung
clc
% clear all
close all

%% Data
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';

A2 = t;
A1 = ones(7,1);
B = h;

%% Step 1
% n1 = n2 = d = 1
H = [A1, A2, B];

[Q,R] = qr(H);
R11 = R(1,1);
R12 = R(1,2:3);
R22 = R(2:3,2:3);

%% Step 2
[U,S,V] = svd(R22);
Z2 = V(:,2);
Z1 = R11 \ (-R12) * Z2;

%% Step 3
% Rz = [Z1;Z2] / norm([Z1;Z2]);
Rz = [Z1;Z2] * 400;
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

figure
hold on
title("Vergleich")
plot(t,h,"ok");
plot(t_svd,h_svd,"*-","color",'b')
plot(t_dach_gh1,h_dach_gh1,'x-',"color",'r');
legend('points','svd','gauss helmert')
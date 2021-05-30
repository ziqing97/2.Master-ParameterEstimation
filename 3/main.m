%% Initialisierung
clc
clear all
close all

%% Data
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';

A1 = t;
A2 = ones(7,1);
B = h;

%% Step 1
% n1 = n2 = d = 1
H = [A2, A1, B];

[Q,R] = qr(H);
R11 = R(1,1);
R12 = R(1,2:3);
R22 = R(2:7,2:3);

% Delta = ones(7,2);
% Delta(:,2) = 2 * ones(7,1);

% Delta = [1,4,9,4,5,6,7;
%          2,3,4,5,6,7,8]';
% 
% C = Delta' * Delta;
% Rc = chol(C);
% 
% D = Delta * Delta';
% e = eig(D);
% Rd = chol(D);
C = eye(2);
D = eye(7);
Rc = chol(C);
Rd = chol(D);

%% Step 2
% [T,W,Z,S1,S2] = gsvd(R22,Rc);
% r = rank([R22;Rc]);
% sigma1 = S1(1,1) / S2(1,1);
% sigma2 = S1(2,2) / S2(2,2);
% 
% Z(:,1) = Z(:,1) / norm(Z(:,1));
% Z(:,2) = Z(:,2) / norm(Z(:,2));
% 
% E22 = eye(2);

[U,S,V] = svd(R22);
Z2 = V(:,2);
Z1 = inv(R11) * (-R12) * Z2;

%% Step 3
[Qz, Rz] = qr([Z1;Z2]);
Rz = [Z1;Z2] / norm([Z1;Z2]);
Y = Rz(1:2);
Tao = Rz(3);

X_dach = -Y * inv(Tao);
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
H = [A1, A2, B];

[Q,R] = qr(H);
R11 = R(1,1);
R12 = R(1,2:3);
R22 = R(2:7,2:3);

C = eye(2);
Rc = chol(C);

D = eye(7);
Rd = chol(D);

%% Step 2
[T,W,Z,S1,S2] = gsvd(R22,Rc);
r = rank([R22;Rc]);
sigma1 = S1(1,1) / S2(1,1);
sigma2 = S1(2,2) / S2(2,2);

Z(:,1) = Z(:,1) / norm(Z(:,1));
Z(:,2) = Z(:,2) / norm(Z(:,2));

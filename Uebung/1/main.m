%% Initialisierung
clc
clear all
close all

% load data
load('t.mat')
load('l.mat')
T = T/10; % 10 durch Schwingungen
len = length(l);


% Nährungswert g0 aus Mittelwert
T_mean = mean(T);
l_mean = mean(l);
g0 = l_mean * 4 * pi^2 / T_mean^2;

%% GHM ohne Gewicht
% Beobachtungen
y = [T;l];
y_dach = y;

g_dach = g0; % Anfangswert Unbekannte

dg = 100; % Schwellwert fuer erste Schleife
Q = diag(ones(len*2,1)); 
P = inv(Q); % Gewicht

j1 = 0;
% Iteration
while abs(dg) > 1e-8
    dy = -(y - y_dach);
    w0 = y(1:len).^2 - 4 * pi^2 * y(len+1:end).^2 / g_dach;
    A = -4 * pi^2 * y(len+1:end) / g_dach;  % nach g ableiten
    Bt = zeros(len,len*2);
    for i = 1:len
        Bt(i,i) = 2 * y(i); % nach T ableiten
        Bt(i,i+len) = 4 * pi^2 / g_dach; % nach l ableiten
    end
    k1 = [Bt * inv(P) * Bt', -A;
         -A', 0];
    w = w0 + Bt * dy;
    k2 = [w;0];
    L = k1 \ k2;
    lambda = L(1:len);
    dg = lambda(end);
    g_dach = g_dach + dg;
    e = - inv(P) * Bt' * lambda;
    y_dach = y + e;
    j1 = j1+1;
end

% Ergebnis
g1 = g_dach;
% Varianz
v1 = e' * P * e / (len - 1);

%% GHM mit Gewicht
% Beobachtungen
y = [T;l];
y_dach = y;

g_dach = g0; % Anfangswert Unbekannte

dg = 100; % Schwellwert fuer erste Schleife
Q = diag([ones(1,25) * 0.2, ones(1,25) * 2.5e-3]);
P = inv(Q); % Gewicht

% Iteration
j2 = 0;
while abs(dg) > 1e-6
    dy = -(y - y_dach);
    w0 = y(1:len).^2 - 4 * pi^2 * y(len+1:end).^2 / g_dach;
    A = -4 * pi^2 * y(len+1:end) / g_dach;  % nach g ableiten
    Bt = zeros(len,len*2);
    for i = 1:len
        Bt(i,i) = 2 * y(i); % nach T ableiten
        Bt(i,i+len) = 4 * pi^2 / g_dach;
    end
    k1 = [Bt * inv(P) * Bt', -A;
         -A', 0];
    w = w0 + Bt * dy;
    k2 = [w;0];
    L = k1 \ k2;
    lambda = L(1:len);
    dg = lambda(end);
    g_dach = g_dach + dg;
    e = - inv(P) * Bt' * lambda;
    y_dach = y + e;
    j2 = j2+1;
end

% Ergebnis
g2 = g_dach;
% Varianz
v2 = e' * P * e / (len - 1);


%% SVD
m = T.^2;
n = l;

M = [n,m];

Q = M' * M;
[U,S,V] = svd(Q);
[U,S,V] = svd(M);

k = V(1,2) / V(2,2);

g3 = 4 * pi^2 / k;

% plot([0,u(1)],[0,u(2)])

%% SVD 2
[U,S,V] = svd(l,'econ');
k = V * inv(S) * U' * T.^2;
g4 = 4 * pi^2 / k;


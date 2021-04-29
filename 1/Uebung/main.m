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

%% Fehler nicht quantifiziert
% Beobachtungen
y = [T;l];
y_dach = y;

g_dach = g0; % Anfangswert Unbekannte

dg = 100; % Schwellwert fuer erste Schleife
Q = diag(ones(len*2,1)); 
P = inv(Q); % Gewicht

% Iteration
while dg > 10e-3
    dy = y - y_dach;
    w0 = T.^2 - 4 * pi^2 * l.^2 / g_dach;
    A = -4 * pi^2 * l / g_dach;
    Bt = zeros(len,len*2);
    for i = 1:len
        Bt(i,2*i-1) = 2 * T(i); % nach T ableiten
        Bt(i,2*i-1) = 4 * pi^2 / g_dach;
    end
    k1 = [Bt * P * Bt', -A;
         -A', 0];
    w = w0 + Bt * dy;
    k2 = [w;0];
    L = k1 \ k2;
    lambda = L(1:len);
    dg = lambda(end);
    g_dach = g_dach + dg;
    e = - P * Bt' * lambda;
    y_dach = y + e;
end

% Ergebnis
g1 = g_dach;
% Varianz
v1 = e' * P * e / (len - 1);

%% Fehler Stochastische Modell
% Beobachtungen
y = [T;l];
y_dach = y;

g_dach = g0; % Anfangswert Unbekannte

dg = 100; % Schwellwert fuer erste Schleife
Q = diag([ones(1,25) * 0.2, ones(1,25) * 2.5e-3]);
P = inv(Q); % Gewicht

% Iteration
while dg > 10e-3
    dy = y - y_dach;
    w0 = T.^2 - 4 * pi^2 * l.^2 / g_dach;
    A = -4 * pi^2 * l / g_dach;
    Bt = zeros(len,len*2);
    for i = 1:len
        Bt(i,2*i-1) = 2 * T(i); % nach T ableiten
        Bt(i,2*i-1) = 4 * pi^2 / g_dach;
    end
    k1 = [Bt * P * Bt', -A;
         -A', 0];
    w = w0 + Bt * dy;
    k2 = [w;0];
    L = k1 \ k2;
    lambda = L(1:len);
    dg = lambda(end);
    g_dach = g_dach + dg;
    e = - P * Bt' * lambda;
    y_dach = y + e;
end

% Ergebnis
g2 = g_dach;
% Varianz
v2 = e' * P * e / (len - 1);
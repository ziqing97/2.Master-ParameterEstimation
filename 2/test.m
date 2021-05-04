% Beobachtungen
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';

% Design Matrix
A = [t,ones(7,1)];

% Nährungswerte Unbekanntern
x0 = (A' * A) \ A' * h; % a0,b0
y = [h;t];
y0 = [h;t];

x_dach = x0;
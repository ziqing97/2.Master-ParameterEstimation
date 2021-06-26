close all
clc
clearvars

%% SVD
X = [3,4,5,6,7,8,9];
Y = [7,7,11,11,15,16,19];
G = ones(1,7) * 1/7;

figure
scatter(X,Y)
% axis equal
% xlim([0,12])
% ylim([0,20])
X_mean = mean(X);
Y_mean = mean(Y);
X = X - X_mean;
Y = Y - Y_mean;

M = [X', Y'];
Q = M' * M;
[U, S, V] = svd(Q);

% U = V
u1 = U(1,1);
u2 = U(2,1);
phi = atan2(u2 ,u1);

a = tan(phi);
b = Y_mean - a * X_mean; 
hold on
plot([0,12],[b,a*12+b])

xsvd = [a;b];

%% GHM ohne
% Naehrungswerte
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
P = diag(ones(14,1));

H = [t,ones(7,1)];
x0 = (H' * H) \ H' * h; % a0,b0
y = [h;t];
y0 = [h;t];
x_dach = x0;

dx_dach = 1;
j = 1;

xohne = zeros(2,8);
xohne(:,1) = x0;
while norm(dx_dach) > 10e-8
    dy = y - y0;
    w0 = h - x_dach(1) * t - x_dach(2);
    A = [-t, -ones(7,1)];
    Bt = zeros(7,14);
    for i = 1:7
        Bt(i,2*i-1) = 1; % nach h
        Bt(i,2*i)  = -x_dach(1);
    end
    N = [Bt * inv(P) * Bt' ,-A;
    -A', zeros(2)];

    w = w0 + Bt * dy;
    k = [w;zeros(2,1)];
    L = N \ k;
    lambda = L(1:7);
    dx_dach = L(8:9);
    x_dach = x0 + dx_dach;
    
    xohne(:,j+1) = x_dach;
    
    e = -inv(P) * Bt' * lambda;
    y0 = y0 + e;
    j = j+1;
end

% figure
% scatter(t,h);
hold on
plot([0,12],[x_dach(2),x_dach(1)*12+x_dach(2)])
% axis equal
% xlim([0,12])
% ylim([0,20])
x_ohne = x_dach;


%% GHM mit
% Naehrungswerte
t = [3,4,5,6,7,8,9]';
h = [7,7,11,11,15,16,19]';
P = diag([1,2,3,4,5,6,7,7,6,15,14,68,4,9]);
% P = diag([5,6,7,8,9,1,2,3,4,10,11,12,13,14]);

H = [t,ones(7,1)];
x0 = (H' * H) \ H' * h; % a0,b0
y = [h;t];
y0 = [h;t];
x_dach = x0;

xmit = zeros(2,8);
xmit(:,1) = x0;

dx_dach = 1;
j = 1;
while norm(dx_dach) > 10e-3
    dy = y - y0;
    w0 = h - x_dach(1) * t - x_dach(2);
    A = [-t, -ones(7,1)];
    Bt = zeros(7,14);
    for i = 1:7
        Bt(i,2*i-1) = 1; % nach h
        Bt(i,2*i)  = -x_dach(1);
    end
    N = [Bt * inv(P) * Bt' ,-A;
    -A', zeros(2)];

    w = w0 + Bt * dy;
    k = [w;zeros(2,1)];
    L = N \ k;
    lambda = L(1:7);
    dx_dach = L(8:9)
    x_dach = x0 + dx_dach;
    
    xmit(:,j+1) = x_dach;
    
    e = -inv(P) * Bt' * lambda;
    y0 = y0 + e;
    j = j+1;
end

% figure
% scatter(t,h);
hold on
plot([0,12],[x_dach(2),x_dach(1)*12+x_dach(2)])
% axis equal
% xlim([0,12])
% ylim([0,20])
x_mit = x_dach;


x_all = [[a;b],x_mit,x_ohne];

legend('measurement', 'svd', 'ohne', 'mit')

X = [3,4,5,6,7,8,9];


phi1 = atan(xsvd(1));
l1 = [cos(phi1);sin(phi1)];
Y = [7,7,11,11,15,16,19];
Y = Y - xsvd(2);
M = [X', Y'];
e1 = l1' * M' * M * l1;

phi2 = atan(xohne(1,2));
l2 = [cos(phi2);sin(phi2)];
Y = [7,7,11,11,15,16,19];
Y = Y - xohne(2,2);
M = [X', Y'];
e2 = l2' * M' * M * l2;

phi3 = atan(xmit(1,3));
l3 = [cos(phi3);sin(phi3)];
% Y = [7,7,11,11,15,16,19];
% Y = Y - xmit(2,3);
% M = [X', Y'];
e3 = l3' * M' * M * l3;


H = [t,ones(7,1)];
P = diag([1,2,3,4,5,6,7]);
x00 = (H' * inv(P)* H) \ H' * inv(P) * h; % a0,b0
hold on
plot([0,12],[x00(2),x00(1)*12+x00(2)])


t1 = [H, h];
[U1, S1, V1] = svd(t1,'econ');
E = U1(:,end) * S1(end,end) * V1(:,end)'
e_dach = E(:,end)
h_dach = h - e_dach
hold on
plot(t,h_dach,'*-')
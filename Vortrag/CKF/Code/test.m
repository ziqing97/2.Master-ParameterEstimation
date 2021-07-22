clc 
close all
clear


%% Conditions
n = 2;      %number of state
q = 0.1;    %std of process 
r = 0.1;    %std of measurement
Q = q^2 * eye(n); % covariance of process

f = @(x)[sin(2*x(2));cos(2*x(1))];  % nonlinear state equations
h = @(x)[x(1);x(2)];                               % measurement equation
R = eye(2) * r^2;

s = [0;0];                                % initial state
x = [0.0581;-0.2192]; %initial state          % initial state with noise

P = eye(n)*0.000001;                               % initial state covraiance
P = eye(n);

z = h(s) + r * randn;                     % measurment

% point test
% ckf
S = chol(P);
S = S';
m = 2*n;
xi = zeros(n,m);
p_ckf = zeros(n,m);
for i = 1:n
    xi(i,i) = 1 * sqrt(m/2);
    xi(i,i+n) = -1 * sqrt(m/2);
end
for i = 1:m
    p_ckf(:,i) = S * xi(:,i) + x;
end
% ukf
alpha=1e-3;                                 %default, tunable
L=numel(x); 
ki=0; 
lambda=alpha^2*(L+ki)-L;                    %scaling factor
c=L+lambda;    
c = sqrt(c);
%scaling factor
A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
p_ukf = [x Y+A Y-A]; 

[x2, P2] = ukf(f,x,P,h,z,Q,R); 
[x3, P3] = ckf(f,x,P,h,z,Q,R); 

ckff = trackingCKF(f,h,x,'stateCovariance',P,'ProcessNoise',Q,'MeasurementNoise',R);
[x4,p4] = correct(ckff,z);


figure
hold on
plot(p_ckf(1,:),p_ckf(2,:),'og')
plot(p_ukf(1,:),p_ukf(2,:),'or')
scatter(x(1),x(2))
axis equal

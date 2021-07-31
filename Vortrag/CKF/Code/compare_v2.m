clc
clear
close all
%% comparing scaled unscented transformation vs cubature transformation
% dimension 6 

n = 6;      %number of state
q = 0.1;    %std of process 
r = 0.1;    %std of measurement
Q = q^2 * eye(n); % covariance of process


f = @(x)[sin(2*x(6));cos(2*x(1));sin(2*x(2));cos(2*x(3));sin(2*x(4));cos(2*x(5))];  % nonlinear state equations
f = @(x)[sin(2*x(6));cos(2*x(5));sin(2*x(4));cos(2*x(3));sin(2*x(2));cos(2*x(1))];
% f = @(x)[x(1)/2+1;x(2)/2+1;x(3)/2+1;x(4)/2+1;x(5)/2+1;x(6)/2+1];

h1 = @(x)[x(1);x(2)];                               % measurement equation
R1 = eye(2) * r^2;
 

s = [0;0;0;0;0;0];                                % initial state
x = s + q * randn(n,1); %initial state          % initial state with noise
P = eye(n)*0.1;                               % initial state covraiance
N = 100;  

% total dynamic steps
xV = zeros(n,N);          %estimate       % allocate memory
xV1 = zeros(n,N);          %estimate       % allocate memory
xC = zeros(n,N);
xC1 = zeros(n,N);

sV = zeros(n,N);          % actual
zV = zeros(1,N);
zV1 = zeros(2,N);


% start status
x1 = x;
x2 = x;
x3 = x;
x4 = x;
P1 = P;
P2 = P;
P3 = P;
P4 = P;

for k=1:N
  % simulate the real state
  z1 = h1(s) + r*randn;                     % measurment
  zV1(:,k)  = z1;                             % save measurment
  s = f(s) + q*randn(n,1);                % update process 

  
  % ukf

  [x1, P1] = ukf(f,x1,P1,h1,z1,Q,R1);            % ukf
  xV1(:,k) = x1;                            % save estimate
  
  % ckf
  [x3, P3] = ckf(f,x3,P3,h1,z1,Q,R1); 

  xC1(:,k) = x3;

  
  sV(:,k)= s;                             % save actual state
end
figure
for k=1:n                                 % plot results
  subplot(n,1,k)
  plot(1:N, sV(k,:),'.-b', 'linewidth',2);
  hold on;
%   plot(1:N, xV1(k,:),'.-r', 'linewidth',2) 
  plot(1:N,xV1(k,:), '.-r', 'linewidth',2)
  plot(1:N,xC1(k,:), '.-g', 'linewidth',2)
%   plot(1:N,xC1(k,:), '.-m', 'linewidth',2)
%   legend('actual state','UKF estimated state with measurements of x(1)',...
%       'UKF estimated state with measurements of x(1) and x(2)',...
%       'UKF estimated state with measurements of x(1)',...
%       'UKF estimated state with measurements of x(1) and x(2)');
end

figure
for k=1:n                                 % plot results
  subplot(n,1,k)
  plot(1:N, sV(k,:)-xV1(k,:), 'linewidth',2); hold on
%   plot(1:N, sV(k,:)-xV1(k,:), 'linewidth',2)
  plot(1:N, sV(k,:)-xC1(k,:), 'linewidth',2)
%   plot(1:N, sV(k,:)-xC1(k,:), 'linewidth',2)
end

legend('error, case:  ukf','error, case: ckf')

% figure
%                                  % plot results
% plot(sV(1,:),sV(2,:))
% hold on
% plot(xV1(1,:),xV1(2,:))
% plot(xC1(1,:),xC1(2,:))

% k = [0:N;0:N];
% for i=1:length(k)
% t(i,:)=f(k(:,i));
% end
% plot(t(:,1),t(:,2))
function[x_hat,P_hat] = ckf(f,x,P,h,z,Q,R)
% [x, P] = ckf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           hmeas: fanction handle for h(x)
%           z: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance
%
%% 
n = length(x); % dimension
l = length(z);
m = n * 2;

%% Time Update
S = chol(P);
S = S';

xi = zeros(n,m);
x_cub = zeros(n,m);

% sample points
for i = 1:n
    xi(i,i) = 1 * sqrt(m/2);
    xi(i,i+n) = -1 * sqrt(m/2);
end

% cubature points
for i = 1:m
    x_cub(:,i) = S * xi(:,i) + x;
end

% propogated cubature points
x_cub_pro = zeros(n,m);
for i = 1:m
    x_cub_pro(:,i) = f(x_cub(:,i));
end

% predicted state & covariance
x_nnp = sum(x_cub_pro,2) / m;
P_nnp = zeros(n);
for i = 1:m
    P_nnp = P_nnp + x_cub_pro(:,i) * x_cub_pro(:,i)';
end
P_nnp = P_nnp / m;
P_nnp = P_nnp - x_nnp * x_nnp' + Q;

%% Measurement Update
% Factorize
S_nnp = chol(P_nnp);
S_nnp = S_nnp';

% cubature points
x_nnp_cub = zeros(n,m);
for i = 1:m
    x_nnp_cub(:,i) = S_nnp * xi(:,i) + x_nnp;
end

% propagated cubature points
z_nnp_cub_pro = zeros(l,m);
for i = 1:m
    z_nnp_cub_pro(:,i) = h(x_nnp_cub(:,i));
end

% predicted measurement & innovation covariance matrix
z_nnp = sum(z_nnp_cub_pro,2) / m;

P_zz_nnp = zeros(l);
for i = 1:m
    P_zz_nnp = P_zz_nnp + z_nnp_cub_pro(:,i) * z_nnp_cub_pro(:,i)';
end
P_zz_nnp = P_zz_nnp / m;
P_zz_nnp = P_zz_nnp - z_nnp * z_nnp' + R;

% cross-covariance matrix
P_xz_nnp = zeros(n,l);
for i = 1:m
    P_xz_nnp = P_xz_nnp + x_cub_pro(:,i) * z_nnp_cub_pro(:,i)';
end
P_xz_nnp = P_xz_nnp / m - x_nnp * z_nnp';

% Kalman gain
K = P_xz_nnp / P_zz_nnp;

% updated state
x_hat = x_nnp + K * (z - z_nnp);
P_hat = P_nnp - K * P_zz_nnp * K';

end

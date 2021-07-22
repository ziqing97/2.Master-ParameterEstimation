clc
clear
% close all

%% This program is written for the undertanding of Baysian model regresssion

% May 2020, Mohammad Tourian


%% generating syntetic data
a0 = -0.3;
a1 =  0.5;

sigma_y = (0.2)^2; % the uncertainty of observations

sigma_x = (0.5)^2; % uncertainty of unknown

sigma2p = 3; % a priori sigma_0 

varSigma2p = 0.01; % a priori variance for the given  a priori sigma_0
% varSigma2p = 1;

mu0 = [0 ; 0];

x = [a0; a1]; % unknown vector


N = 100; % length of syntetick sample


t = 2 * rand(N,1) -1 ; % generating random number of the t  between [-1 1]

A = [ones(N,1) t];   % Design matrix

y= A * x + sqrt(sigma_y) * randn(N,1);  % generating synthetic observations with a Gaussian noise with std sqrt(sigma_y)

%% prior p(x|sigma_x^-1) = N (x|0,sigma_x I)


RandomUnknown = mvnrnd(zeros(1,2),sigma_x * eye(2,2),10^6); % 10^6 values are generated for a0 and a1 using the prior

pts = linspace(-1, 1, 50);

[Prior, xe, ye] = histcounts2(RandomUnknown(:,2), RandomUnknown(:,1), pts, pts);


% taking 10 data points for visualisation
random_x = datasample(RandomUnknown,10);  
% random_x = datasample(RandomUnknown,1000); 


% plotting Prior
Fig = figure;

subplot(2,3,2);

imagesc(xe,ye,Prior); pbaspect([1,1,1]);
view(0,-90); % so that x and y access 

set(gca,'fontsize',20)
xlabel('a_0');ylabel('a_1')
title('prior')


subplot(2,3,3);
for i=1:length(random_x)
    
    plot(t,A*random_x(i,:)');
    
    hold on
    
end

set(gca,'fontsize',20)
xlabel('t');ylabel('y')
title('data space')



%% Selecting data from synthetic generated data sets and updating prioir


Nr_sel = input('please input number of data to be selected:   ');


[y_select, id] = datasample(y,Nr_sel); % randomly selecting one data

t_select = t(id); % finding the correspoing t value



%% likelihood
% P(y|t,x,sigma_y^-1) = \product_i=1^{Nr_sel} N (y_i|f(t_i,x), sigma_y) 

Yartificial = -3:.001:3;

for i=1:Nr_sel

P_y(i,:) = pdf('Normal',Yartificial,y_select(i),sqrt(sigma_y));


end


% taking product of all likelihoods
if Nr_sel~=1
    PYgiven_t_and_x = prod(P_y)';
else
    PYgiven_t_and_x = P_y;
end



% producing likelihood function p(y|t,x) for the selected data points as a function of x.

a0pts = linspace(-1, 1, 49);
a1pts = linspace(-1, 1, 49);


% generating likelihoods based on each data as a function of x
for iCountrNrselPoints=1:Nr_sel
    
    likelihood{iCountrNrselPoints} = nan(49,49); 
    
    for iCounterOnX=1:length(a0pts)
        
        a1tmp = (Yartificial - a0pts(iCounterOnX) ) / t_select(iCountrNrselPoints);
        
        likelihood{iCountrNrselPoints}(:,iCounterOnX) = interp1( a1tmp, P_y(iCountrNrselPoints,:), a1pts)';
        
    end
end


% generating likelihood function as the product of all likelihoods 
LikelihoodFunction = ones(49,49);

for iCountrNrselPoints = 1 : Nr_sel
    
LikelihoodFunction = LikelihoodFunction .* likelihood{iCountrNrselPoints};

end


subplot(2,3,4);

imagesc(a0pts,a1pts,LikelihoodFunction); pbaspect([1,1,1]);

view(0,-90)

hold on

plot(a0,a1,'r+')
set(gca,'fontsize',20)
xlabel('a_0');ylabel('a_1')
title('likelihood')

%% searching for posterior


Posterior = Prior.*LikelihoodFunction / sum(sum(Prior.*LikelihoodFunction));

% the maximum of Posterior should be the solution, however for the case of
% normal distribution of likelihood and Prior in the lecture on 06.05.2020
% we have derived the folloing formulae, which give the analytical solution


% we can also
A_select = [ones(Nr_sel,1) t_select]; % Design matrix for the selected data points

%Lambda
sigma_post = inv(A_select' *sigma_y^-1 * eye(Nr_sel,Nr_sel)* A_select + sigma_x^-1);

mu =  sigma_post * A_select' *sigma_y^-1 * eye(Nr_sel,Nr_sel)* y_select; 

sigma2B = (Nr_sel + 2 * sigma2p ^ 2/varSigma2p +2)^-1 *(2 * (sigma2p ^ 2/varSigma2p + 1) * sigma2p + ...
   ( mu-mu0)'  * sigma_x^-1*eye(length(x),length(x)) * ( mu-mu0) + (y_select - A_select * mu)' * sigma_y^-1 * eye(Nr_sel,Nr_sel) * (y_select-A_select*mu));

sigmahat_post = sigma2B * sigma_post;


% Based on the obtained solution, we can now generate Posterior (PosteriorGenerated)
R = mvnrnd(mu,sigma_post,10^6);
pts = linspace(-1, 1, 50);
PosteriorGenerated = histcounts2(R(:,2), R(:,1), pts, pts);


%---------------------

subplot(2,3,5);
imagesc(pts,pts,Posterior); pbaspect([1,1,1]); view(0,-90); hold on
plot(a0,a1,'r+')
set(gca,'fontsize',20)
xlabel('a_0');ylabel('a_1')
title('posterior')

plot(mu(1),mu(2),'m+'); % ploting the x_MAP on top of posterior


% plotting 10 lines from simulated x
random_x = datasample(R,10);

subplot(2,3,6);
for i=1:length(random_x)
    plot(t,A*random_x(i,:)');
    hold on
end
hold on

plot(t_select,y_select,'bo'); % plotting selected data points

set(gca,'fontsize',20)
xlabel('t');ylabel('y')

%%%
%=====comparing Posterior with generated Posterior out of mu and sigma============
figure;

subplot(1,2,1);imagesc(pts, pts, Posterior); view(0,-90);pbaspect([1,1,1]);title('Posterior: Likelihood x Prior')

subplot(1,2,2);imagesc(pts, pts, PosteriorGenerated); view(0,-90); pbaspect([1,1,1]); title('Posterior: based on obtained solution')



%% Comparing with weighted least squares case

P = eye(length(t_select),length(t_select)) * sigma_y^-1; % weight matrix

x_hat = inv(A_select' * P * A_select) * A_select' * P* y_select; % least squares result

Q_xhat = inv(A_select' * P * A_select); % Covariance matrix

e_hat = y_select - A_select * x_hat; % e_hat

sigma02 = e_hat' * P * e_hat / (length(y_select)-2); % Estimation of variance of unit weight

Q_hat_xhat = sigma02 * Q_xhat; % updating covariance matrix



% adding the result of WLS to the figure  to be able to compare it with the
% MAP solution
figure(Fig);

subplot(2,3,1);
Pos = get(gca,'position');
subplot(2,3,5)

plot(x_hat(1), x_hat(2), 'k+')


leg = legend('true values',['MAP, error=' num2str(sqrt((mu(1)-a0)^2+(mu(2)-a1)^2))] ,...
    ['WLS, error= ' num2str(sqrt((x_hat(1)-a0)^2+(x_hat(2)-a1)^2))]);

set(leg,'position',Pos)

subplot(2,3,1);
axis off



%%%
%=====Plotting fitted lines using different approaches============

figure;

plot(t_select,y_select,'bo')

hold on

T = linspace(-1,1,100);

A_T = [ones(length(T),1) T'];

Y_MAP = mu(1) + mu(2)*T;
Q_y_MAP = sigma_y*eye(length(Y_MAP),length(Y_MAP)) + A_T* sigmahat_post* A_T';

uncer_MAP = diag(Q_y_MAP).^.5;

Y_WLS = x_hat(1) + x_hat(2)*T;
Q_y_WLS = (A_T*Q_hat_xhat*A_T');

uncer_WLS = diag(Q_y_WLS).^.5;

Y_real = a0 + a1*T;

plot(T,Y_real,'r','linewidth',2);

plot(T,Y_MAP,'m','linewidth',2);


plot(T,Y_WLS,'k','linewidth',2);



legend('observations','true mmodel',['MAP, RMSE=' num2str(rms(Y_MAP-Y_real))] ,...
    ['WLS, RMSE= ' num2str(rms(Y_WLS-Y_real))])

plot(T,Y_MAP-uncer_MAP','.-m');hold on;plot(T,Y_MAP+uncer_MAP','.-m')

plot(T,Y_WLS-uncer_WLS','.-k');hold on;plot(T,Y_WLS+uncer_WLS','.-k')






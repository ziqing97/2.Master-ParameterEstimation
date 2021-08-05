%% Part 2
clc
clear
close all
load("Discharge_Data.mat")
load("TWSAerr_Data.mat")
load("ETdata_Ue6.mat")
load("PrecData_Ue6.mat")

% unit m^3
t_twsa = TWSA_data.time;
t_dis = Discharge.date;
t_p = Prec_data.time;
t_et = ET_data.time;

TWSA = TWSA_data.value;
Discharge = Discharge.value;
ET = ET_data.value * 30;
Pre = Prec_data.value * 30;

std_ET = ET_data.std * 30;
std_P = Prec_data.std * 30;
std_TWSA = TWSA_data.std;

%% Data preprocessing
for i = 1:length(t_dis)    
    t_dis(i) = datetime(t_dis(i).Year, t_dis(i).Month, t_dis(i).Day);
end

for i = 1:length(t_twsa)
    t_twsa(i) = datetime(t_twsa(i).Year, t_twsa(i).Month, t_twsa(i).Day);
end

R = zeros(204,1);
std_R = zeros(204,1);
i = 1;
for year = 2003:2019
    for month = 1:12
        t_i = datetime(year,month,16);
        idx_R = find(t_dis == t_i);
        if ~isempty(idx_R)
            R(i) = mean(Discharge(idx_R-15:idx_R+15),'omitnan')*3600*24*31;  % m^3/month 
            std_R(i) = R(i)/10;
        else 
            R(i) = NaN;
        end
        t_r(i) = t_i;
        i = i+1;
    end
end

t = t_r';
R = R';

%% 
% load ('parameter.mat')
% tau = parameter(1);
% S0 = parameter(2);
% omega = parameter(3);

tau = 36.03/30; % month^-1
omega = sqrt(2.385e-4)*31; % month^-1
S0 = 1.802e12; % m^3


x0 = [TWSA(1);R(1)];

para1 = omega^2*tau*(exp(1/tau)-1);
para2 = exp(1/tau);

%% Matlab dynamic system
% A = [0,-1;omega^2,-1/tau];
% B = [0,1,-1;omega^2,0,0];
%  
% C = [1,0;0,1];
% D = [1,0,0;0,0,0];
% 
% sys = ss(A,B,C,D,1,'StateName',{'TWSA','Discharge'},'InputName',{'S0','P','ET'});
% sys.TimeUnit = 'month';


%% self writing functions
code = 1;
if code == 1
    % case 1
    A = [0,-1;omega^2,-1/tau];
    B = [0,1,-1;omega^2,0,0];
    C = eye(2);
    G = expm([A,B;zeros(3,5)]);
    Ad = G(1:2,1:2);
    Bd = G(1:2,3:5);
    A = Ad;
    B = Bd;
else 
    % case 2
    A = [1,-1;para1,para2];
    B = [0,1,-1;para1,0,0];
    C = eye(2);
end

P = cell(204,1);
Pnnp = cell(204,1);
P{1} = [1,0;0,1];
Pnnp{1} = [1,0;0,1];
f = @(x,u)[A*x+B*u];
h1 = @(x,u)[C*x];
h2 = @(x,u)[x(2)];

x = zeros(2,204);
xnnp = zeros(2,204);

x(:,1) = [TWSA(1);R(1)];
xnnp(:,1) = [TWSA(1);R(1)];

std_S0 = S0/10;

for i = 1:203
    if isnan(R(i+1))
        i
        u = [S0;Pre(i);ET(i)];
        P{i+1} = P{i};
        x(:,i+1) = f(x(:,i),u);
        xnnp(:,i+1) = x(:,i+1);
        Pnnp{i+1} = P{i};
    else
        u = [S0;Pre(i);ET(i)];
        check = t(i);
        index = find(check == t_twsa);
        if isempty(index)
            z = R(i+1);
            R_pp = std_R(i)^2;
            Q_pp = (B * diag([std_S0;std_P(i);std_ET(i)].^2) * B');
            [xnnp(:,i+1),Pnnp{i+1},x(:,i+1), P{i+1}] = ukf(f,x(:,i),u,P{i},h2,z,Q_pp,R_pp);
        else
            z = [TWSA(index+1);R(i+1)];
            R_pp = [std_TWSA(index)^2,0;0,std_R(i)^2];
            Q_pp = (B * diag([std_S0;std_P(i);std_ET(i)].^2) * B');
            [xnnp(:,i+1),Pnnp{i+1},x(:,i+1), P{i+1}] = ukf(f,x(:,i),u,P{i},h1,z,Q_pp,R_pp);
        end
    end
end
x = x';
xnnp = xnnp';

std2_TWSA = zeros(203,1);
std2_R = zeros(203,1);
std_xnnp_TWSA =zeros(203,1);
std_xnnp_R = zeros(203,1);
for i = 1:203
    std2_TWSA(i) = sqrt(P{i+1}(1,1));
    std2_R(i) = sqrt(P{i+1}(2,2));
    
    std_xnnp_TWSA(i) = sqrt(Pnnp{i+1}(1,1));
    std_xnnp_R(i) = sqrt(Pnnp{i+1}(2,2));
end

t1 = t;
t2 = t_twsa;
t = datenum(t);
t_twsa = datenum(t_twsa);


figure
hold on
title('TWSA')

err_plot(t(2:end),x(2:end,1),std2_TWSA,"blue")
err_plot(t(1:end-1),xnnp(1:end-1,1),std_xnnp_TWSA,"green")
err_plot(t_twsa(1:147),TWSA(1:147),std_TWSA(1:147),"red")
err_plot(t_twsa(148:end),TWSA(148:end),std_TWSA(148:end),"red")

datetick('x')
legend('Estimate','','Prior','','Likelihood GRACE','','Likelihood GRACE-FO','')
set(gca,'FontSize', 20);
pbaspect([3 1 1])

figure
hold on
title('Discharge')
R = R';
t_r = datenum(t_r);
t_r = t_r';
err_plot(t(2:end),x(2:end,2),std2_R,"blue")
err_plot(t(1:end-1),xnnp(1:end-1,2),std_xnnp_R,"green")
err_plot(t_r,R,std_R,"red")

datetick('x')
legend('Estimate','','Prior','','Likelihood','')
set(gca,'FontSize', 20);
pbaspect([3 1 1])
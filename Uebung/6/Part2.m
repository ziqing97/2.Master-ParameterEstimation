%% Part 2
clc
clear
close all
load("Discharge_Data.mat")
load("TWSA_Data.mat")
load("ETdata_Ue6.mat")
load("PrecData_Ue6.mat")

% unit m^3
t_twsa = TWSA.date;
t_dis = Discharge.date;
t_p = Prec_data.time;
t_et = ET_data.time;

TWSA = TWSA.value;
Discharge = Discharge.value;
ET = ET_data.value * 30;
Pre = Prec_data.value * 30;

std_ET = ET_data.std * 30;
std_P = Prec_data.std * 30;

%% Data preprocessing
for i = 1:length(t_dis)    
    t_dis(i) = datetime(t_dis(i).Year, t_dis(i).Month, t_dis(i).Day);
end

for i = 1:length(t_twsa)
    t_twsa(i) = datetime(t_twsa(i).Year, t_twsa(i).Month, t_twsa(i).Day);
end

R = zeros(204,1);
i = 1;
for year = 2003:2019
    for month = 1:12
        t_i = datetime(year,month,16);
        idx_R = find(t_dis == t_i);
        if ~isempty(idx_R)
            R(i) = mean(Discharge(idx_R-15:idx_R+15),'omitnan')*3600*24*30;  % m^3/month 
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
load ('parameter.mat')
tau = parameter(1);
S0 = parameter(2);
omega = parameter(3);

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
A = [1,-1;para1,para2];
B = [0,1,-1;para1,0,0];
C = eye(2);

P = [1,0;0,1];
f = @(x,u)[A*x+B*u];
h1 = @(x,u)[C*x];
h2 = @(x,u)[x(2)];

x = zeros(2,204);
x(:,1) = [TWSA(1);R(1)];

Q_pp = eye(2);
R_pp = eye(2);
R_pp2 = 1;

for i = 1:204
    if isnan(R(i))
        i
    else
        u = [S0;Pre(i);ET(i)];
        check = t(i);
        index = find(check == t_twsa);
        if isempty(index)
            z = R(i+1);
            [x(:,i+1), P] = ukf(f,x(:,i),u,P,h2,z,Q_pp,R_pp2);
        else
            z = [TWSA(index+1);R(i+1)];
            [x(:,i+1), P] = ukf(f,x(:,i),u,P,h1,z,Q_pp,R_pp);
        end
    end
end
figure
hold on
title('TWSA')
plot(t(1:i),x(1,1:i))
plot(t_twsa,TWSA)


figure
hold on
title('Discharge')
plot(t(1:i),x(2,1:i))
plot(t(1:i),R(1:i))


% for i = 1:147
%     u = [S0;Pre(i);ET(i)];
%     z = [TWSA(i+1);R(i+1)];
%     [x(:,i+1), P] = ukf(f,x(:,i),u,P,h1,z,Q_pp,R_pp);
% end
% figure
% hold on
% plot(t,x(1,:))
% plot(t_twsa,TWSA)
% 
% for i = 148:190
%     u = [S0;Pre(i);ET(i)];
%     z = R(i+1);
%     [x(:,i+1), P] = ukf(f,x(:,i),u,P,h2,z,Q_pp,R_pp2);
% end
% 
% for i = 191:203
%     u = [S0;Pre(i);ET(i)];
%     z = [TWSA(i-43);R(i+1)];
%     [x(:,i+1), P] = ukf(f,x(:,i),u,P,h1,z,Q_pp,R_pp);
% end
clc
clear
close all
load("Discharge_Data.mat")
load("TWSAerr_Data.mat")
load("ETdata_Ue6.mat")
load("PrecData_Ue6.mat")

% unit m^3
t1 = TWSA_data.time;
t2 = Discharge.date;
t3 = Prec_data.time;
t4 = ET_data.time;

TWSA = TWSA_data.value;
Discharge = Discharge.value*24*3600*30;
ET = ET_data.value*30;
P = Prec_data.value*30;
%% Part 1
% Data preprocessing
for i = 1:length(t2)    
    t2(i) = datetime(t2(i).Year, t2(i).Month, t2(i).Day);
end

Dis = zeros(length(t1),1);

for i = 1:length(t1)
    t1(i) = datetime(t1(i).Year, t1(i).Month, t1(i).Day);
    idx = find(t2 == t1(i));
    if isempty(idx)
        Dis(i) = NaN;
    else
        Dis(i) = mean(Discharge(idx-15:idx+14),'omitnan');
    end
end

idx_t3 = t3(t3>=t1(1) & t3<=t1(end));

t = t1;
clear t1 t2 Discharge 

% not using GRACE-FO
R = Dis(1:147);
dS = TWSA(1:147);
t = t(1:147);

% Least Square(or TLS ?)
% Design matrix
A = [R(1:end-1), dS(1:end-1), ones(length(R)-1,1)];
para = (A' * A) \ A' * R(2:end);

% parameter
tau = 1 / log(para(1));
S0 = para(3) / para(2);
omega = sqrt(para(2) / tau / (exp(1/tau) -1));

parameter = [tau,S0,omega];
save('parameter.mat','parameter')
save('sigma_S0.mat','sigma_S0')
clc
clear
close all
load("Discharge_Data.mat")
load("TWSA_Data.mat")


% unit m^3
t1 = TWSA.date;
t2 = Discharge.date; 

TWSA = TWSA.value;
Discharge = Discharge.value;

TWSA2 = TWSA_data.value;
t3 = TWSA_data.time;

figure
plot(t3,TWSA2)
hold on 
plot(t1,TWSA)

% plot(t1, TWSA)
% figure
% plot(t2,Discharge)

% pre
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
t = t1;
% not using GRACE-FO
Dis = Dis(1:147);
TWSA = TWSA(1:147);
t = t(1:147);

figure
subplot(2,1,1)
plot(t,Dis)
title('Discharge')
xlabel('time')
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)

subplot(2,1,2)
plot(t,TWSA)
title('TWSA')
xlabel('time')
ylabel('')
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)

figure
plot(Dis,TWSA)
xlabel('Discharge')
ylabel('\Delta S')
set(gca,'YGrid','on')
set(gcf,'color','w')
set(gca,'fontsize',20)


k = length(Dis);
Dis = Dis(1:k);
TWSA = TWSA(1:k);

R = Dis(2:end);
dS = TWSA(2:end);

% A = [Dis(1:end-1), ones(length(R),1), dS];
% para = (A' * A) \ A' * R;
% y_n = A * para;
% alpha = para(1);
% beta = para(2);
% gamma = para(3);

A = [R,ones(k-1,1)];
para = (A'*A)\A'*dS;
dS_n = A * para;


figure
plot(dS)
hold on
plot(dS_n)

% S0 = beta/gamma;
% omega = sqrt(gamma / alpha);
% tau = alpha / (1-alpha);


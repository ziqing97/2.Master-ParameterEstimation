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
t = t1;
clear t1 t2 Discharge 

% not using GRACE-FO
R = Dis(1:147);
dS = TWSA(1:147);
t = t(1:147);

% delta t = 1 month
% desigh matrix
A = [R(1:end-1), dS(1:end-1) + dS(2:end), ones(length(t)-1,1)];
para = (A' * A) \ A' * R(2:end);
y = A * para;

figure
hold on
plot(R(2:end))
plot(y)

S0 = para(3)/para(2)/2;
omega = sqrt(para(2)) / para(1);

%% Part 2
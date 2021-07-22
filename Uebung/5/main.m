clc
close all
clear

load('sim_TWS_Greenland.mat')
t = data_real(:,1); 
s = data_real(:,2);
len = length(t);

% c
A = [t, ones(len,1), cos(2 * pi * t / 12), sin(2 * pi * t / 12)];
x_dach = (A'* A) \ A'* s;
s_dach = A * x_dach;

e_dach = s - A * x_dach;
sigma = e_dach' * e_dach / (len - length(x_dach));
Sigma_xx = sigma * inv(A'* A);
Sigma_yy = A * Sigma_xx * A';


t_new = [46;47;115;116;117;118];
len2 = length(t_new);
A_new = [t_new, ones(len2,1), cos(2 * pi * t_new / 12), sin(2 * pi * t_new / 12)];
s_new = A_new * x_dach;
Sigma_ynyn = A_new * Sigma_xx * A_new';

figure
hold on
scatter(t,s,'o')
scatter(t,s_dach,'+')
scatter(t_new,s_new,'*')
xlabel('time')
ylabel('TWSC')
set(gca,'FontSize', 20);
legend('Observation','Adjusted Observation','Filled Gap')
saveas(gca,'gapfilling.png')

%% 
Sigma_yy2 = diag(diag(Sigma_yy));
Sigma_xx2 = diag(diag(Sigma_xx));
D = inv(A' * inv(Sigma_yy2) * A + inv(Sigma_xx2));
mu = D * (A' * inv(Sigma_yy2)  * s + inv(Sigma_xx2) * x_dach);
std_xx = diag(sqrt(D));
pts = cell(4,1);
Prior = cell(4,1);

figure
lab = {'a','b','C','S'};
for i = 1:4
    pts{i} = linspace(mu(i) - 3 * std_xx(i),mu(i) + 3 * std_xx(i),999);
    Prior{i} = pdf('Normal',pts{i},mu(i),std_xx(i));
    subplot(2,2,i)
    plot(pts{i},Prior{i})
    title(lab{i})
    grid on
    set(gca,'FontSize', 20);
end
set(gca,'FontSize', 20);
saveas(gca,'paradist.png')


pts_s = cell(6,1);
Prob_yn = cell(6,1);
Prob_yn_x = cell(6,1);
std_ynyn = sqrt(diag(Sigma_ynyn));
figure
for i = 1:6
    pts_s{i} = linspace(s_new(i) - 3 * std_ynyn(i), s_new(i) + 3 * std_ynyn(i), 999);
    Prob_yn{i} = pdf('Normal',pts_s{i},s_new(i),std_ynyn(i));
    Prob_yn_x{i} = Prob_yn{i} .* Prior{1} + Prob_yn{i} .* Prior{2} + Prob_yn{i} .* Prior{3} + Prob_yn{i} .* Prior{4};
    subplot(3,2,i)
    plot(pts_s{i},Prob_yn_x{i})
    title(t_new(i))
    set(gca,'FontSize', 20);
    grid on
end
saveas(gca,'paradist2.png')


%% 
A2 = [t.^2, t, ones(len,1), cos(2 * pi * t / 12), sin(2 * pi * t / 12)];
x_dach = (A'* A) \ A'* s;
s_dach = A * x_dach;

e_dach = s - A * x_dach;
sigma = e_dach' * e_dach / (len - length(x_dach));
Sigma_xx = sigma * inv(A'* A);
Sigma_yy = A * Sigma_xx * A';


t_new = [46;47;115;116;117;118];
len2 = length(t_new);
A_new = [t_new, ones(len2,1), cos(2 * pi * t_new / 12), sin(2 * pi * t_new / 12)];
s_new = A_new * x_dach;
Sigma_ynyn = A_new * Sigma_xx * A_new';

figure
hold on
scatter(t,s,'o')
scatter(t,s_dach,'+')
scatter(t_new,s_new,'*')
xlabel('time')
ylabel('TWSC')
set(gca,'FontSize', 20);
legend('Observation','Adjusted Observation','Filled Gap')
saveas(gca,'gapfilling.png')

%% 
Sigma_yy2 = diag(diag(Sigma_yy));
Sigma_xx2 = diag(diag(Sigma_xx));
D = inv(A' * inv(Sigma_yy2) * A + inv(Sigma_xx2));
mu = D * (A' * inv(Sigma_yy2)  * s + inv(Sigma_xx2) * x_dach);
std_xx = diag(sqrt(D));
pts = cell(4,1);
Prior = cell(4,1);

figure
lab = {'a','b','C','S'};
for i = 1:4
    pts{i} = linspace(mu(i) - 3 * std_xx(i),mu(i) + 3 * std_xx(i),999);
    Prior{i} = pdf('Normal',pts{i},mu(i),std_xx(i));
    subplot(2,2,i)
    plot(pts{i},Prior{i})
    title(lab{i})
    grid on
    set(gca,'FontSize', 20);
end
set(gca,'FontSize', 20);



pts_s = cell(6,1);
Prob_yn2 = cell(6,1);
Prob_yn_x2 = cell(6,1);
std_ynyn2 = sqrt(diag(Sigma_ynyn));
figure
for i = 1:6
    pts_s{i} = linspace(s_new(i) - 3 * std_ynyn2(i), s_new(i) + 3 * std_ynyn2(i), 999);
    Prob_yn2{i} = pdf('Normal',pts_s{i},s_new(i),std_ynyn2(i));
    Prob_yn_x2{i} = Prob_yn2{i} .* Prior{1} + Prob_yn2{i} .* Prior{2} + Prob_yn2{i} .* Prior{3} + Prob_yn2{i} .* Prior{4};
    subplot(3,2,i)
    plot(pts_s{i},Prob_yn_x2{i})
    title(t_new(i))
    set(gca,'FontSize', 20);
    grid on
end

% vergleich = cell(1,6);
% for i = 1:6
%     vergleich{i} = Prob_yn_x2{i} ./ Prob_yn_x{i};
%     plot()
% end
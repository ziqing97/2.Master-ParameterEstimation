clc
clearvars
% close all

load("cointoss.mat")
count_head = length(find(data == 1));
count_tail = (length(data) - count_head);

H = 0:0.01:1;

% either '1' or other
pdf_choice = '2';

if pdf_choice == '1'
    pdfH = 1;
else
    pdfH = H.^3;
end

figure
hold on
for N = 1:4
    toss = data(1:N);
    toss_head = find(toss == 1);
    R = length(toss_head);
    L = nchoosek(N,R) .* H.^R .* (1-H).^(N-R);
    prob_H_results = L .* pdfH;
    plot(H,prob_H_results)
end

N = 200;
toss = data(1:N);
toss_head = find(toss == 1);
R = length(toss_head);
L = nchoosek(N,R) .* H.^R .* (1-H).^(N-R);
prob_H_results = L .* pdfH;
plot(H,prob_H_results)

legend('1','2','3','4','200')
hold off

%%
pdf1 = 1;
pdf2 = H.^3;
N = 200;
toss = data(1:N);
toss_head = find(toss == 1);
R = length(toss_head);
L = nchoosek(N,R) .* H.^R .* (1-H).^(N-R);
prob_H_results1 = L .* pdf1;
prob_H_results2 = L .* pdf2;
figure
hold on
plot(H,prob_H_results1)
plot(H,prob_H_results2)
legend('case 1','case 2')

clear all;
load ("a_1.mat");load("a_2.mat");load("a_3.mat");
% 3 oscillators
a=plot(a_1(1,:),a_1(2,:),'r') %load(a_1.mat)
yline(1/2,'--r','LineWidth',2)
%%4 oscillators
hold on;
b=plot(a_2(1,:),a_2(2,:),'b') %load(a_2.mat)
yline(1/3,'--b','LineWidth',2)
%%5 oscillators
hold on;
c=plot(a_3(1,:),a_3(2,:),'g') %load(a_3.mat)
yline(1/4,'--g','LineWidth',2)
legend([a,b,c],'3 oscillators','4 oscillators','5 oscillators')
xlabel("k")
ylabel("Inf -log(\rho)/log(\lambda)")
title("stabilizability vs k")
ylim([0 0.6])

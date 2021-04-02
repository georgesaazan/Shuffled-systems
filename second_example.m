%this provides the results of the second example of the paper,
%this requires the JSR toolbox (https://www.mathworks.com/matlabcentral/fileexchange/33202-the-jsr-toolbox),
%this requires a Matlab toolbox for solving optimization problems: Yalmip (https://yalmip.github.io/),
%this requires an semidefinite programming solver: SeDuMi(https://github.com/SQLP/SeDuMi),
%this requires the Econometrics Toolbox in Matlab.(https://fr.mathworks.com/products/econometrics.html)
%% In the first part we give the design figure introduced in the example
%remark 1: the simulation may pause multiple times with "Press any key to
%proceed", to override this, comment out a "pause" in the JSR toolbox:
%JSR_louvain/Methods/jsr.m line 357.
%remark 2: the simulation may take a long time, however the results are saved
%in gamma_result.m with the plots.
clear all;
C = {'r','b','g'}; %define colors used in the figure
for n=3:1:5 %for 3,4 and 5 oscillators
resultx=[]; %initialize the result every iteration on the number of oscillators
for gamma=0:0.01:0.9 %we change the oscillators gain \gamma from 0 to 0.9
    if(gamma<=0.8) %we verify that for \gamma<=0.8 the jsr is 1.02 for all n
        rho=1.02;
    else
    t=jsr(oscillators(n,gamma)); % for each \gamma and n, we generate a set of matrices using oscillators.m
rho=t.bounds(1); %if \gamma>0.8 we calculate the jsr using the jsr toolbox
    end
b=bounds_lambda_2(oscillators(n,gamma),rho,10^-4) %upper bounds on \lambda 
x=-log(rho)/log(b); %we verify numerically that the infimum is reached at the JSR when x<1/(n-1)
if(x<1/(n-1)&&b~=1)
else x=1/(n-1); %we verify that when x>=1/(n-1) the infimum is 1/(n-1).
end 
resultx=[resultx ,x];
end
a(n-2)=plot(0:0.01:0.9,resultx,'color',C{n-2});
yline(1/(n-1),'--','LineWidth',2,'color',C{n-2});
xlabel("k");
ylabel("Inf -log(\rho)/log(\lambda)");
title("stabilizability vs k");
ylim([0 0.6]);
hold on;
end
legend([a(1),a(2),a(3)],'3 oscillators','4 oscillators','5 oscillators');

%% In the second part we give 2 simulations for the example
N=500; %number of steps
o=oscillators(3,0.4); 
F1=o{1};
F2=o{2};
for p12=[1/10,1/70]
x=[-1.5;-0.5;2;-1]; %initial state
traj=x; %traj=[x11 x12;x21 x22]
p21=p12; %probability of going from state 1 to state 2 and vice versa
P=[1-p12,p12;p21,1-p21];
mc=dtmc(P); %define Markov chain with transition matrix P
sigma(1:N)=simulate(mc,N-1); %generate N random variables from Markov chain
q=1/p12;
ka=kappa(sigma,2); %this counts the number of shuffles using kappa.m (Definition 2.1)
for i=1:N
    if sigma(i)==1
       
        x=F1*x;
        traj=[traj x];
    else
       
        x=F2*x;
        traj=[traj x];
    end
end
figure();
subplot(4,1,1);
plot(0:length(ka)-1,ka./(0:length(ka)-1))
xlabel('t');
ylabel('\kappa^{\theta}(t)/t');
subplot(4,1,3);
plot([0:N-1],traj(1,1:N),[0:N-1],traj(2,1:N),'--')
xlabel('t')
ylabel('x_1_1(t), x_1_2(t)')
legend('x_1_1(t)','x_1_2(t)')
%ylim([-3,2])
subplot(4,1,4);
plot([0:N-1],traj(3,1:N),[0:N-1],traj(4,1:N),'--')
xlabel('t')
ylabel('x_2_1(t), x_2_2(t)')
legend('x_2_1(t)','x_2_2(t)')
%ylim([-3,2])
subplot(4,1,2);
stairs([0:N-1],sigma)
axis([0 N 0.9 2.1])
xlabel('t')
ylabel('\theta(t)')
end



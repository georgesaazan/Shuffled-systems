%this provides the results of the first example of the paper,
%this requires the JSR toolbox (https://www.mathworks.com/matlabcentral/fileexchange/33202-the-jsr-toolbox),
%this requires a Matlab toolbox for solving optimization problems: Yalmip (https://yalmip.github.io/),
%this requires an semidefinite programming solver: SeDuMi (https://github.com/SQLP/SeDuMi).
%% In the first part we give the expression of \lambda(\mathcal A,\rho)
clear all;
p=[];
lambda1=0.9;
lambda2=0.2;
phi1=pi/6;
phi2=pi/3;
A1=[1 0 0; 0 lambda1*cos(phi1) -lambda1*sin(phi1); 0 lambda1*sin(phi1) lambda1*cos(phi1)];
A2=[lambda2*cos(phi2) -lambda2*sin(phi2) 0; lambda2*sin(phi2) lambda2*cos(phi2) 0; 0 0 1];
A={A1,A2}
m=length(A);
for i=1:m
A{i}=cell2mat((A(i)));
end
n=size(A{1},1);
P=perms(1:m);
[x,y]=size(P);
B={};
for i=1:x
B{i}= eye(n);
end

for i=1:x
 for j=1:y
B{i}=B{i}*A{P(i,j)};%construct B=\mathcal N_{I}={A1*A2,A2*A1}
 end
end
a=jsr(B); %This gives tight lower and upper bounds on \rho(A1*A2,A2*A1).
Q=R_norm(B,a(2));% This checks if there exsist a solution to the LMIs (28) and (29) using R_norm.m, and returns the solution matrix P in the positive case.
if isempty(Q)==0 %when norm * exists (Q exists).
for i=1:m
p=[p,norm(Q^(0.5)*A{i}*Q^(-0.5))]; %||A_i||_*=||Q^0.5*A_i*Q^-0.5||
end
R=max(p); % R=max||Q^0.5*A_i*Q^-0.5||=1.3278
end
a,Q,R,%a=[\underline{\rho},overline{\rho}],Q is the solution of the LMIs, in the interval [R,+\infty] we have the expression of \lambda{\mathcal A,\rho}.
%Now for \rho \in ]1,R[, we don't have the right expression but we can have %approximations on \lambda(\mathcal A,\rho):



%% In the second part we give approximations on \lambda(\mathcal A,\rho) using 2 Lyapunov methods
rho0=1;
rho1=1.35; %approximate \lambda(\mathcal A,\rho) for \rho \in ]rho0,rho1]
pas=0.05; %pitch for \rho
lambda_1=0; %initiate \lambda(\mathcal A,\rho) for 1st Lyapunov method
lambda_2=0; %initiate \lambda(\mathcal A,\rho) for 2nd Lyapunov method
resultl1=[5]; %initiate the vector values of lambda_1
resultl2=[5]; %initiate the vector values of lambda_2, the 5 will be removed later.
for rho=rho0:pas:rho1 %interval of \rho
lambda_1=bounds_lambda_1({A1,A2},rho,10^-4); %find upper bounds using 1st Lyapunov method (bounds_lambda_1.m)
lambda_2=bounds_lambda_2({A1,A2},rho,10^-4); %find upper bounds using 2nd Lyapunov method (bounds_lambda_2.m)
resultl1=[resultl1,min(lambda_1,((rho-pas)/rho)^2*resultl1(end))];%\lambda is the minimum of the one obtained from the LMIs and the one from the inequality (4)
resultl2=[resultl2,min(lambda_2,((rho-pas)/rho)^2*resultl2(end))];%\lambda is the minimum of the one obtained from the LMIs and the one from the inequality (4)
end
resultl1 = resultl1(resultl1~=5);
resultl2 = resultl2(resultl2~=5);
resultl=[resultl1;resultl2];
figure();
plot(rho0:pas:rho1,resultl(1,:),'-.','color','b') %upper bounds on lambda_i(\mathcal A,\rho)
hold on;
plot(rho0:pas:rho1,resultl(2,:),'color','b')
hold on;
xlim([rho0,rho1]);
hold on
plot(rho0:pas:rho1,a(1)./(rho0:pas:rho1).^2,'--','color','r') %lower bounds on lambda_i(\mathcal A,\rho) from Proposition 5.2
ylim([0.4 1]);
xlim([rho0,rho1]);
xlabel('\rho')
ylabel('\lambda(A,\rho)')
legend('upper bound (Theorem 5.7)','upper bound (Theorem 5.8)','lower bound (Proposition 5.2)')
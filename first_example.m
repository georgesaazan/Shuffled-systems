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
t=rho0:pas:rho1;
for rho=t %interval of \rho
lambda_1=bounds_lambda_1({A1,A2},rho,10^-4); %find upper bounds using 1st Lyapunov method (bounds_lambda_1.m)
lambda_2=bounds_lambda_2({A1,A2},rho,10^-4); %find upper bounds using 2nd Lyapunov method (bounds_lambda_2.m)
resultl1=[resultl1,min(lambda_1,((rho-pas)/rho)^2*resultl1(end))];%\lambda is the minimum of the one obtained from the LMIs and the one from the inequality (5)
resultl2=[resultl2,min(lambda_2,((rho-pas)/rho)^2*resultl2(end))];%\lambda is the minimum of the one obtained from the LMIs and the one from the inequality (5)
end
resultl1 = resultl1(resultl1~=5);
resultl2 = resultl2(resultl2~=5);
y1=zeros(length(t)-1,1); %upper bounds using Theorem 3
y2=zeros(length(t)-1,1); %upper bounds using Theorem 4
y1(1)=resultl1(1); 
y2(1)=resultl2(1); 
for i=1:1:length(t)-2
y1(i+1)=min(y1(i),t(i+1)^2*resultl1(i+1));
y2(i+1)=min(y2(i),t(i+1)^2*resultl2(i+1)); %applying Equation (5) written in the example: min...
end
syms z
y=piecewise(z>=1 &z<1.05,y1(1)/z^2,z>=1.05 &z<1.1, y1(2)/z^2,z>=1.1 &z<1.15,y1(3)/z^2 , z>=1.15 &z<1.2, y1(4)/z^2, z>=1.2 &z<1.25,y1(5)/z^2 , z>=1.25 &z<1.3,y1(6)/z^2 , z>=1.3 &z<1.35,y1(7)/z^2);
figure();
fplot(y,'-.','color','b')%upper bounds on lambda_i(\mathcal A,\rho) using 1st method
hold on;
y=piecewise(z>=1 &z<1.05,y2(1)/z^2,z>=1.05 &z<1.1, y2(2)/z^2,z>=1.1 &z<1.15,y2(3)/z^2 , z>=1.15 &z<1.2, y2(4)/z^2, z>=1.2 &z<1.25,y2(5)/z^2 , z>=1.25 &z<1.3,y2(6)/z^2 , z>=1.3 &z<1.35,y2(7)/z^2);
fplot(y,'color','b')%upper bounds on lambda_i(\mathcal A,\rho) using 2nd method
hold on;
y=a(1)/z^2;
fplot(y,'--','color','r')%lower bounds on lambda_i(\mathcal A,\rho) from Proposition 5
xlim([rho0,rho1])
ylim([0.4,1])
xlabel('\rho')
ylabel('\lambda(A,\rho)')
legend('upper bound (Theorem 3)','upper bound (Theorem 4)','lower bound (Proposition 5)')

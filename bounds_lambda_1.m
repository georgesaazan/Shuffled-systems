%this function outputs bounds on the upper bounds of lambda(rho) using the 1st Lyapunov method (Theorem 3),
%it takes the following arguments respectively: a set of matrices (ceil), rho, the fixed error for the dichotomy.
%ex: bounds_lambda_1({eye(2),ones(2,2)},1,10^-3).
%this requires a Matlab toolbox for solving optimization problems: Yalmip (https://yalmip.github.io/),
%this requires an semidefinite programming solver: SeDuMi (https://github.com/SQLP/SeDuMi).
function c=bounds_lambda_1(B,rhoo,eps)
m=length(B);
for i=1:m
A{i}=cell2mat((B(i)));
end
s=size(A{1});n=s(1);
a=0;
b=1;
err=b-a;
while(err>eps)%dichotomy algorithm
lambda=(a+b)/2;
F=[];
for i=1:m
M{i}=sdpvar(n,n);
F=[F M{i}>=eye(n)]; %step to construct the LMIs
end
for i=1:m
    for j=1:m
        if(j~=i)
            F=[F,A{j}'*M{i}*A{j} <=rhoo^2*M{i}];%step to construct the LMIs
        end
    end
end
for i=1:m-1
    F=[F A{i}'*M{i+1}*A{i} <=rhoo^2*M{i}];%step to construct the LMIs
    
end
F=[F,A{m}'*M{1}*A{m} <=rhoo^2*lambda^(2*m)*M{m}];%step to construct the LMIs

ops = sdpsettings('solver','sedumi');
if m==1 %case of 1 matrix
    su=trace(M{1});
else su=sum(trace(M{i}));
end
optimize(F,su,ops);%solve the LMIs
[primalfeas,dualfeas] = check(F);
check(F);
if(min(primalfeas)<-1e-7)
   a=lambda
   err=abs(lambda-b);
else %feasible LMIs
    b=lambda
    err=abs(lambda-a);
end 
end
c=b;%output bounds


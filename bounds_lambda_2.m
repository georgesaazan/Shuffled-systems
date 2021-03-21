%this function outputs bounds on the upper bounds of lambda(rho) using the 2nd Lyapunov method (Theorem 5.8),
%it takes the following arguments respectively: a set of matrices (ceil), rho, the fixed error for the dichotomy.
%ex: bounds_lambda_2({2*eye(2),ones(2,2)},1,10^-3).
%this requires a Matlab toolbox for solving optimization problems: Yalmip (https://yalmip.github.io/),
%this requires an semidefinite programming solver: SeDuMi (https://github.com/SQLP/SeDuMi).
function c=bounds_lambda_2(B,rhoo,eps)
m=length(B);
C=PowerSet(1:m)%powerset of I (PowerSet.m)
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
for i=1:2^m-1
M{i}=sdpvar(n,n);
F=[F M{i}>=eye(n)];  %step to construct the LMIs
end
if(m~=1)
for i=1:m
    for j=1:2^m-1
        if(length(union([i],cell2mat(C(j))))~=m)%check if J union {i} ~= I 
            l=find(cellfun(@(x) isequal(x,reshape(union([i],cell2mat(C(j))),1, [])), C));%find the index of J union {i} in the powerset
            if isempty(l)~=1
            F=[F,A{i}'*M{l}*A{i} <=rhoo^2*M{j}];%step to construct the LMIs
            end
        else
            F=[F,A{i}'*M{1}*A{i} <=rhoo^2*lambda^2*M{j}];%step to construct the LMIs
        end
    end
end
else 
    M{2}=sdpvar(n,n);
F=[F M{2}>=eye(n)]; 
    F=[F,A{1}'*M{1}*A{1} <=rhoo^2*lambda^2*M{1}];
     F=[F,A{1}'*M{1}*A{1} <=rhoo^2*lambda^2*M{2}];
    
end
ops = sdpsettings('solver','sedumi');
if m==1 % if one matrix is present in the set 
    su=trace(M{1})+trace(M{2});
else su=sum(trace(M{i}));
end
optimize(F,su,ops);%solve the LMIs
[primalfeas,dualfeas] = check(F);
check(F);
if(min(primalfeas)<-1e-7)%1e-5
   a=lambda
   err=abs(lambda-b);
else %%feasable lmi
    b=lambda
    err=abs(lambda-a);
end 
end
c=b;%output bounds



%this function takes 2 arguments: a set of matrices (ceil) and a scalar,
%this function checks if the norm * exists for the Theorem 5.6 (2) and if it exists it outputs the corresponding Q by solving the LMIs else it gives the empty set.
%this requires a Matlab toolbox for solving optimization problems: Yalmip (https://yalmip.github.io/),
%this requires an semidefinite programming solver: SeDuMi (https://github.com/SQLP/SeDuMi).
%ex: R_norm({2*eye(2),ones(2,2)},1).
function a=R_norm(B,rho)
n=length(B);
Q=sdpvar(max(size(B{1})),max(size(B{1})));
F=[];
p=[];
for i=1:n
F=[F, B{i}'*Q*B{i}<=rho^2*Q,trace(Q)==1];  %%check if exists a quadratic norm satisfying ||M||_{*}<=\rho(N_I)
end
F=[F ,Q>=0,trace(Q)==1];%the LMIs  
ops = sdpsettings('solver','sedumi');
optimize(F,trace(Q));

[primalfeas,dualfeas] = check(F);

if(min(primalfeas)>-1e-3) %%if the LMIs are feasible
    a=value(Q);
else
    a=[];
end
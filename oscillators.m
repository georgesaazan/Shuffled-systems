%this function outputs the n-1 state matrices for n oscillators
function A=oscillators(n,gamma)
phi=pi/6;
R=1.02*[cos(phi) -sin(phi);sin(phi) cos(phi)];
K=gamma*eye(2,2);
B=[];
for(i=1:n-1)
A=kron(eye(n-1),R);
A(2*i-1:2*i,2*i-1:2*i)=A(2*i-1:2*i,2*i-1:2*i)-2*K;
if(i~=n-1)
A(2*i+1:2*i+2,2*i-1:2*i)=K;
end
if(i~=1)
A(2*i-3:2*i-2,2*i-1:2*i)=K;
end
B=[B {A}];
end
A=B


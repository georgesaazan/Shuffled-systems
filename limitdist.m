function a=limitdist();
%Obtain the stationary probability distribution
%vector p of an irreducible, recurrent Markov
%chain by state reduction. P is the transition
%probabilities matrix of a discrete-time Markov
%chain or the generator matrix Q.
% syms a12;syms a13;syms a21;syms a23;syms a31;syms a32;a11=1-a12-a13;a22=1-a21-a23;a33=1-a31-a32;
a12=0.3;a21=0.4;
%P=[1-a a;q 1-q]
%P=[0 ,a21/(a12+a21) ,a12/(a12+a21);a12 1-a12 0;a21 0 1-a21]
P=[0 0 1-a12 a12;0 0 a21 1-a21;0 a12 1-a12 0;a21 0 0 1-a21]
%P=[0.1,0.2,0.3,0.4;0.1,0.2,0.3,0.4;0.1,0.2,0.3,0.4;0.1,0.2,0.3,0.4];

%{
P=[ 0 0 0 a11 a12 a13 0 0 0 0 0 0;
    0 0 0 a21 a22 a23 0 0 0 0 0 0;
    0 0 0 a31 a32 a33 0 0 0 0 0 0;
    0 0 0 a11 0 0 0 a12 0 0 0 a13;
    0 0 0 0 a22 0 a21 0 0 a23 0 0;
    0 0 0 0 0 a33 0 0 a32 0 a31 0;
    0 0 a13 0 0 0 a11 a12 0 0 0 0;
    0 0 a23 0 0 0 a21 a22 0 0 0 0;
    a21 0 0 0 0 0 0 0 a22 a23 0 0;
    a31 0 0 0 0 0 0 0 a32 a33 0 0;
    0 a12 0 0 0 0 0 0 0 0 a11 a13;
    0 a32 0 0 0 0 0 0 0 0 a31 a33]
 %}
%A=[zeros(6,6),[0 0 1;0 1 0;0 0 1;1 0 0;0 1 0;1 0 0],zeros(6,6)];
%B=[zeros(3,6),diag([a11,a22,a33],0),[a12,a13,0,0,0,0;0 0 a21 a23 0 0;0 0 0 0 a31,a32]];
%C=[diag([a23,a32,a13,a31,a12,a21],0),zeros(6,3),diag([1-a23,1-a32,1-a13,1-a31,1-a12,1-a21],0)];
%P=[A;B;C];
[ns ms]=size(P);
n=ns;
a=zeros(n);
%a=sym(zeros(n));
while n>1
n1=n-1;
s=sum(P(n,1:n1));
P(1:n1,n)=P(1:n1,n)/s;
n2=n1;
while n2>0
P(1:n1,n2)=P(1:n1,n2)+P(1:n1,n)*P(n,n2);
n2=n2-1;
end
n=n-1;
end
%backtracking
a(1)=1;
j=2;
while j<=ns
j1=j-1;
a(j)=sum(a(1:j1).*(P(1:j1,j))');
j=j+1;
end
a=a/(sum(a))
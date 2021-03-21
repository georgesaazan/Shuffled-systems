%This outputs \kappa^{\theta}(t) for a shuffled switched signal according
%to Definition 2.1.
function ka =kappa(A,n)
ka=[];
k=0;
a=0;
j=zeros(1,n);
while(isempty(A)==0)
    for i=1:n
     temp=find(A==i,1,'first');
     if(isempty(temp)~=1)
j(i)=temp;
     else 
a=1;
     end
    end
     if(a~=1)
m=max(j);
     else
         m=length(A);
     end
A=A(m+1:end);
ka=[ka k*ones(1,m)];
k=k+1;
end
end
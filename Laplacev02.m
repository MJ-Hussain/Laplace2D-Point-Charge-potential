% Laplace Equation
clear; clc;


%Nubmer of mesh points
N=51;
%Tolerance for the problem
tolerence=1E-03;

%length of one side of region 2m
l=2; %m
%mesh size
dx=l/N;



%Charge 
q=1E-3;

V_o=zeros(N);

%Boundary Condition

n=1; s=0; e=1; w=0; 
%N BC
V_o(1,:)=n;
%S BC
V_o(N,:)=s;
%E BC
V_o(:,1)=e;
%W BC
V_o(:,N)=w;

%Point charge potential
V_o((N-1)/2,(N-1)/2+1)=9*10^9*q/dx;
V_o((N-1)/2+1,(N-1)/2)=9*10^9*q/dx;
V_o((N-1)/2+1,(N-1)/2+2)=9*10^9*q/dx;
V_o((N-1)/2+2,(N-1)/2+1)=9*10^9*q/dx;




% call for solution using Relaxation method--> Jacobi method
[V_j,it_j,t_j]=Jacobi(V_o,tolerence,N);
fprintf("Jacobi method took %i iterations to converge in %6.4f seconds.\n",it_j,t_j)
figure(1)
contourf(V_j,100,'LineStyle','none')
title('Potential using Jacobi method')

% call for solution using Guass-Seidel method
%relaxation parameter=1 for GS method
Wp=1;
[V_R,it_R,t_R]=GS_SOR(V_o,tolerence,Wp,N);
fprintf("Guass-Seidel method took %i iterations to converge in %6.4f seconds.\n",it_R,t_R)
figure(2)
contourf(V_R,100,'LineStyle','none')
title('Potential using Guass-Seidel method')


%Finding optimum Relaxation parameter
Wp=1:0.1:1.9;
t_mat=zeros(1,length(Wp));
it_mat=zeros(1,length(Wp));
for wi=1:length(Wp)
[V_R,it_R,t_R]=GS_SOR(V_o,tolerence,Wp(wi),N);
t_mat(wi)=t_R;
it_mat(wi)=it_R;
end
[min_val,indx]=min(t_mat);
%Optimun value of W is Wp(indx)=> for which time taken is minimun
%call for solution at optimum Wp
[V_R,it_R,t_R]=GS_SOR(V_o,tolerence,Wp(indx),N);
fprintf("Successive Over-Relaxation method took %i iterations to converge in %6.4f seconds.\n",it_R,t_R)
figure(3)
contourf(V_R,100,'LineStyle','none')
title('Potential using Successive Over-Relaxation method')


% Relaxation method--> Jacobi method
function [V,it,time]=Jacobi(V_in,tol,N)
t_start=cputime;


diff=999;
it=0;

V=V_in;
while diff>tol
  it=it+1;
  V_prev=V;
  for i=2:N-1
    for j=2:N-1
      V(i,j)=0.25*(V_prev(i+1,j)+V_prev(i-1,j)+V_prev(i,j+1)+V_prev(i,j-1));
    end
  end
  diff=max(abs(sum(V)-sum(V_prev)));
  
end
t_stop=cputime;
time=t_stop-t_start;

end

% Gauss-seidel method and over relaxation method
function [V,it,time]=GS_SOR(V_in,tol,Wp,N)
t_start=cputime;

diff=999;
it=0;

V=V_in;
while diff>tol
  it=it+1;
  V_prev=V;
  for i=2:N-1
    for j=2:N-1
      V(i,j)=(1-Wp)*V_prev(i,j)+0.25*Wp*(V_prev(i+1,j)+V(i-1,j)+V_prev(i,j+1)+V(i,j-1));
    end
  end
  diff=max(abs(sum(V)-sum(V_prev)));
  
end
t_stop=cputime;
time=t_stop-t_start;

end
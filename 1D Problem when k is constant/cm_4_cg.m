function u = cm_4_cg(u1,uend,f,k,N)

%% This is a function file to solve Helmholtz equation in the form u''+ku = f ; 
%% with boundries at x=0 and x=pi
% u1 = boundary condition at x=0 % example: u1=1; 
% uend = boundary condition at x=pi % example: uend=-1; 
% f = function f on RHS 
        %example for constant wavenumber: f=@(x) 3599*cos(x);
% k = wavenumber 
        %example for constant wavenumber: k= 3600;
% N= the number of grid points (number of intervals+1); so that i=1,2,...N
        %example: N=11; 

h=(pi-0)/(N-1);
xi=0:h:pi;
A=sparse(N,N);
a = 5*k/6 - 2/(h^2);         
b = 1/(h^2) + k/12;            
A = diag(a * ones(N,1)) + diag(b * ones(N-1,1), -1)  + diag(b * ones(N-1,1), 1);
A(1,:) = 0; A(1,1) = 1;
A(end,:) = 0; A(end,end) = 1;
B=sparse(N,1);
for i=2:N-1
B(i)=(1/(12))*f(xi(i+1))+(5/(6))*f(xi(i))+(1/(12))*f(xi(i-1));
end
B(1)=u1; 
B(end)=uend;
u=cgs(A,B,10^(-15),100000);
end
function u = cm_6_dm(u1,uend,f,k,N)

%% This is uses compact dicretization - Sixth order accurate - 
%% This assumes all the required additional boundary conditions are known excalty
%% This is a function file to solve Helmholtz equation in the form u''+ku = f ; 
%% with boundries at x=0 and x=pi
% u1 = boundary condition at x=0 % example: u1=1; 
% uend = boundary condition at x=pi % example: uend=-1; 
% f = function f on RHS 
        %example for constant wavenumber: f=@(x) 3599*cos(x);
        %example for variable wavenumber: f=@(x) -cos(x) + 3600*(1+x^2)*cos(x); 
% k = wavenumber 
        %example for constant wavenumber: k=@(x) 3600
        %example for variable wavenumber: k=@(x) 3600*(1+x^2);
% N= the number of grid points (number of intervals+1); so that i=1,2,...N
        %example: N=11; 

h=(pi-0)/(N-1);
%need extended domain
xi=0-h:h:pi+h;
A=sparse(N+2,N+2);
A(1,1)=1;
A(end,end)=1;
A(2,2)=1;
A(end-1,end-1)=1;
for i=3:N
   A(i,i) =  11*k(xi(i))/15 - 34/(20*(h^2));
   A(i,i-1) = 4/(5*(h^2))+2*k(xi(i-1))/15;
   A(i,i+1) = 4/(5*(h^2))+2*k(xi(i+1))/15;
   A(i,i-2) = 1/(20*(h^2));
   A(i,i+2) = 1/(20*(h^2));
end
B=sparse(N+2,1);
B(2)=u1; 
B(end-1)=uend;
%Additional boundary from exact solution
B(1) = cos(xi(1));
B(end) = cos(xi(end));
for i=3:N
B(i)=(2/(15))*f(xi(i+1))...
    +(11/(15))*f(xi(i))...
    +(2/(15))*f(xi(i-1));
end
u=A\B; u= u(2:end-1);
end

function u = ccm_6_bc4_dm(u1,uend,f,k,N)

%% This is uses compact combined dicretization - Sixth order accurate - 
%% When required additional boundary conditions are unknown, 
   %% this method calculates the additional required boundary using one sided formula
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
P=sparse(N+2,N+2);
    P(1,1)=1; P(N+2,N+2)=1;
    P(2,2)=1; P(N+1,N+1)=1;
    for i=3:N
       P(i,i) =  1;
       P(i,i-1) = 2/11;
       P(i,i+1) = 2/11;
    end
Q=sparse(N+2,N+2);
    Q(1,1)=k(xi(1));
    Q(2,2)=k(xi(2));
    Q(N+1,N+1)=k(xi(N+1));
    Q(N+2,N+2)=k(xi(N+2));
    for i=3:N
       Q(i,i) =  -51/(22*h^2);
       Q(i,i-1) = 12/(11*h^2);
       Q(i,i+1) = 12/(11*h^2);
       Q(i,i-2) = 3/(44*h^2);
       Q(i,i+2) = 3/(44*h^2);
    end
D=sparse(N+2,N+2);
    for i=1:N+2
       D(i,i) = k(xi(i)); 
    end
F=sparse(N+2,1);
    for i=1:N+2
        F(i)=f(xi(i));
    end
R=sparse(N+2,1);
    R(1) = f(xi(1)); R(2) = f(xi(2)); 
    R(N+1) = f(xi(N+1)); R(N+2) = f(xi(N+2));
A = Q+P*D; 
A(1,1) = 1; A(2,2) = 1; 
A(N+1,N+1) = 1; A(N+2,N+2) = 1;
B = P*F-R; 
B(2) = u1; B(N+1) = uend;
%Additional boundary from fourth order solution
u4=cm_4_dm(u1,uend,f,k,N);
B(1) = 7*u4(1) - 21*u4(2) + 35*u4(3) - 35*u4(4) + 21*u4(5) -7*u4(6) +u4(7);
B(end) = 7*u4(end) - 21*u4(end-1) + 35*u4(end-2) - 35*u4(end-3) + 21*u4(end-4) -7*u4(end-5) +u4(end-6);
u=A\B;u= u(2:end-1);
end

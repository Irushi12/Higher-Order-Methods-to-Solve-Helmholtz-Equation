function u = ccm_4_dm(u1,uend,f,k,N)

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
xi=0:h:pi;
P=sparse(N,N);
    P(1,1)=1; P(N,N)=1;
    for i=2:N-1
       P(i,i) =  1;
       P(i,i-1) = 1/10;
       P(i,i+1) = 1/10;
    end
Q=sparse(N,N);
    Q(1,1)=k(xi(1));
    Q(N,N)=k(xi(N));
    for i=2:N-1
       Q(i,i) =  -12/(5*h^2);
       Q(i,i-1) = 6/(5*h^2);
       Q(i,i+1) = 6/(5*h^2);
    end
D=sparse(N,N);
    for i=1:N
       D(i,i) = k(xi(i)); 
    end
F=sparse(N,1);
    for i=1:N
        F(i)=f(xi(i));
    end
R=sparse(N,1);
    R(1) = f(xi(1)); R(N) = f(xi(N));
A = Q+P*D; A(1,1) = 1; A(N,N) = 1;
B = P*F-R; B(1) = u1; B(N) = uend;
u=A\B;
end
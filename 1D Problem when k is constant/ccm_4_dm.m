function u = ccm_4_dm(u1,uend,f,k,N)

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
P=sparse(N,N);
    P = diag(ones(N,1)) + diag((1/10)*ones(N-1,1), -1) + diag((1/10)*ones(N-1,1), 1);
    P(1,:) = 0; P(1,1) = 1;
    P(N,:) = 0; P(N,N) = 1; 
Q=sparse(N,N);
    Q = diag((-12/(5*h^2))*ones(N,1)) + diag((6/(5*h^2))*ones(N-1,1), -1) + diag((6/(5*h^2))*ones(N-1,1), 1);
    Q(1,:) = 0;
    Q(N,:) = 0;
    Q(1,1)=-k;Q(N,N)=-k;
F=sparse(N,1);
    for i=1:N
        F(i)=f(xi(i));
    end
R=sparse(N,1);
    R(1) = f(xi(1)); R(N) = f(xi(N));
A = Q+P*k; A(1,1) = 1; A(N,N) = 1;
B = P*F-R; B(1) = u1; B(N) = uend;
u=A\B;
end
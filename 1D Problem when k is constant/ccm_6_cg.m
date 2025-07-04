﻿function u = ccm_6_cg(u1,uend,f,k,N)

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
%need extended domain
xi=0-h:h:pi+h;
P=sparse(N+2,N+2);
P = diag(ones(N+2,1)) + diag((2/11)*ones(N+1,1), -1) + diag((2/11)*ones(N+1,1), 1);
    P([1 2 N+1 N+2], :) = 0;
    P(1,1)=1; P(N+2,N+2)=1;
    P(2,2)=1; P(N+1,N+1)=1;
Q=sparse(N+2,N+2);
Q = diag((-51/(22*h^2)) * ones(N+2,1)) + diag((12/(11*h^2)) * ones(N+1,1), -1) ...
  + diag((12/(11*h^2)) * ones(N+1,1), 1) + diag((3/(44*h^2))  * ones(N,1), -2) ...
  + diag((3/(44*h^2))  * ones(N,1), 2);
    Q([1 2 N+1 N+2], :) = 0;
    Q(1,1)=-k;
    Q(2,2)=-k;
    Q(N+1,N+1)=-k;
    Q(N+2,N+2)=-k;
F=sparse(N+2,1);
    for i=1:N+2
        F(i)=f(xi(i));
    end
R=sparse(N+2,1);
    R(1) = f(xi(1)); R(2) = f(xi(2)); 
    R(N+1) = f(xi(N+1)); R(N+2) = f(xi(N+2));
A = Q+P*k; 
A(1,1) = 1; A(2,2) = 1; 
A(N+1,N+1) = 1; A(N+2,N+2) = 1;
B = P*F-R; 
B(1) = cos(xi(1)); B(N+2) = cos(xi(N+2));
B(2) = u1; B(N+1) = uend;
u=cgs(A,B,10^(-15),100000); u= u(2:end-1);
end
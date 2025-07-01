function [u] = adi6_bc4_dm_cg(N,ux1,uxend,uy1,uyend,f,k)

%% This is a function file to solve Helmholtz equation 
%% in the form u_xx+u_yy+ku = f ; 
%% with boundries at x=0, x=pi/2, y=0 and y=pi/2
% ux1 = boundary condition at y=0 % example: ux1 = @(x) sin(x); 
% uxend = boundary condition at y=p/2i % example: uxend = @(x) -sin(x); 
% uy1 = boundary condition at x=0 % example: uy1 = @(y) 0; 
% uyend = boundary condition at x=pi/2 % example: uyend = @(y) cos(2*y);; 
% f = function f on RHS 
        %example for constant wavenumber: f=@(x,y) 25000*sin(x)*cos(2*y);
% k = wavenumber 
        %example for constant wavenumber: k= 25000;
% N= the number of grid points in one dorection  (number of intervals+1); 
        %so that i,j=1,2,...N example: N=11; 
%---------------------------------------------------------------------------------
%% Improved CM-6-CG-2D by using an initial guess
%% This uses compact dicretization - Sixth order accurate
%% When required additional boundary conditions are unknown, 
    %% this method calculates the additional required boundary using one sided formula
%% Then use adi6_bc4_dm_extendedsol to calculate an initial guess
%% The initial guess is in the required extended domain
%% Solve Au=b by using CG algorithm

h = pi/(2*(N-1));
% need extended domain
x = (0-h:h:pi/2+h);  
y = (0-h:h:pi/2+h); 
%% additional boundary conditions using fourth order solutions
u4 = cm_4_dm_2d(N,ux1,uxend,uy1,uyend,f,k);
ux0 = 7*u4(:,1) -21*u4(:,2) + 35*u4(:,3) -35*u4(:,4) + 21*u4(:,5) -7*u4(:,6) +u4(:,7);
uxend2 = 7*u4(:,end) -21*u4(:,end-1) + 35*u4(:,end-2) -35*u4(:,end-3) + 21*u4(:,end-4) -7*u4(:,end-5) +u4(:,end-6);
u4extended = [ux0 u4 uxend2];
uy0 = 7*u4extended(1,:) -21*u4extended(2,:) + 35*u4extended(3,:) ...
    -35*u4extended(4,:) + 21*u4extended(5,:) -7*u4extended(6,:) +u4extended(7,:);
uyend2 = 7*u4extended(end,:) -21*u4extended(end-1,:) + 35*u4extended(end-2,:) ...
    -35*u4extended(end-3,:) + 21*u4extended(end-4,:) -7*u4extended(end-5,:) +u4extended(end-6,:);
ux0 = [uy0(1) ux0' uyend2(1)];
uxend2 = [uy0(end) uxend2' uyend2(end)];
%% set up the numbering
num = 1;
for i=1:N+2
    for j=1:N+2   
        Number(i,j) = num;     
        num = num + 1;  
    end
end
%% Matrix A
A = sparse((N+2)*(N+2),(N+2)*(N+2));
% set up the  boundary
for i=1:N+2
    ii = Number(i,1);       A(ii,ii) = 1;
    ii = Number(i,N+2);     A(ii,ii) = 1;
    ii = Number(1,i);       A(ii,ii) = 1;
    ii = Number(N+2,i);     A(ii,ii) = 1;
    ii = Number(i,2);       A(ii,ii) = 1;
    ii = Number(i,N+1);     A(ii,ii) = 1;
    ii = Number(2,i);       A(ii,ii) = 1;
    ii = Number(N+1,i);     A(ii,ii) = 1;
end
for i=3:N
    for j=3:N
        ii = Number(i,j);
        A(ii,ii) = (121*h^2*k/225) - (187/(75)); 
        A(ii,Number(i+1,j)) = (22*h^2*k/225) + (9/(25)); 
        A(ii,Number(i,j+1)) = (22*h^2*k/225) + (9/(25));  
        A(ii,Number(i,j-1)) = (22*h^2*k/225) + (9/(25));        
        A(ii,Number(i-1,j)) = (22*h^2*k/225) + (9/(25));  
        A(ii,Number(i+2,j)) = (11/(300));  
        A(ii,Number(i,j+2)) = (11/(300));  
        A(ii,Number(i,j-2)) = (11/(300));        
        A(ii,Number(i-2,j)) = (11/(300));  
        A(ii,Number(i+1,j+1)) = (4*h^2*k/225) + (16/(75)); 
        A(ii,Number(i+1,j-1)) = (4*h^2*k/225) + (16/(75)); 
        A(ii,Number(i-1,j+1)) = (4*h^2*k/225) + (16/(75)); 
        A(ii,Number(i-1,j-1)) = (4*h^2*k/225) + (16/(75)); 
        A(ii,Number(i+2,j+1)) = 1/(150);
        A(ii,Number(i+2,j-1)) = 1/(150);
        A(ii,Number(i-2,j+1)) = 1/(150);
        A(ii,Number(i-2,j-1)) = 1/(150);
        A(ii,Number(i+1,j+2)) = 1/(150);
        A(ii,Number(i+1,j-2)) = 1/(150);
        A(ii,Number(i-1,j+2)) = 1/(150);
        A(ii,Number(i-1,j-2)) = 1/(150);
    end
end

%% Matrix B
b=sparse(N+2,N+2);
for i=1:N+2
    for j=1:N+2
        b(i,j)= f(x(i),y(j));
    end
end
B = sparse((N+2)*(N+2),1);
% set up the boundary
for i=1:N+2
    ii = Number(i,2);     B(ii,1) = ux1(x(i));
    ii = Number(i,N+1);   B(ii,1) = uxend(x(i));
    jj = Number(2,i);     B(jj,1) = uy1(y(i));
    jj = Number(N+1,i);   B(jj,1) = uyend(y(i));
    ii = Number(i,1);     B(ii,1) = ux0(i);
    ii = Number(i,N+2);   B(ii,1) = uxend2(i);
    jj = Number(1,i);     B(jj,1) = uy0(i);
    jj = Number(N+2,i);   B(jj,1) = uyend2(i);
end
for i=3:N
    for j=3:N
        ii = Number(i,j);
        B(ii,1)= (22*h*h/225)*b(i+1,j) ...
            +(22*h*h/225)*b(i-1,j)  ...
            +(22*h*h/225)*b(i,j+1)  ...
            +(22*h*h/225)*b(i,j-1)  ...
            +(4*h*h/225)*b(i+1,j+1)   ...
            +(4*h*h/225)*b(i-1,j+1)   ...
            +(4*h*h/225)*b(i+1,j-1)   ...
            +(4*h*h/225)*b(i-1,j-1)   ...
            +(121*h*h/225)*b(i,j);
    end
end
%initial guess from ADI
U = adi6_dm_extendedsol(N,ux1,uxend,uy1,uyend,f,k);
for i=1:N+2
    for j=1:N+2
        ii = Number(i,j); 
        Ui(ii) = U(i,j);  
    end
end
U = Ui';
%Ui=cgs(A,B,10^(-15),100000,spdiags([0 1 0],-1:1,(N+3)*(N+3),(N+3)*(N+3)),spdiags([0 1 0],-1:1,(N+3)*(N+3),(N+3)*(N+3)),U);
Ui = cgs(A,B,10^(-15),100000,[],[],U);
%conert back to N by N
for i=1:N+2
    for j=1:N+2
        ii = Number(i,j); 
        u(i,j) = Ui(ii);  
    end
end
u = u(2: end-1 , 2:end-1);
end

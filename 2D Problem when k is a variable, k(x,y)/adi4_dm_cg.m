function [u] = adi4_dm_cg(N,ux1,uxend,uy1,uyend,f,k)

%% This is a function file to solve Helmholtz equation 
%% in the form u_xx+u_yy+ku = f ; 
%% with boundries at x=0, x=pi/2, y=0 and y=pi/2
% ux1 = boundary condition at y=0 % example: ux1 = @(x) sin(x); 
% uxend = boundary condition at y=p/2i % example: uxend = @(x) -sin(x); 
% uy1 = boundary condition at x=0 % example: uy1 = @(y) 0; 
% uyend = boundary condition at x=pi/2 % example: uyend = @(y) cos(2*y);; 
% f = function f on RHS 
        %example for constant wavenumber: f=@(x,y) 25000*sin(x)*cos(2*y);
        %example for variable wavenumber:
            %f=@(x,y) -5*sin(x)*cos(2*y)+(25000)*(1+x^2/2+y^2/2)*sin(x)*cos(2*y); 
% k = wavenumber 
        %example for constant wavenumber: k=@(x,y) 25000
        %example for variable wavenumber: k = @(x,y) 25000*(1+x^2/2+y^2/2);
% N= the number of grid points in one dorection  (number of intervals+1); 
        %so that i,j=1,2,...N example: N=11; 
%---------------------------------------------------------------------------------
%% Improved CM-4-CG-2D by using an initial guess
%% This uses compact dicretization - Fourth order accurate
%% Then use adi4_dm to calculate an initial guess
%% Solve Au=b by using CG algorithm

h = pi/(2*(N-1));
x = (0:h:pi/2);  
y = (0:h:pi/2); 

%% set up the numbering
num = 1;
for i=1:N
    for j=1:N 
        Number(i,j) = num;     
        num = num + 1;  
    end
end
K = sparse(N,N);
for i=1:N
    for j=1:N
        K(i,j)= k(x(i),y(j));
    end
end
A = sparse((N)*(N),(N)*(N));
% set up the  boundary
for i=1:N
    ii = Number(i,1);     A(ii,ii) = 1;
    ii = Number(i,N);     A(ii,ii) = 1;
    ii = Number(1,i);     A(ii,ii) = 1;
    ii = Number(N,i);     A(ii,ii) = 1;
end
for i=2:N-1
    for j=2:N-1
        ii = Number(i,j);
        A(ii,ii) = (25*h*h*K(i,j)/36) -(10/3);
        A(ii,Number(i+1,j)) = (2/3) + (5*h*h*K(i+1,j)/72); 
        A(ii,Number(i,j+1)) = (2/3) + (5*h*h*K(i,j+1)/72); 
        A(ii,Number(i,j-1)) = (2/3) + (5*h*h*K(i,j-1)/72);
        A(ii,Number(i-1,j)) = (2/3) + (5*h*h*K(i-1,j)/72); 
        A(ii,Number(i+1,j+1)) =(1/6) + (h*h*K(i+1,j+1)/144);
        A(ii,Number(i+1,j-1)) =(1/6) + (h*h*K(i+1,j-1)/144); 
        A(ii,Number(i-1,j+1)) =(1/6) + (h*h*K(i-1,j+1)/144);
        A(ii,Number(i-1,j-1)) =(1/6) + (h*h*K(i-1,j-1)/144);
    end
end
B = sparse((N)*(N),1);
% set up the F for boundary
for i=1:N
    ii = Number(i,1);     B(ii,1) = ux1(x(i));
    ii = Number(i,N);     B(ii,1) = uxend(x(i));
    jj = Number(1,i);     B(jj,1) = uy1(y(i));
    jj = Number(N,i);     B(jj,1) = uyend(y(i));
end
% set up inner values
f1=sparse((N),(N));
for i=1:N
    for j=1:N
        f1(i,j)= f(x(i),y(j));
    end
end
for i=2:N-1
    for j=2:N-1
        ii = Number(i,j);
        B(ii,1)= +(5*h*h/72)*f1(i+1,j) ...
            +(5*h*h/72)*f1(i-1,j)  ...
            +(5*h*h/72)*f1(i,j+1)  ...
            +(5*h*h/72)*f1(i,j-1)  ...
            +(h*h/144)*f1(i+1,j+1)   ...
            +(h*h/144)*f1(i-1,j+1)   ...
            +(h*h/144)*f1(i+1,j-1)   ...
            +(h*h/144)*f1(i-1,j-1)   ...
            +(25*h*h/36)*f1(i,j);
    end
end
%initial guess from ADI
U = adi4_dm(N,ux1,uxend,uy1,uyend,f,k);
for i=1:N
    for j=1:N
        ii = Number(i,j); 
        Ui(ii) = U(i,j);  
    end
end
U = Ui';
%Ui=cgs(A,B,10^(-15),100000,spdiags([0 1 0],-1:1,(N)*(N),(N)*(N)),spdiags([0 1 0],-1:1,(N)*(N),(N)*(N)),U); %this works with MATLAB
Ui=cgs(A,B,10^(-15),100000,[],[],U);

%conert back to N by N
for i=1:N
    for j=1:N
        ii = Number(i,j); 
        u(i,j) = Ui(ii);  
    end
end
end

function [u] = adi4_dm(N,ux1,uxend,uy1,uyend,f,k)

%% This is a function file to solve Helmholtz equation 
%% in the form u_xx+u_yy+ku = f ; 
%% with boundries at x=0, x=pi/2, y=0 and y=pi/2
% ux1 = boundary condition at y=0 % example: ux1 = @(x) sin(x); 
% uxend = boundary condition at y=p/2i % example: uxend = @(x) -sin(x); 
% uy1 = boundary condition at x=0 % example: uy1 = @(y) 0; 
% uyend = boundary condition at x=pi/2 % example: uyend = @(y) cos(2*y);; 
% f = function f on RHS 
        %example for variable wavenumber:
            %f=@(x,y) -5*sin(x)*cos(2*y)+(25000)*(1+x^2/2+y^2/2)*sin(x)*cos(2*y); 
% k = wavenumber 
        %example for variable wavenumber: k = @(x,y) 25000*(1+x^2/2+y^2/2);
% N= the number of grid points in one dorection  (number of intervals+1); 
        %so that i,j=1,2,...N example: N=11; 
%---------------------------------------------------------------------------------
%% This is a supplementry function for ADI4-DM-CG
%% This is based on fourth order compact discretization
%% The solution is not accurate since we have ignored terms in the derivation

h = pi/(2*(N-1));
x = (0:h:pi/2);  
y = (0:h:pi/2); 
% set up the numbering
num = 1;
for i=1:N
    for j=1:N  
        Number(i,j) = num;     
        num = num + 1;  
    end
end
% RHS matrix B
B1=sparse(N,N);
for i=1:N
    for j=1:N
        f2(i,j)= (1/12)*f(x(i),y(j)+h)*k(x(i),y(j)+h)...
                +(5/6)*f(x(i),y(j))*k(x(i),y(j))...
                +(1/12)*f(x(i),y(j)-h)*k(x(i),y(j)-h);
    end
end
for i=2:N-1
    for j=1:N
        B1(i,j)= h*h*h*h*((1/12)*f2(i+1,j)+(5/6)*f2(i,j)+(1/12)*f2(i-1,j));
    end
end
% Converting B(N by N) to B (N*N by 1)
for i=1:N
    for j=1:N
        ii = Number(i,j); 
        B(ii) = B1(i,j);  
    end
end

for i=1:N
    ii = Number(1,i);       B(ii) = (5*h*h*k(x(1),y(i))/6-2)*uy1(y(i))...
                                    +(1+ k(x(1),y(i)+h)*h*h/12)*uy1(y(i)+h)...
                                    +(1+ k(x(1),y(i)-h)*h*h/12)*uy1(y(i)-h);
    ii = Number(N,i);       B(ii) = (5*h*h*k(x(N),y(i))/6-2 )*uyend(y(i))...
                                    +(1+ k(x(N),y(i)+h)*h*h/12)*uyend(y(i)+h)...
                                    +(1+ k(x(N),y(i)-h)*h*h/12)*uyend(y(i)-h); %using given BC and next step2
end

% %  solving for U_star
A1 = sparse(N*N,N*N);
% set up the  boundary
for i=1:N
    ii = Number(1,i);       A1(ii,ii) = 1;
    ii = Number(N,i);       A1(ii,ii) = 1;
    
end

for i=2:N-1
    for j=1:N
        ii = Number(i,j);
        A1(ii,ii) = 5*(k(x(i),y(j))*h^2)/6 - 2;
        A1(ii,Number(i+1,j)) = 1+ (k(x(i+1),y(j))*h^2)/12;      
        A1(ii,Number(i-1,j)) = 1+ (k(x(i-1),y(j))*h^2)/12; 
    end
end

U_star=A1\B';

% solving for U
A = sparse(N*N,N*N);
% set up the  boundary
for i=1:N
    ii = Number(i,1);       A(ii,ii) = 1;
    ii = Number(i,N);       A(ii,ii) = 1;
    ii = Number(1,i);       A(ii,ii) = 1;
    ii = Number(N,i);       A(ii,ii) = 1;
end
for i=2:N-1
    for j=2:N-1
        ii = Number(i,j);
        A(ii,ii) = 5*(k(x(i),y(j))*h^2)/6-2 ;
        A(ii,Number(i,j+1)) = 1+ (k(x(i),y(j+1))*h^2)/12 ;      
        A(ii,Number(i,j-1)) = 1+ (k(x(i),y(j-1))*h^2)/12 ; 
    end
end
B_col = U_star;
for i=1:N
    ii = Number(i,1);       B_col(ii) = ux1(x(i));
    ii = Number(i,N);       B_col(ii) = uxend(x(i));
    ii = Number(1,i);       B_col(ii) = uy1(y(i));
    ii = Number(N,i);       B_col(ii) = uyend(y(i));
end

U = A\B_col;
%conert back to N by N
for i=1:N
    for j=1:N
        ii = Number(i,j); 
        u(i,j) = U(ii);  
    end
end
end

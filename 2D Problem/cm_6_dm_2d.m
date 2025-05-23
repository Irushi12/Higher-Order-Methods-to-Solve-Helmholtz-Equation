function [u] = cm_6_dm_2d(N,ux1,uxend,uy1,uyend,f,k)

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

h = pi/(2*(N-1));
x = [0-h:h: (pi/2)+h]; % extended x
y = [0-h:h: (pi/2)+h]; % extended y
%% set up the numbering
num = 1;
for i=1:N+2
    for j=1:N+2   
        Number(i,j) = num;     
        num = num + 1;  
    end
end
K = sparse(N+2,N+2);
for i=1:N+2
    for j=1:N+2
        K(i,j)= k(x(i),y(j));
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
        A(ii,ii) = (121*h^2*K(i,j)/225) - (187/(75)); 
        A(ii,Number(i+1,j)) = (22*h^2*K(i+1,j)/225) + (9/(25)); 
        A(ii,Number(i,j+1)) = (22*h^2*K(i,j+1)/225) + (9/(25));  
        A(ii,Number(i,j-1)) = (22*h^2*K(i,j-1)/225) + (9/(25));        
        A(ii,Number(i-1,j)) = (22*h^2*K(i-1,j)/225) + (9/(25));  
        A(ii,Number(i+2,j)) = (11/(300));  
        A(ii,Number(i,j+2)) = (11/(300));  
        A(ii,Number(i,j-2)) = (11/(300));        
        A(ii,Number(i-2,j)) = (11/(300));  
        A(ii,Number(i+1,j+1)) = (4*h^2*K(i+1,j+1)/225) + (16/(75)); 
        A(ii,Number(i+1,j-1)) = (4*h^2*K(i+1,j-1)/225) + (16/(75)); 
        A(ii,Number(i-1,j+1)) = (4*h^2*K(i-1,j+1)/225) + (16/(75)); 
        A(ii,Number(i-1,j-1)) = (4*h^2*K(i-1,j-1)/225) + (16/(75)); 
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
    % additional boundary condition from exact solution
    ii = Number(i,1);     B(ii,1) = sin(x(i))*cos(2*(-h)); 
    ii = Number(i,N+2);   B(ii,1) = sin(x(i))*cos(2*(pi/2+h));
    jj = Number(1,i);     B(jj,1) = sin(-h)*cos(2*y(i));
    jj = Number(N+2,i);   B(jj,1) = sin(pi/2+h)*cos(2*y(i));
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
%% solving
U=A\B;
% convert u back to 2d vector unew
for i=1:N+2
    for j=1:N+2
        ii = Number(i,j); 
        u(i,j) = U(ii);  
    end
end
u = u(2:end-1,2:end-1);
end
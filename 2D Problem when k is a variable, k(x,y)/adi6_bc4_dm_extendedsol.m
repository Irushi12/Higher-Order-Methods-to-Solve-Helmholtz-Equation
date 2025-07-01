function [u] = adi6_bc4_dm_extendedsol(N,ux1,uxend,uy1,uyend,f,k)

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
%% This is a supplementry function for ADI6-DM-CG
%% This is based on sixth order compact discretization
%% When required additional boundary conditions are unknown, 
    %% this method calculates the additional required boundary using one sided formula
%% The solution is not accurate since we have ignored terms in the derivation
%% The solution is in the extended domain, since it is required in ADI6-DM-CG

h = pi/(2*(N-1));
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
%% RHS matrix B
B1=sparse(N+2,N+2);
for i=1:N+2
    for j=1:N+2
        f2(i,j)= (2/15)*f(x(i),y(j)+h)+(11/15)*f(x(i),y(j))+(2/15)*f(x(i),y(j)-h);
    end
end
for i=2:N+1
    for j=1:N+2
        B1(i,j)= (2/15)*f2(i+1,j)/k(x(i+1),y(j))+(11/15)*f2(i,j)/k(x(i),y(j))+(2/15)*f2(i-1,j)/k(x(i-1),y(j));
    end
end
% Converting B(N by N) to B (N*N by 1)
for i=1:N+2
    for j=1:N+2
        ii = Number(i,j); 
        B(ii) = B1(i,j);  
    end
end
for i=3:N
    ii = Number(1,i);       
    B(ii) = (1-2*((1/(k(x(1),y(i))*h^2)) + (2/15))+6*(1/(20*h^2*k(x(1),y(i)))))*uy0(i)...
        +(((1/(k(x(1),y(i)+h)*h^2)) + (2/15)) ...
        - 4*(1/(20*h^2*k(x(1),y(i)+h))))*uy0(i+1)+(((1/(k(x(1),y(i)-h)*h^2)) + (2/15)) ...
        - 4*(1/(20*h^2*k(x(1),y(i)-h))))*uy0(i-1)...
        +(1/(20*h^2*k(x(1),y(i)+2*h)))*uy0(i+2)...
        +(1/(20*h^2*k(x(1),y(i)+2*h)))*uy0(i-2);
    ii = Number(N+2,i);     
    B(ii) = (1-2*((1/(k(x(i),y(i))*h^2)) + (2/15))+6*(1/(20*h^2*k(x(i),y(i)))))*uyend2(i)...
        +(((1/(k(x(i),y(i)+h)*h^2)) + (2/15)) ...
        - 4*(1/(20*h^2*k(x(i),y(i)+h))))*uyend2(i+1)+(((1/(k(x(i),y(i)-h)*h^2)) + (2/15)) ...
        - 4*(1/(20*h^2*k(x(i),y(i)-h))))*uyend2(i-1) ...
        +(1/(20*h^2*k(x(i),y(i)+2*h)))*uyend2(i+2) ...
        +(1/(20*h^2*k(x(i),y(i)+2*h)))*uyend2(i-2);
end
% solving for U_star
A1 = sparse(N+2*N+2,N+2*N+2);
% set up the  boundary
for i=1:N+2
    ii = Number(1,i);       A1(ii,ii) = 1;
    ii = Number(N+2,i);     A1(ii,ii) = 1;
    ii = Number(i,1);       A1(ii,ii) = 1;
    ii = Number(i,N+2);     A1(ii,ii) = 1;
    ii = Number(2,i);       A1(ii,ii) = 1;
    ii = Number(N+1,i);     A1(ii,ii) = 1;
    ii = Number(i,2);       A1(ii,ii) = 1;
    ii = Number(i,N+1);     A1(ii,ii) = 1;
end

for i=3:N
    for j=3:N
        ii = Number(i,j);
        A1(ii,ii) = 1-2*((1/(k(x(i),y(j))*h^2)) + (2/15))+6*(1/(20*h^2*k(x(i),y(j))));
        A1(ii,Number(i+1,j)) = ((1/(k(x(i+1),y(j))*h^2)) + (2/15)) - 4*(1/(20*h^2*k(x(i+1),y(j))));      
        A1(ii,Number(i-1,j)) = ((1/(k(x(i-1),y(j))*h^2)) + (2/15)) - 4*(1/(20*h^2*k(x(i-1),y(j))));
        A1(ii,Number(i+2,j)) = (1/(20*h^2*k(x(i+2),y(j))));      
        A1(ii,Number(i-2,j)) = (1/(20*h^2*k(x(i-2),y(j))));
    end
end
U_star=A1\B';
% solving for U
A = sparse(N+2*N+2,N+2*N+2);
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
        A(ii,ii) = 1-2*((1/(k(x(i),y(j))*h^2)) + (2/15))+6*(1/(20*h^2*k(x(i),y(j))));
        A(ii,Number(i,j+1)) = ((1/(k(x(i),y(j+1))*h^2)) + (2/15)) -4*(1/(20*h^2*k(x(i),y(j+1))));      
        A(ii,Number(i,j-1)) = ((1/(k(x(i),y(j-1))*h^2)) + (2/15)) -4*(1/(20*h^2*k(x(i),y(j-1))));
        A(ii,Number(i,j+2)) = (1/(20*h^2*k(x(i),y(j+2))));      
        A(ii,Number(i,j-2)) = (1/(20*h^2*k(x(i),y(j+2))));
    end
end
B_col = U_star;
for i=1:N+2
    ii = Number(i,1);       B_col(ii) = ux0(i); 
    ii = Number(i,N+2);     B_col(ii) = uxend2(i); 
    ii = Number(1,i);       B_col(ii) = uy0(i); 
    ii = Number(N+2,i);     B_col(ii) = uyend2(i); 
    ii = Number(i,2);       B_col(ii) = ux1(x(i)); 
    ii = Number(i,N+1);     B_col(ii) = uxend(x(i));
    ii = Number(2,i);       B_col(ii) = uy1(y(i));
    ii = Number(N+1,i);     B_col(ii) = uyend(y(i));
end
U = A\B_col;
%conert back to N by N
for i=1:N+2
    for j=1:N+2
        ii = Number(i,j); 
        u(i,j) = U(ii);  
    end
end
end

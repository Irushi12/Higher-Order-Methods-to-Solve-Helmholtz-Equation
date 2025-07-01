% Exact Solution of example 2D problem

function u = u_exact2d(N)
h=(pi/2-0)/(N-1);
xi=0:h:pi/2;
yj=0:h:pi/2;
%% Exact Solution
for i=1:N
    for j=1:N
    u(i,j)=sin(xi(i))*cos(2*yj(j));
    end
end
end

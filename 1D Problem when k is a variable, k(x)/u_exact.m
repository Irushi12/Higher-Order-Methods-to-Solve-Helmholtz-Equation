% Exact Solution of example 1D problem
function u = u_exact(N)
h=(pi-0)/(N-1);
xi=0:h:pi;
%% Exact Solution
for i=1:N
    u(i)=cos(xi(i));
end
end

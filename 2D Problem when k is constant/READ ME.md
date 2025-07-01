
This folder contains the MATLAB files of function files to solve the 2D Helmholtz problem in the form:
## u_xx+u_yy+k^2u(x,y) = f(x,y); (x,y)\in [(a,b) \times (c,d)]; with boundary conditions u(a,y)=u_a(y); u(b,y)=u_b(y); u(x,c)=u_c(x);  u(x,d)=u_d(x);
### Direct Methods - Solve the Au=b by u=A\b directly
1. CM-4-DM-2D: - This uses compact dicretization - Fourth order accurate
2. CM-6-DM-2D - This uses compact dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
3. CM-6-BC4-DM-2D - This uses compact dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula
### Iterative Methods - Solve Au=b by using CG algorithm
1. CM-4-CG-2D - This uses compact dicretization - Fourth order accurate
2. CM-6-CG-2D - This uses compact dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
3. CM-6-BC4-CG-2D - This uses compact dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula
4. ADI4-DM-CG - Improved CM-4-CG-2D by using an initial guess. Calculates the initial guess using an ADI method derived from fourth order discretization - ADI4-DM is needed
5. ADI6-DM-CG - Improved CM-6-CG-2D by using an initial guess. Calculates the initial guess using an ADI method derived from sixth order discretization - ADI6-DM_extendedsol is needed
6. ADI6-BC4-DM-CG - Improved CM-6-BC4-CG-2D by using an initial guess. Calculates the initial guess using an ADI method derived from sixth order discretization - ADI6-BC4-DM-extendedsol is needed

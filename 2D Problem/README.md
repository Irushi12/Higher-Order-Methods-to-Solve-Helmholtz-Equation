This folder contains the MATLAB files of function files to solve the 2D Helmholtz problem in the form:
## u_xx+u_yy+k^2(x,y)u(x,y) = f(x,y); (x,y)\in [(a,b) \times (c,d)]; with boundary conditions u(a,y)=u_a(y); u(b,y)=u_b(y); u(x,c)=u_c(x);  u(x,d)=u_d(x);
Direct Methods - Solve the Au=b by u=A\b directly
1. CM-4-DM: - This is uses compact dicretization - Fourth order accurate
2. CM-6-DM - This is uses compact dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
3. CM-6-BC4-DM - This is uses compact dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula
4. CCM-4-DM - This is uses compact combined dicretization - Fourth order accurate
5. CCM-6-DM - This is uses compact combined dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
6. CCM-6-BC4-DM - This is uses compact combined dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula

Iterative Methods - Solve Au=b by using CG algorithm
1. CM-4-CG: - This is uses compact dicretization - Fourth order accurate
2. CM-6-CG - This is uses compact dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
3. CM-6-BC4-CG - This is uses compact dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula
4. CCM-4-CG - This is uses compact combined dicretization - Fourth order accurate
5. CCM-6-CG - This is uses compact combined dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
6. CCM-6-BC4-CG - This is uses compact combined dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula

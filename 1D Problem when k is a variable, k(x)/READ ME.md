
This folder contains the MATLAB files of function files to solve the 1D Helmholtz problem in the form:
## u''+k^2(x)u(x) = f(x); x\in [a,b]; with boundary conditions u(a)=u0; u(b)=u_N
Direct Methods - Solve the Au=b by u=A\b directly
1. CM-4-DM: - This uses compact dicretization - Fourth order accurate
2. CM-6-DM - This uses compact dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
3. CM-6-BC4-DM - This uses compact dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula
4. CCM-4-DM - This uses compact combined dicretization - Fourth order accurate
5. CCM-6-DM - This uses compact combined dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
6. CCM-6-BC4-DM - This uses compact combined dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula

Iterative Methods - Solve Au=b by using CG algorithm
1. CM-4-CG: - This uses compact dicretization - Fourth order accurate
2. CM-6-CG - This uses compact dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
3. CM-6-BC4-CG - This uses compact dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula
4. CCM-4-CG - This uses compact combined dicretization - Fourth order accurate
5. CCM-6-CG - This uses compact combined dicretization - Sixth order accurate - This assumes all the required additional boundary conditions are known excalty
6. CCM-6-BC4-CG - This uses compact combined dicretization - Sixth order accurate - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula

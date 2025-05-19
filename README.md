# Higher-Order-Methods-to-Solve-Helmholtz-Equation
This Repository contains my master's thesis work on Fourth Order and Sixth Order Methods to Solve 1D and 2D Helmholtz Equations with large wavenumbers.
## 1D Helmholtz Equation
To solve 1D Helmholtz equation, six direct methods and six iterative methods are developed.
### Direct Methods
  1) CM-4-DM:
    - This is uses compact dicretization
    - Fourth order accurate
  2) CM-6-DM
    - This is uses compact dicretization
    - Sixth order accurate
    - This assumes all the required additional boundary conditions are known excalty 
  3) CM-6-BC4-DM
    - This is uses compact dicretization
    - Sixth order accurate
    - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula
  4) CCM-4-DM
    - This is uses compact combined dicretization
    - Fourth order accurate
  5) CCM-6-DM
    - This is uses compact combined dicretization
    - Sixth order accurate
    - This assumes all the required additional boundary conditions are known excalty 
  6) CCM-6-BC4-DM
    - This is uses compact combined dicretization
    - Sixth order accurate
    - When required additional boundary conditions are unknown, this method calculates the additional required boundary using one sided formula


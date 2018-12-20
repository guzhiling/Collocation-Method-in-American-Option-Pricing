# Collocation-Method-in-American-Option-Pricing

Welcome to my project: Radial Basis Function (RBF) Method in Solving Fee-boundary Partial Differential Equation. The aim of this project is to prove the convergence behaviour of RBF method in solving such problems.
RBF method has been very popular in the empirical world, yet the proof of the convergence still remain unsolved. 

## Previous works include:
- Binomial method: Cox, Ross, and Rubinstein (1979)
- Finite difference method: Brennan and Schwartz (1977)
  - Artificial boundary conditions: Han and Wu (1985)
  - Projected SOR method: Wilmott (1993)
  - Front-fixing technique: Wu and Kwok (1997)
  - Penalty method: Zvan, Forsyth, and Vetzal(1998)
  - Far field boundary conditions: Kangro and Nicolaides (2000)
- Radial basis function method: Kansa (1990)
  - Error bounds: Franke and Schaback (1998)
  - Numerical Analysis: Hon and Mao (1999), Hon and Schaback (2001), Larsson and Fornberg (2003)
  - H2 Convergence of Lease-Square Kernel Collocation Method: Cheung, Ling, and Schaback (2017)

## Preliminary Numerical Study: American Call Under Black-Scholes Equation
- Rewrite the original BS equation into finite difference for $\frac{\partial C}{\partial \tau}$ and apply Radial Basis Function Method on BS differential operator. 
- See the presentation slides for details.
- See '1d-BS-hybrid' for Numerical test.

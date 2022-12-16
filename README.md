# Description
In this repository, I present a simple coupled steady-state, ocean-atmosphere 4-box model. It is of a Sarmiento and Toggweiler-type (1985) that can be used to understand the qualitative feedbacks between the ocean and atmosphere systems. 

# Model Equations
I present the simple model equations in this section, which are modified from Emerson and Hamme (2022). We start with the phosphate equations:
```math
\begin{align}
  V_h \frac{\partial P_h}{\partial t} &= T(P_l - P_h) + f_{hd}(P_d - P_h) - J_hS_h, \\
  V_l \frac{\partial P_l}{\partial t} &= T(P_d - P_l) - J_lS_l.
\end{align}
```
Now I present the alkalinity equations ::
```math
\begin{align}
  V_h \frac{\partial A_h}{\partial t} &= T(A_l - A_h) + f_{hd}(P_d - P_h) - J_hS_hr_{a:p}, \\
  V_l \frac{\partial A_l}{\partial t} &= T(P_d - P_l) - J_lS_lr_{a:p}, \\
  \Sigma A &= V_lA_l + V_dA_d + V_hA_h.
\end{align}
```
Next I cover the equations for DIC ::
```math
\begin{align}
  V_h \frac{\partial C_h}{\partial t} &= T(C_l - C_h) + f_{hd}(C_d - C_h) - J_hS_hr_{c:p} + kK_{H,h}S_h(f\text{CO}_2^a - f\text{CO}_2^h), \\
  V_l \frac{\partial C_l}{\partial t} &= T(C_d - C_l) - J_lS_lr_{c:p} + kK_{H,l}S_l(f\text{CO}_2^a - f\text{CO}_2^l), \\
  V_a \frac{\partial C_a}{\partial t} &= kK_{H,h}S_h(f\text{CO}_2^a - f\text{CO}_2^h) + kK_{H,l}S_l(f\text{CO}_2^a - f\text{CO}_2^l), \\
  \Sigma C &= V_lC_l + V_dC_d + V_hC_h + M_a\cdot f\text{CO}_2^a.
\end{align}
```
We continue to follow Emerson and Hamme (2022) and linear the relationship between $A, C$, and $f\text{CO}_2^{sw}$:
```math
\begin{align}
  (A_h - C_h) &= \beta_hf\text{CO}_2^h + \gamma_h, \\
  (A_l - C_l) &= \beta_lf\text{CO}_2^l + \gamma_l. \\
\end{align}
```
# Model Formulation
We assume steady-state, instead of using an interative solver as is done by Sarmiento and Toggweiler (1985). This allows for a simulataneously solution of matrices of the form:
```math
Ax=b
```
which has the solution, assuming that $A$ is invertible, of the form:
```math
x^\ast = A^{-1}b.
```
Therefore, we have the following set of matrix equations, starting with those for phosphate, which are modified from those presented in Emerson and Hamme (2022) to let the concentrations in the high and low latitude boxes be the unknowns that are determined through inversion:
$$ 
M =
\begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}
$$

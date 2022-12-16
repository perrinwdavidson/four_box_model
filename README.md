# Description
In this repository, I present a simple coupled steady-state, ocean-atmosphere 4-box model. It is of a Sarmiento and Toggweiler-type (1985) that can be used to understand the qualitative feedbacks between the ocean and atmosphere systems. 

# Model Equations
I present the simple model equations in this section. We start with the phosphate equations:
\begin{align}
  V_h \frac{\partial P_h}{\partial t} &= T(P_l - P_h) + f_{hd}(P_d - P_h) - J_hS_h \\
  V_l \frac{\partial P_l}{\partial t} &= T(P_d - P_l) - J_lS_l
\end{align}

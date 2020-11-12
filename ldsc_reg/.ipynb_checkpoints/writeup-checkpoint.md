---
title: SIB LDSC
---

For SNP $i$, let $\hat{\theta_i} \sim \mathcal{N}(0, S_i + V)$

Where

$$
S_i = \frac{1}{M}\begin{bmatrix}
\sigma_{i1}^2 & r_{is}\sigma_{i1}\sigma_{i2} \\
r_{is}\sigma_{i1}\sigma_{i2} & \sigma_{i2}^2
\end{bmatrix}
$$
$$
V = \frac{1}{M} \begin{bmatrix}
v_1 & r\sqrt{v_1 v_2}\\
r\sqrt{v_1 v_2} & v_2
\end{bmatrix} 
$$

We drop the $i$ subscript from now on.

Define $D := \sqrt{M}\begin{bmatrix}
\frac{1}{\sigma_1} & 0\\
0 & \frac{1}{\sigma_2}
\end{bmatrix} $

Let $z := D \hat{\theta}$

Then 

$z \sim \mathcal{N} \left (0, \begin{bmatrix}
1 & r_s\\
r_s & 1
\end{bmatrix} + 
 \begin{bmatrix}
\frac{v_1}{\sigma_1^2} & r \frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2}\\
r \frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2} & \frac{v_2}{\sigma_2^2}
\end{bmatrix} \right )$

Let $\Sigma := \begin{bmatrix}
1 & r_s\\
r_s & 1
\end{bmatrix} + 
 \begin{bmatrix}
\frac{v_1}{\sigma_1^2} & r \frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2}\\
r \frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2} & \frac{v_2}{\sigma_2^2}
\end{bmatrix} $

Then we can define the log likelihood $l$ as:

$$
l = - log(2 \pi) - \frac{1}{2} |\Sigma| - \frac{1}{2} z' \Sigma^{-1} z
$$


Let $T := z_1^2 \left( 1 + \frac{v_2}{\sigma_2^2}\right) - 2 \left( r_s + r\frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2}\right) z_1 z_2 + z_2^2 \left( 1 + \frac{v_1}{\sigma_1^2}\right) $

The gradient with respect to the parameters of interest are:

$$
\frac{\partial l}{\partial v_1} = \frac{1}{|\Sigma|} \frac{\partial |\Sigma|}{\partial v_1} + \frac{1}{|\Sigma|} \left(\frac{\partial T}{\partial v_1}  - \frac{T}{|\Sigma|} \frac{\partial |\Sigma|}{\partial v_1} \right)
$$

Where $\frac{\partial |\Sigma|}{\partial v_1} = \frac{1}{\sigma_1^2} \left(1 + \frac{v_2}{\sigma_2^2} \right) - \frac{r}{\sigma_1 \sigma_2} \left(r_s \sqrt{\frac{v_2}{v_1}} + \frac{r v_2}{\sigma_1 \sigma_2} \right)$

and $\frac{\partial T}{\partial v_1} = \frac{z_2^2}{\sigma_1^2} - \frac{r z_1 z_2}{\sigma_1 \sigma_2} \sqrt{\frac{v_2}{v_1}}$


$$
\frac{\partial l}{\partial v_2} = \frac{1}{|\Sigma|} \frac{\partial |\Sigma|}{\partial v_2} + \frac{1}{|\Sigma|} \left(\frac{\partial T}{\partial v_2}  - \frac{T}{|\Sigma|} \frac{\partial |\Sigma|}{\partial v_2} \right)
$$

Where $\frac{\partial |\Sigma|}{\partial v_2} = \frac{1}{\sigma_2^2} \left(1 + \frac{v_1}{\sigma_1^2} \right) - \frac{r}{\sigma_1 \sigma_2} \left(r_s \sqrt{\frac{v_1}{v_2}} + \frac{r v_1}{\sigma_1 \sigma_2} \right)$

and $\frac{\partial T}{\partial v_2} = \frac{z_1^2}{\sigma_2^2} - \frac{r z_1 z_2}{\sigma_1 \sigma_2} \sqrt{\frac{v_1}{v_2}}$

$$
\frac{\partial l}{\partial r} = \frac{1}{|\Sigma|} \frac{\partial |\Sigma|}{\partial r} + \frac{1}{|\Sigma|} \left(\frac{\partial T}{\partial r} - \frac{T}{|\Sigma|}  \frac{\partial |\Sigma|}{\partial r} \right)
$$

Where $\frac{\partial |\Sigma|}{\partial r} = -2 \frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2} \left(r_s + r\frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2} \right)$

and $\frac{\partial T}{\partial r} = -2 \frac{\sqrt{v_1 v_2}}{\sigma_1 \sigma_2} z_1 z_2$

# Conventions comparison

See for Tichy: https://arxiv.org/abs/1410.7687, and for Shchesnovich: https://arxiv.org/abs/1410.1506

## Operators change

* Tichy: ```\hat{A}_{j,|\Phi\rangle}^{\dagger} \rightarrow \hat{U} \hat{A}_{j,|\Phi\rangle}^{\dagger} \hat{U}^{-1}=\sum_{k=1}^{m} U_{j, k} \hat{B}_{k,|\Phi\rangle}^{\dagger}```

* Shchesnovich: ```A spatial unitary network can be defined by an unitary transformation between input $a_{k, s}^{\dagger}(\omega)$ and output $b_{k, s}^{\dagger}(\omega)$ photon creation operators, we set $a_{k, s}^{\dagger}(\omega)=\sum_{l=1}^{M} U_{k l} b_{l, s}^{\dagger}(\omega)$, where $U_{k l}$ is the unitary matrix describing such an optical network.```

This means that if Tichy uses U, then Shchesnovich has in place U^†.

## Scattering matrix

* Tichy: ```Using the mode assignment list $\vec{d}(\vec{s})=\left(d_{1}, \ldots d_{n}\right)$ [49], which indicates the mode in which the $j$ th particle resides, the effective scattering matrix becomes
$$
M=U_{\vec{d}(\vec{r}), \vec{d}(\vec{s})}
$$
where our convention identifies the $j$ th row (column) with the $j$ th input (output) mode, as illustrated in Fig. 1)(a).```

* Shchesnovich: identical

# Global Argument-based Multidimensional Equation Solver

An updated version of a global algorithm is available: [Global Equation Solver](https://github.com/Drakicy/GES)

## Definition

Global Argument-based Multidimensional Equation Solver (GAMES) is a global multidimensional complex zeros and poles finding algorithm based on the argument analysis. The proposed algorithm can be used to solve a parametric equation
$$f(z,{\bf p}) = 0,$$
where $z = x + iy \in \mathbb{C}$ and ${\bf p} \in \mathbb{R}^d$ with $d \in \mathbb{N} \cup \\{ 0 \\}$. The function $f$ is assumed to be meromorphic for any ${\bf p}$ in the inspected region, however, in special cases, the applicability conditions can be broadened further.

The solver is represented as a hybrid class in MATLAB for $d \in \\{ 0, 1\\}$. The selection of the appropriate methods is performed automatically based on the input function $f$.

## Detailed description

For details, see **GAMES.m** and the following article:
1. Viktor A. Frantsuzov, Anton V. Artemyev, "A global argument-based algorithm for finding complex zeros and poles to investigate plasma kinetic instabilities", Journal of Computational and Applied Mathematics, [link](http://dx.doi.org/10.1016/j.cam.2024.116217)

Refer to this publication if the algorithm is used in a scientific work.

# Global Argument-based Multidimensional Equation Solver

## Definition

Global Argument-based Multidimensional Equation Solver (GAMES) is a derivative-free multidimensional complex zeros and poles finding algorithm based on the argument analysis. The proposed algorithm can be used to solve a parametric equation
$$f(z,{\bf p}) = 0,$$
where $z = x + iy \in \mathbb{C}$ and ${\bf p} \in \mathbb{R}^d$ with $d \in \mathbb{N} \cup \\{ 0 \\}$. The function $f$ is assumed to be meromorphic for any ${\bf p}$ in the inspected region, however, in special cases the applicability conditions can be broadened further.

The solver is represented as a hybrid class in MATLAB for $d \in \\{ 0, 1\\}$. The selection of the appropriate methods is performed automatically based on the arguments list of the input function $f$.

## Detailed description

A detailed description of the solver is presented in an article:

**Insert DOI Here**

## Input parameters

### FuncHandle

A handle of the function $f$

*Possible values: a function handle with 1 or 2 arguments*

### RegBound

A boundary of the inspected region for $x$, $y$ and $p$ (optionally)

*Possible values: a real array with 2 or 3 rows of 2 ascending elements each*

### Tol

A global relative tolerance of the approximation

*Possible values: a real number less than 1*

### FuncEvalLim

A limitation for the number of the function evaluation

*Possible values: an integer greater than 8*

## Output parameters

### GAMES handle

A handle of a GAMES class

## Class properties

### Dim

Dimensions of the problem

*Possible values: 2 or 3*

### FuncHandle

A handle of the function $f$

*Possible values: a function handle with 1 or 2 arguments*

### GridNorm

Grid linear normalization coefficients (slope and rise)

*Possible values: a real array with 2 or 3 rows of 2 elements each*

### DT

Delaunay triangulation class

*Possible values: MATLAB Delaunay triangulation class*

### NodesNum

Number of triangulation nodes

*Possible values: an integer greater than 8*

### ApproxNodes

An approximation of the solutions and its classification: 1 is a zero, -1 is a pole and 0 is a noise point.

*Possible values: a complex array with 2 or 3 coordinates columns and 1 classification column*

### Tol

A global relative tolerance of the approximation

*Possible values: a real number less than 1*

## Class methods

### fitTriang

Fits the triangulation with the input parameters

*Method parameters: none*

### visSol

Visualizes the solution

*Method parameters: visParam (a subarray of* $[-1,0,1]$*, optional) selects the type of approximation points to visualize*

## Examples

Every example from the article is linked above is present in the directory. See the article for more information.

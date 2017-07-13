# fenics-rigid-motions
----------------------

This repository accompanies the [paper](https://arxiv.org/abs/1609.09425) **On the Singular Neumann Problem in Linear Elasticity** by _Miroslav Kuchta_, _Kent-Andre Mardal_ and _Mikael Mortensen_ (currently under review for being published in Numerical Linear Algebra with Applications). It contains implementation of three different formulations of the problem in [FEniCS](https://fenicsproject.org/) and [cbc.block](https://bitbucket.org/fenics-apps/cbc.block), namely

  1. Formulation where the unknowns are displacement and Lagrange multiplier
  2. Formulation where the unknowns are displacement, solid pressure and Lagrange multiplier
  3. Formulation based on the complemented energy functional (see, Section 5. of the paper) with the only unknown being the displacement
  
Linear systems arising from the formulations are to be solved iteratively and we provide implementations of mesh independent (and for 2. also parameter independent) preconditioners analyzed in the paper. The code is tested with FEniCS stack version 2017.1.0 and cbc.block commit [8447b45](https://bitbucket.org/fenics-apps/cbc.block/commits/8447b459156459b5a6f29ce8129fc07c2dc644cd?at=master).

----------------------------

Cloning the repo and executing `./test.sh NPROCS`, where `NPROCS` is the number of processes you wish to engage in the computation (4 and higher is recommended for reasonable simulation time), you should get plots similar to those below. These graphs show iteration counts of the solvers for different problem sizes. You can observe that the iteration count is bounded in the problem size. For 2. also the robustness in $\lambda$ parameter is demonstrated. We note that the setup, i.e. computational domain, right hand sides, tolerances, is different from the examples used in the paper.

-----------------------------

<p align="center">
  <img src="https://github.com/MiroK/fenics-rigid-motions/blob/master/.img/primal.png" width="600" height="400" alt="Formulation 2."/>
</p>

<p align="center">
  <img src="https://github.com/MiroK/fenics-rigid-motions/blob/master/.img/mixed.png" width="600" height="400" alt="Formulation 2."/>
</p>

<p align="center">
  <img src="https://github.com/MiroK/fenics-rigid-motions/blob/master/.img/energy.png" width="600" height="400" alt="Formulation 3."/>
</p>


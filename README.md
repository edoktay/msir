# tsir
MATLAB codes for performing three-stage mixed precision iterative refinement


## Included MATLAB files
* **_chop.m, float_params.m, lutx_chop.m,  trisol.m, and roundit.m_** are functions in chop library that simulate half precision. The library and associated functions are available at https://github.com/higham/chop and https://github.com/SrikaraPranesh/LowPrecision\_Simulation.

* **_scale_diag_2side.m_** is a function that performs two-sided diagonal scaling of matrix. It is available at https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels.

* **_lp_matvec.m_** is a function that performs matvec in low precision. It is available at https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels.

* **_sir.m_** is a function that performs LU-based iterative refinement with three precisions.

* **_sgmresir.m_** is a function that performs GMRES-based iterative refinement in three precisions.

* **_gmresir.m_** is a function that performs GMRES-based iterative refinement in three precisions with an extra higher precision in factorization and preconditioner steps of GMRES.

* **_tsir.m_** is a function that performs three-stage iterative refinement in three precisions, switching from SIR to SGMRES-IR to GMRES-IR based on stopping criteria.

* **_tsir_sir.m, tsir_sgmresir.m, tsir_gmresir.m_** are the modified versions of SIR, SGMRES-IR, and GMRES-IR, respectively, called from within the tsir.m function. These functions are not designed to be called alone; if you want to just run SIR, SGMRES-IR, or GMRES-IR, use the functions without the "tsir_" prefix. 

* **_gmres_hs.m, gmres_sd.m, and gmres_dq.m_** are functions that run left-preconditioned GMRES using precisions half/single, single/double, and double/quad, resp. Application of the preconditioned coefficient matrix to a vector and the preconditioner to the right-hand-side vector are performed in the higher precision; other computations performed all use the lower precision.  

* **_gmres_hh.m, gmres_ss.m, and gmres_dd.m_** are functions that run left-preconditioned GMRES using precisions half/half, single/single, and double/double, resp. Application of the preconditioned coefficient matrix to a vector and the preconditioner to the right-hand-side vector are performed in the same precision with other computations.

* **_tsir_test.m_** is an example script for comparing SIR, SGMRES-IR, GMRES-IR, and TSIR (with 3 precisions) on random dense matrices having various condition numbers and various number of nonzero singular values. 

* **_tsir_test_ss.m_** is an example script for comparing SIR, SGMRES-IR, GMRES-IR, and TSIR (with 3 precisions) on matrices in SuiteSparse collection.


## Requirements
* The codes have been developed and tested with MATLAB 2020a.
* tsir_test_ss.m requires ssget MATLAB interface for testing the algorithm on matrices in SuiteSparse collection.
* The codes require the Advanpix Multiprecision Computing Toolbox for extended precision computations. 
A free trial of Advanpix is available for download from https://www.advanpix.com/.



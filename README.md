# RandomizedGMRES

We present the randomized flexible GMRES method with/without deflated restarting. In the construction of Krylov basis vectors, the randomized Gram-Schmidt (RGS) algorithm is employed based on *random sketching*. The RGS process leads the reduction of dimension to yeild decent stability as well as cheap computational costs. Consequently, we can develop the randomized FGMRES with some advantages in stability and computational complexity. Also, we develop the method of deflation in two ways; (i) computing harmonic Ritz eigenpairs and (ii) computing singular vectors. Both deflation in GMRES could improve the convergence of GMRES.

Here, we provide some numerical experiments (written in 'matlab')to show the comparison between our proposed method and others, and numerical performances to solve linear systems arising in CFD simulations. 
1. Comparison of CGS, MGS and RGS
2. Solving Euler equations with the randomized GRMES
3. Solving RANS equations with the randomized FGMRES

For more details of algorithms and the results, please see our paper (TBA). If you have any enquiry, please contact **Yongseok Jang (yongseok.jang@onera.fr)**.

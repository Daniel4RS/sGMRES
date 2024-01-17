# sGMRES
Implementation of the sketched Generalized minimal residual method (sGMRES) for solving nonsingular linear algebraic systems. This program supports Krylov subspace basis generation with (incomplete/ truncated) modified Gram-Schmidt orthogonalization and Chebyshev polynomials. For Chebyshev polynomials, a rectangle containing the spectrum of the matrix must be provided. The overdetermined least squares problem in the GMRES routine can be solved by a thin or full QR factorization. Numerical experiments suggest that the thin QR factorization can lead to increasing residual norms and that the full QR factorization is more stable. The dimension reduction of the QR factorization with sketching can be done with the subsampled random cosine/ Fourier transform as a fast Johnson-Lindenstrauss transformation or, for testing purposes, with a Gaussian essemble.

# Reference:
Yuji Nakatsukasa and Joel A. Tropp. “Fast & Accurate Randomized Algorithms for Linear
Systems and Eigenvalue Problems”. In: arXiv e-prints (Oct. 2021). doi: 10.48550/arXiv.
2111.00113.

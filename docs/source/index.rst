.. role:: underline
    :class: underline

Welcome to DAFMM2Dlib's Documentation!
**************************************

About :math:`\texttt{DAFMM2Dlib}`:
==================================

DAFMM (Directional Algebraic FMM) is an FMM implementation for oscillatory kernels. It differs from the usual FMM, in the sense that the far-field region of boxes is divided into conical regions. The sub-matrices corresponding to the interactions in these conical regions are of low rank.

:math:`\texttt{DAFMM2Dlib}` is a library consisting of DAFMM for kernels in 2D. It can also be used when the underlying kernel is not known, but the matrix entries are known explicitly.

Low-rank approximation of the appropriate blocks is obtained using ACA. The algorithm has been parallelized using OpenMP.

The code is written in C++ and features an easy-to-use interface, where the user provides the following inputs:

- The matrix that needs to be worked that is defined through the method ``getMatrixEntry`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix

- the vector ``b`` to be multiplied to the matrix


Obtains :math:`A x` at a cost of :math:`\mathcal{O}\left(N\log(N)\right)`

Doc Contents
============
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial

Other Links
===========

Learn more about :math:`\texttt{DAFMM2Dlib}` by visiting the

* Code Repository: http://github.com/sivaramambikasaran/DAFMM2D
* Documentation: http://afmm2d.rtfd.io

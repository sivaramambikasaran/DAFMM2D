********
Tutorial
********

For the sake of this tutorial, we are going to be using the ``tutorial.cpp`` file that is listed under ``examples/`` since it demonstrates the features of this library. For the most part, comments in the file demonstrate intended functionality. However, we go over the main functions that may be of interest to a user on this page.

**NOTE**: It is assumed that you have already completed the installation process of getting the dependencies.

Setting Parameters in Makefile2D.mk
-----------------------------------

There are some variables that need to be set by the user at the top of the ``CMakeLists.txt`` file:

- ``SOURCES``: This is the input ``.cpp`` file that needs to be compiled. For this tutorial, it's going to be set to ``examples/tutorial.cpp``.
- ``EXECUTABLE``: This is the name that the final build executable is given. Here we are just going to set is as ``tutorial``.

Running the program::
---------------------

For this particular tutorial, the problem parameters are passed to the executable during runtime. We have the lines::

  nCones_LFR	=	atoi(argv[1]); // Number of cones into which the far-field of each box in the lowest level of High Frequency Regime is divided
  nChebNodes	=	atoi(argv[2]); // nChebNodes*nChebNodes is the size of leaf box
  L						=	atof(argv[3]); // L is the half size of the square domain centered at origin
  kappa				= atof(argv[4]); // wave number of the problem
  TOL_POW			= atoi(argv[5]); // Tolerance of the problem

This means that the first argument would be the Number of cones into which the far-field of each box in the lowest level of High Frequency Regime is divided, the second one would be the square root of the leaf size, the third one would be the  half size of the square domain centered at origin, the fourth one would be the wave number of the problem, and the final argument is approximately the number of digits of accuracy we want. For instance, running ``./tutorial 10000 32 12`` would correspond to solving the problem with parameters :math:`N=1000, MinParticlesInLeaf=32, \epsilon=10^{-12}`.


Creating a Derived Class of ``kernel``:
---------------------------------------

The matrix that needs to be worked with is defined through the method ``getMatrixEntry`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix. For instance, for the :math:`exp(I*kappa*R)/R` kernel, this would be set as::

  kernel_dtype FMM_Matrix::getMatrixEntry(const unsigned i, const unsigned j) {
  	pts2D ri = particles_X[i];
  	pts2D rj = particles_X[j];
  	double R2 = (ri.x-rj.x)*(ri.x-rj.x) + (ri.y-rj.y)*(ri.y-rj.y);
  	double R = sqrt(R2);
  	kernel_dtype out = exp(I*kappa*R)/R;
  	if (R < machinePrecision) {
  		return R;
  	}
  	else {
  		return out;
  	}
  }

The location of nodes is set to be the chebyshev nodes of each leaf box.

Defining vector ``b``:
----------------------

In this tutorial, we have defined ``b``, the vector that is to be multiplied to the matrix, as a random set of points::

  Eigen::VectorXd b = Eigen::VectorXd::Random(N);

Creating the Instance of ``DAFMM2D``:
-------------------------------------

The main operations of this library are carried out through the ``DAFMM2D`` class. The parameters that are taken for the constructor are N, MinParticlesInLeaf, TOL_POW, loc::

  DAFMM2D *dafmm2d = new DAFMM2D(inputs);

where inputs is an instance of struct inputsToDFMM::

  struct inputsToDFMM {
    int nCones_LFR;
    int nChebNodes;
    double L;
    int yes2DFMM;
    int TOL_POW;
  };

where yes2DFMM is a flag which can be set to 0 if one wishes to perform AFMM2D; set it to 1 for DAFMM2D. By default it is 1.

We will now proceed to demonstrate the individual methods available under this class.

``assemble``
^^^^^^^^^^^^

Assemble the matrix in DAFMM structure; i.e. it finds the low rank representation of the appropriate matrix sub-blocks::

  dafmm2d->assemble();

``Mat-Vec product``
^^^^^^^^^^^^^^^^^^^

Multiplies the matrix that is defined through object ``FMM_Matrix`` with the vector ``b``::

  dafmm2d->MatVecProduct(b, DAFMM_Ab);

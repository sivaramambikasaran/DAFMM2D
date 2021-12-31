*************************
Installation and Building
*************************

Downloading the Source
-----------------------

:math:`\texttt{DAFMM2Dlib}` is distributed using the git version control system, and is hosted on Github. The repository can be cloned using::

    git clone https://github.com/sivaramambikasaran/DAFMM2D.git --recursive

The ``--recursive`` flag is argument to ensure that all submodules are also checked out.

Dependencies
-------------

- Eigen Linear Algebra Library (get it `here <https://bitbucket.org/eigen/eigen/>`_)
- (optional) An OpenMP enabled compiler (e.g. gcc4.2 or above) is required to use shared-memory parallelism.

**NOTE**: On MacOS, the default compiler is `clang` which doesn't have OpenMP support. You will have to use g++ to make use of the speedups from OpenMP::

    user@computer DAFMM2D$ brew install g++-8
    user@computer DAFMM2D$ export CXX=g++

Installation
-------------

Manually Installing
^^^^^^^^^^^^^^^^^^^

Move the downloaded ``eigen`` library to the location usr/local/include

Testing
-------

Now, we need to ensure that all the functions of the libraries function as intended. For this purpose, we will be running the script ``examples/testDAFMM2D.cpp``.
Key in the file ```examples/testDAFMM2D.cpp`` as input under ``SOURCES`` in ``DAFMM2Dlib/Makefile2D.mk``. Here you also set the name of the output executable, say ``test_DAFMM2D``, under ``EXECUTABLE``.
Run the ``Makefile`` to get your executable.
To check this on your computer, run the following lines::

    user@computer DAFMM2D$ make -f Makefile2D.mk
    user@computer DAFMM2D$ ./test_DAFMM2D

Building and Executing
----------------------

Key in the required ``.cpp`` to be used as input under ``SOURCES`` in ``DAFMM2Dlib/Makefile2D.mk``. Here you also set the name of the output executable under ``EXECUTABLE``. Then run the ``Makefile`` to get your executable::

  user@computer DAFMM2D$ make -f Makefile2D.mk
  user@computer DAFMM2D$ ./executable

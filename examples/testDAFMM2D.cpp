#include "../DAFMM2D.cpp"

kernel_dtype FMM_Matrix::getMatrixEntry(const unsigned i, const unsigned j) {
	pts2D ri = particles_X[i];//particles_X is a member of base class FMM_Matrix
	pts2D rj = particles_X[j];//particles_X is a member of base class FMM_Matrix
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

int main(int argc, char* argv[]) {
	int nCones_LFR;
	int nChebNodes;
	double L;
  int TOL_POW;
	if(argc < 4)
	{
			std::cout << "All arguments weren't passed to executable!" << std::endl;
			std::cout << "Using Default Arguments:" << std::endl;
			nCones_LFR	=	16;
			nChebNodes	=	6;
			L						=	1.0;
			kappa				= 100.0;
		  TOL_POW			= 8;
	}

	else
	{
		nCones_LFR	=	atoi(argv[1]);
		nChebNodes	=	atoi(argv[2]);
		L						=	atof(argv[3]);
		kappa				= atof(argv[4]);
		TOL_POW			= atoi(argv[5]);
	}

	int yes2DFMM				=	1;
  inputsToDFMM inputs;
  inputs.nCones_LFR = nCones_LFR;
  inputs.nChebNodes = nChebNodes;
  inputs.L = L;
  inputs.yes2DFMM = yes2DFMM;
  inputs.TOL_POW = TOL_POW;
  Vec Phi;
	double start, end;

  DAFMM2D *dafmm2d = new DAFMM2D(inputs);

  int N = dafmm2d->K->N;
  Vec b = Vec::Random(N);//incidence field

  std::cout << "========================= Problem Parameters =========================" << std::endl;
  std::cout << "Matrix Size                        :" << N << std::endl;
  std::cout << "Tolerance                          :" << pow(10,-TOL_POW) << std::endl << std::endl;

  start		=	omp_get_wtime();
	dafmm2d->assemble();
	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);
  std::cout << "========================= Assembly Time =========================" << std::endl;
  std::cout << "Time for assemble in DAFMM form    :" << timeAssemble << std::endl;

	Vec DAFMM_Ab;
  start		=	omp_get_wtime();
	dafmm2d->MatVecProduct(b, DAFMM_Ab);
  end		=	omp_get_wtime();
  double timeMatVecProduct =	(end-start);
  std::cout << "========================= Matrix-Vector Multiplication =========================" << std::endl;
  std::cout << "Time for MatVec in DAFMM form      :" << timeMatVecProduct << std::endl;

	delete dafmm2d;
}

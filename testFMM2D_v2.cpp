#include "FMM2DTree_v2.hpp"

class userkernel{
public:
	userkernel() {
	};
	double ContrastFunction(const pts2D r) {
		double R = r.x*r.x + r.y*r.y;
		double q = 1.5*exp(-160.0*R);
		//double q = 1.5*sin(R);
		return q;
	};
	std::complex<double> IncidenceFunction(const pts2D r) {
		std::complex<double> q = exp(I*kappa*r.x);
		return q;
	};
	std::complex<double> RHSFunction(const pts2D r) {
		std::complex<double> q = -1.0*kappa*kappa*ContrastFunction(r)*IncidenceFunction(r);
		return q;
	};
	~userkernel() {};
};

int main(int argc, char* argv[]) {
	int nCones_LFR			=	atoi(argv[1]);
	int nChebNodes			=	atoi(argv[2]);
	double L						=	atof(argv[3]);
	kappa 							= atof(argv[4]);
	int yes2DFMM				=	atoi(argv[5]);
	int degreeOfBases 	= atoi(argv[6]);
	int treeAdaptivity	=	atoi(argv[7]);

	double TOL_POW = 9;
	cout << "Wavenumber:		" << kappa << endl;
	cout << "Wavelength:		" << 2*PI/kappa << endl;
	cout << "no. of full cycles:	" << L*kappa/2*PI << endl;
	double start, end;

	std::vector<pts2D> particles_X;//locations
	std::vector<pts2D> particles_Y;//dummy values
	std::vector<int> row_indices;
	std::vector<int> col_indices;
	userkernel* mykernel		=	new userkernel();
	FMM2DTree<userkernel>* C	=	new FMM2DTree<userkernel>(mykernel, nCones_LFR, nChebNodes, L, yes2DFMM, TOL_POW, particles_X, particles_Y, row_indices, col_indices, degreeOfBases, treeAdaptivity);
	std::vector<pts2D> gridPoints;
	double timeIn_getMatrixEntry_Offset = 0.0;

	start	=	omp_get_wtime();
	C->set_Standard_Cheb_Nodes();
	C->get_Transfer_Matrix();
	C->createAdaptiveTree();
	C->makeLevelRestriction_Tree();
	C->createCones();
	C->assign_Tree_Interactions();
	C->assignLeafChargeLocations(gridPoints);//chargeLocations
	C->assignNonLeafChargeLocations();
	C->assignLeafChebNodes();

	std::string filename2 = "tree.tex";
	C->outputAdaptiveGrid(filename2);
	std::string filename3 = "M/M";


	C->evaluatePrecomputations();
	C->writeMToFile(filename3);
	C->getNeighborInteractions();


	C->readMFromFile(filename3);
	C->getNeighbors();
	C->readNeighborInteractions();

	C->getNodes_LFR();
	C->getNodes_HFR();

	C->getUserCheckPoints();
	C->Assemble_HFR_M2M();
	C->Assemble_HFR_M2L();
	C->Assemble_LFR_M2L();
	C->Assemble_LFR_L2L();
	C->Assemble_HFR_L2L();

	end		=	omp_get_wtime();
	double timeAssemble =	(end-start);
	std::cout << std::endl << "Time to assemble: " << timeAssemble << std::endl;
	std::cout << "No of matrix entries calculated: " << C->NoIn_getMatrixEntry << std::endl;

	//Vec charges = Vec::Random(C->N);
	Vec charges(C->N);
	for (size_t i = 0; i < C->N; i++) {
		charges(i) = mykernel->RHSFunction(gridPoints[i]); //exp(I*kappa*S->gridPoints[i].x);
	}

	C->assignLeafCharges(charges);
	C->assign_NonLeaf_Charges();

	start		=	omp_get_wtime();
	C->LFR_M2M();
	C->HFR_M2M();
	end		=	omp_get_wtime();
	double timeM2M =	(end-start);
	std::cout << std::endl << "Time taken by M2M: " << timeM2M << std::endl;
	timeIn_getMatrixEntry_Offset = C->timeIn_getMatrixEntry;

	start		=	omp_get_wtime();
	C->LFR_M2L();
	C->HFR_M2L();
	end		=	omp_get_wtime();
	double timeM2L =	(end-start);
	std::cout << std::endl << "Time taken by M2L: " << timeM2L << std::endl;
	timeIn_getMatrixEntry_Offset = C->timeIn_getMatrixEntry;

	start		=	omp_get_wtime();
	C->HFR_L2L();
	C->LFR_L2L();
	end		=	omp_get_wtime();
	double timeL2L =	(end-start);
	std::cout << std::endl << "Time taken by L2L: " << timeL2L << std::endl;
	timeIn_getMatrixEntry_Offset = C->timeIn_getMatrixEntry;

	//C->evaluate_NearField_Using_Precomputations_H();
	start		=	omp_get_wtime();
	C->evaluate_NearField_Using_Precomputations_LS();
	//C->evaluate_NearField();
	end		=	omp_get_wtime();
	double timeNearField =	(end-start);
	std::cout << std::endl << "Time taken by NearField: " << timeNearField << std::endl;

	Vec potential;
	C->collectPotential(potential);
	cout << "potential.n(): " << potential.norm() << endl;

	string filename;
	filename = "result/potential_v2";
	std::ofstream myfile1;
	myfile1.open(filename.c_str());
	myfile1 << potential << endl;
	myfile1.close();

	double totalTime = timeAssemble + timeM2M + timeM2L + timeL2L + timeNearField;
	double MatVecProduct = timeM2M + timeM2L + timeL2L + timeNearField;
	cout << endl << "time for MatVec Product: " << MatVecProduct << endl;

	srand(time(NULL));

	for (size_t i = 0; i < 20; i++) {
		int t	=	rand()%C->leafNodes.size();
		int j = C->leafNodes[t].x;
		int k = C->leafNodes[t].y;
		std::cout << endl << "j: " << j << "	k: " << k << endl;
		double errC = C->perform_Error_Check2(t);
		cout << "		Error is: " << errC << endl;
	}
	std::cout << std::endl;
	delete mykernel;
	delete C;
}

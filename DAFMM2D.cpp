#include "DAFMM2D.hpp"

struct inputsToDFMM {
  int nCones_LFR;
	int nChebNodes;
	double L;
	int yes2DFMM;
	int TOL_POW;
};

class DAFMM2D {
public:
  inputsToDFMM inputs;
  std::vector<pts2D> particles_X;//locations
	std::vector<pts2D> particles_Y;//dummy values
	std::vector<pts2D> gridPoints;//same as particles_X, particles_Y
	std::vector<int> row_indices;
	std::vector<int> col_indices;
	FMM2DTree *K;
  std::string filename;

  DAFMM2D(inputsToDFMM inputs) {
    this->inputs = inputs;
    filename = "DFMMtree.tex";
    K	=	new FMM2DTree(inputs.nCones_LFR, inputs.nChebNodes, inputs.L, inputs.yes2DFMM, inputs.TOL_POW, particles_X, particles_Y, row_indices, col_indices);
    K->set_Standard_Cheb_Nodes();
  	K->get_Transfer_Matrix();
		K->createAdaptiveTree();
		K->makeLevelRestriction_Tree();
		string filename = "tree.tex";
		K->outputAdaptiveGrid(filename);
  	K->createCones();
  	K->assign_Tree_Interactions();
		K->assignLeafChargeLocations(gridPoints);//chargeLocations
    K->assignNonLeafChargeLocations();
    K->assignLeafChebNodes();
  };

  ~DAFMM2D() {};

  // void FMMTreeInitialisation() {
  //   K->set_Standard_Cheb_Nodes();
  // 	K->get_Transfer_Matrix();
  // 	// K->createUniformTree();
	// 	K->createAdaptiveTree();
	// 	K->makeLevelRestriction_Tree();
	// 	string filename = "tree.tex";
	// 	K->outputAdaptiveGrid(filename);
  // 	K->createCones();
  // 	K->assign_Tree_Interactions();
	// 	K->assignLeafChargeLocations(gridPoints);//chargeLocations
  //   K->assignNonLeafChargeLocations();
  //   K->assignLeafChebNodes();
  // }

  void assemble() {
		K->getNodes_LFR();
		K->getNodes_HFR();
		K->getUserCheckPoints();
		K->Assemble_HFR_M2M();
		K->Assemble_HFR_M2L();
		K->Assemble_LFR_M2L();
		K->Assemble_LFR_L2L();
		K->Assemble_HFR_L2L();
		// K->collectBoundary();//to plot
  }

  void MatVecProduct(Vec &b, Vec &output) {
		K->assignLeafCharges(b);
    K->assign_NonLeaf_Charges();
		K->LFR_M2M();
		K->HFR_M2M();
		K->LFR_M2L();
		K->HFR_M2L();
		K->HFR_L2L();
		K->LFR_L2L();
		K->evaluate_NearField();//without precomputations
		K->collectPotential(output);
  }
};

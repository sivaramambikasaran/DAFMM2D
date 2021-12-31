#ifndef _G_FMM2DTree_HPP__
#define _G_FMM2DTree_HPP__
double kappa;
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <chrono>
#include <filesystem>
struct orderedPair {
	int x,y;
};
#include "ACA.hpp"
string currentDirectory;

double arctan(double y, double x) {//returns atan2 in range (0,2*PI)
	double temp = atan2(y, x);
	if (temp < 0.0)
		return temp + 2*PI;
	else
		return temp;
}

const double machinePrecision = pow(10,-16);
using namespace std::chrono;

class FMM2DCone {
public:
	bool outgoing_ILActive;
	bool incoming_ILActive;
	//Directional multipoles and locals of the box in this cone direction
	Vec multipoles;

	//defined only in HFR
	//std::vector<pts2D> &outgoing_chargePoints;//equivalent points {y_{k}^{B,o,l}}
	Vec outgoing_charges;//equivalent densities {f_{k}^{B,o,l}}
	//std::vector<pts2D> &outgoing_checkPoints;//check points {x_{k}^{B,o,l}}
	Vec outgoing_potential;//check potentials {u_{k}^{B,o,l}}

	//outgoing_chargePoints and outgoing_checkPoints are considered to be same as incoming_checkPoints and incoming_chargePoints respectively; hence are not declared
	//This is because of the assumption that the potentials are evaluated at the charge locations
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i,l}}
	Vec incoming_charges;//equivalent densities {f_{k}^{B,i,l}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i,l}}
	std::vector<int> outgoing_chargePoints;//equivalent points {y_{k}^{B,i,l}}
	std::vector<int> outgoing_checkPoints;//check points {x_{k}^{B,i,l}}
	std::vector<int> user_checkPoints;
	Vec incoming_potential;//check potentials {u_{k}^{B,i,l}}

	ColPivHouseholderQR<Mat> outgoing_Atilde_dec;
	ColPivHouseholderQR<Mat> incoming_Atilde_dec;

	std::vector<orderedPair > InteractionList;
	double angle;
	Mat outgoing_Ar;
	Mat L2L[4], M2M, Ktilde;
	std::vector<Mat> M2L;
	FMM2DCone () {}

};

class FMM2DBox {
public:
	double radius;
	int active;
	bool outgoing_ILActive;
	bool incoming_ILActive;
	int level;
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	bool isLeaf;

	std::vector<orderedPair > neighborNumbers;

	int fineNeighbors[12];//12
	int coarseNeighbors[12];//12
	int separatedFineNeighbors[20];//20
	int colleagueNeighbors[9];//9

	std::vector<orderedPair > InteractionList;
	std::vector<int > MFR_IL;
	//defined only in LFR
	//std::vector<pts2D> &outgoing_chargePoints;//equivalent points {y_{k}^{B,o}}
	Vec outgoing_charges;//equivalent densities {f_{k}^{B,o}}
	//std::vector<pts2D> &outgoing_checkPoints;//check points {x_{k}^{B,o}}
	Vec outgoing_potential;//check potentials {u_{k}^{B,o}}

	//outgoing_chargePoints and outgoing_checkPoints are considered to be same as incoming_checkPoints and incoming_chargePoints respectively; hence are not declared
	//This is because of the assumption that the potentials are evaluated at the charge locations
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	Vec incoming_charges;//equivalent densities {f_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	std::vector<int> user_checkPoints;
	std::vector<int> outgoing_chargePoints;//equivalent points {y_{k}^{B,i,l}}
	std::vector<int> outgoing_checkPoints;//check points {x_{k}^{B,i,l}}
	Vec incoming_potential;//check potentials {u_{k}^{B,i}}

	ColPivHouseholderQR<Mat> outgoing_Atilde_dec;
	ColPivHouseholderQR<Mat> incoming_Atilde_dec;

	Mat outgoing_Ar;
	Mat L2L, Ktilde;
	std::vector<Mat> M2L;

	FMM2DBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<4; ++l) {
			childrenNumbers[l]	=	-1;
		}
		isLeaf = false;
		active = true;
		for (size_t i = 0; i < 12; i++) {
			fineNeighbors[i] = -1;
		}
		for (size_t i = 0; i < 12; i++) {
			coarseNeighbors[i] = -1;
		}
		for (size_t i = 0; i < 20; i++) {
			separatedFineNeighbors[i] = -1;
		}
		for (size_t i = 0; i < 9; i++) {
			colleagueNeighbors[i] = -1;
		}
	}

	Vec multipoles;

	pts2D center;
	std::vector<int> chargeLocations;

	std::vector<pts2D> chebNodes;
	std::vector<FMM2DCone> ConeTree;
};

class FMM2DTree: public LowRank {
public:
	double timeIn_getMatrixEntry;//for time profiling
	long NoIn_getMatrixEntry;
	int nLevels;			//	Number of levels in the tree.
	int nChebNodes;			//	Number of Chebyshev nodes along one direction.
	int rank;				//	Rank of interaction, i.e., rank = nChebNodes*nChebNodes.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).
	double a;				//	Cut-off for self-interaction. This is less than the length of the smallest box size.
	const double A = 3.0; //higher A, higher level_FarField
	const double B = 0.35;//4.0; //smaller B, higher level_LFR
	int level_LFR;
	int level_FarField;
	int nCones_LFR;
	Vec chargesAll;
	double SFN_angles[20];
	double N_angles[9];
	double FN_angles[12];
	double CN_angles[12];

	std::vector<pts2D> bottomTopBoundary;//used for plotting field
	std::vector<pts2D> leftRightBoundary;//used for plotting field
	bool findPhi;
	string NeighborFilename;
	string MFilename;
	std::vector<pts2D> gridPoints;
	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<double> ConeAperture;
	std::vector<int> nCones;
	std::vector<double> boxHomogRadius;			//	Stores the value of boxRadius^{alpha}
	std::vector<double> boxLogHomogRadius;		//	Stores the value of alpha*log(boxRadius)
	std::vector<std::vector<FMM2DBox> > tree;	//	The tree storing all the information.
	std::vector<std::vector<int> > indexTree;
	int TOL_POW;
	//	Different Operators
	int yes2DFMM;
	std::vector<double> standardChebNodes1D;
	std::vector<pts2D> standardChebNodes;
	std::vector<pts2D> standardChebNodesChild;
	std::vector<pts2D> leafChebNodes;
	std::vector<orderedPair> leafNodes;
	//	Different Operators
	Eigen::MatrixXd M2M[4];					//	Transfer from multipoles of 4 children to multipoles of parent.
	Eigen::MatrixXd L2L[4];					//	Transfer from locals of parent to locals of 4 children.
	std::vector<Mat > M; //M matrices computed during precomputations
	Eigen::MatrixXd Q_pinv;
	Eigen::MatrixXd Q_pinv2;
	Mat colleagueNeighborInteraction[9][9];
	Mat fineNeighborInteraction[9][12];
	Mat separatedFineNeighborInteraction[9][20];
	Mat coarseNeighborInteraction[9][12];

// public:
	FMM2DTree(int nCones_LFR, int nChebNodes, double L, int yes2DFMM, int TOL_POW, std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y, std::vector<int> row_indices, std::vector<int> col_indices):
	LowRank(TOL_POW, particles_X, particles_Y, row_indices, col_indices) {
		this->findPhi = true;
		this->NoIn_getMatrixEntry = 0;
		this->timeIn_getMatrixEntry = 0.0;
		this->TOL_POW = TOL_POW;
		this->nChebNodes		=	nChebNodes;
		this->rank			=	nChebNodes*nChebNodes;
		this->L				=	L;
		this->nCones_LFR		=	nCones_LFR;
		this->yes2DFMM = yes2DFMM;
		nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		this->a			=	smallestBoxSize;
		int k;
		if (yes2DFMM == 1) {
			if (kappa==0)
				this->level_LFR	=	2;	//actually should be 2; but for checking the accuracy of DFMM i am making it 3; so that HFR code runs even for LFR; if that gives good result it means DFMM code is perfect
			else {
				this->level_LFR	=	floor(log(kappa*L/B)/log(2.0));
			}
			if (level_LFR < 2) level_LFR = 2;
		}
		else {
			level_LFR = 2;
		}
		//cout << "level_LFR: " << level_LFR << endl;
	}

	void shift_scale_Nodes(std::vector<pts2D>& Nodes, std::vector<pts2D>& shifted_scaled_Nodes, double xShift, double yShift, double radius) {
		for (int k=0; k < Nodes.size(); ++k) {
			pts2D temp;
			temp.x	=	radius*Nodes[k].x+xShift;
			temp.y	=	radius*Nodes[k].y+yShift;
			shifted_scaled_Nodes.push_back(temp);
		}
	}

	//	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}

	//	get_S
	double get_S(double x, double y, int n) {
		double S	=	0.5;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return 2.0/n*S;
	}
	//	set_Standard_Cheb_Nodes
	void set_Standard_Cheb_Nodes() {
		for (int k=0; k<nChebNodes; ++k) {
			standardChebNodes1D.push_back(-cos((k+0.5)/nChebNodes*PI));
		}
		pts2D temp1;
		for (int j=0; j<nChebNodes; ++j) {
			for (int k=0; k<nChebNodes; ++k) {
				temp1.x	=	standardChebNodes1D[k];
				temp1.y	=	standardChebNodes1D[j];
				standardChebNodes.push_back(temp1);
			}
		}
		//	Left Bottom child, i.e., Child 0
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Bottom child, i.e., Child 1
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Top child, i.e., Child 2
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Left Top child, i.e., Child 3
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
	}

	void get_Transfer_Matrix() {
		for (int l=0; l<4; ++l) {
			L2L[l]	=	Eigen::MatrixXd(rank,rank);
			for (int j=0; j<rank; ++j) {
				for (int k=0; k<rank; ++k) {
					L2L[l](j,k)	=	get_S(standardChebNodes[k].x, standardChebNodesChild[j+l*rank].x, nChebNodes)*get_S(standardChebNodes[k].y, standardChebNodesChild[j+l*rank].y, nChebNodes);
				}
			}
		}
		for (int l=0; l<4; ++l) {
			M2M[l]	=	L2L[l].transpose();
		}
	}

		void createTree(int j, int k) {
			//k: index of box in the vector tree[j]
			//b: boxNumber of box, tree[j][k]
			bool condition3Leaf = false;
			if (kappa*tree[j][k].radius/2/PI <= 1.0) condition3Leaf = true; //each leaf shoould contain less than or equal to 1 wave cycle
			//nChebnodes are the number of points in 1 wavecycle/Wavelength - in 1D.
			if (condition3Leaf) {
				tree[j][k].isLeaf = true;
				orderedPair op;
				op.x = j;
				op.y = k;
				leafNodes.push_back(op);
			}
			else {
				if (int(tree.size()) == j+1) {
					std::vector<FMM2DBox> level;
					tree.push_back(level);
					std::vector<int> index;
				  indexTree.push_back(index);
				}
				int n	=	tree[j+1].size();
				int b = tree[j][k].boxNumber;
				for (size_t c = 0; c < 4; c++) {
					FMM2DBox box;
					box.level = j+1;
					box.boxNumber		=	b*4+c;
					box.parentNumber	=	b;
					box.radius = 0.5*tree[j][k].radius;
					if (c==0) {
						box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
					}
					else if (c==1) {
						box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
					}
					else if (c==2) {
						box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
					}
					else {
						box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
						box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
					}
					tree[j+1].push_back(box);
					indexTree[j+1].push_back(box.boxNumber);
				}
				for (size_t c = 0; c < 4; c++) {
					createTree(j+1, n+c);
				}
			}
		}

		void createAdaptiveTree() {
			FMM2DBox root;
			root.level = 0;
			root.boxNumber		=	0;
			root.parentNumber	=	-1;
			root.radius = L;
			root.center.x = 0.0;
			root.center.y = 0.0;
			std::vector<FMM2DBox> level;
			level.push_back(root);
			tree.push_back(level);

			std::vector<int> index;
			index.push_back(0);
			indexTree.push_back(index);
			createTree(0, 0);
			nLevels = tree.size() - 1;
			// cout << "nLevels: " << nLevels << endl;
			for (size_t j = 1; j <= 10; j++) {
				boxRadius.push_back(boxRadius[j-1]/2.0);
			}
			if (level_LFR >= nLevels) {
				level_LFR = nLevels-1;
			}
		}

		void getFutureGeneration_LFR(int bj, int bk, int pnj, int pni, std::vector<orderedPair>& futureGeneration) {
			if (tree[pnj][pni].isLeaf) {
				if ( fabs(tree[bj][bk].center.x - tree[pnj][pni].center.x) >= 3.0*boxRadius[bj] + boxRadius[pnj]-machinePrecision
				 || fabs(tree[bj][bk].center.y - tree[pnj][pni].center.y) >= 3.0*boxRadius[bj] + boxRadius[pnj]-machinePrecision) {
				}
				else {
					orderedPair op;
					op.x = pnj;
					op.y = pni;
					futureGeneration.push_back(op);
				}
			}
			else {
				for (int nc=0; nc<4; ++nc) { //children of parents neighbors
					int pnn = tree[pnj][pni].boxNumber;//boxNumber
					int boxB = 4*pnn+nc;//its index=?
					std::vector<int>::iterator indx = std::find(indexTree[pnj+1].begin(), indexTree[pnj+1].end(), boxB);
					int boxB_index = indx-indexTree[pnj+1].begin();
					getFutureGeneration_LFR(bj, bk, pnj+1, boxB_index, futureGeneration);
				}
			}
		}

		void assign_Child_Interaction_LFR(int c, int j, int k, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			int parentboxNumber = tree[j][k].boxNumber;
			int boxA = 4*parentboxNumber+c;//child box number; its index=?
			std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxA);
			int boxA_index = indx-indexTree[j+1].begin();
			for (int n=0; n<Tree_Neighbors_LFR[j][k].size(); ++n) {//parents neighbors; so u need its index to access it which is k
				//children of neighbors of parent which are not neighbors to child=IL
				int pnj = Tree_Neighbors_LFR[j][k][n].x;//level
				int pni = Tree_Neighbors_LFR[j][k][n].y;//index
				int pnn = tree[pnj][pni].boxNumber;//boxNumber

				std::vector<orderedPair> futureGeneration;
				getFutureGeneration_LFR(j+1, boxA_index, pnj, pni, futureGeneration);
				for (size_t d = 0; d < futureGeneration.size(); d++) {
					Tree_Neighbors_LFR[j+1][boxA_index].push_back(futureGeneration[d]);
				}
			}
		}

		//	Assigns the interactions for the children of a box
		void assign_Box_Interactions_LFR(int j, int k, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			if (!tree[j][k].isLeaf) {
				//#pragma omp parallel for
				for (int c=0; c<4; ++c) {
					assign_Child_Interaction_LFR(c,j,k, Tree_Neighbors_LFR);
				}
			}
		}

		//	Assigns the interactions for the children all boxes at a given level
		void assign_Level_Interactions_LFR(int j, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			//#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				assign_Box_Interactions_LFR(j,k, Tree_Neighbors_LFR);//k is index number of box in tree[j] vector
			}
		}

		//	Assigns the interactions for the children all boxes in the tree
		void assign_Tree_Interactions_LFR(std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR) {
			int J = 1;
			for (int c=0; c<4; ++c) {
				for (int n=0; n<4; ++n) {
					orderedPair op;
					op.x = J;
					op.y = n;
					Tree_Neighbors_LFR[J][c].push_back(op);
				}
			}
			for (int j=1; j<=nLevels-1; ++j) {
				assign_Level_Interactions_LFR(j, Tree_Neighbors_LFR);
			}
		}

		void makeLevelRestriction_Tree() {
			for (size_t n = 0; n < nLevels; n++) {
				std::vector<std::vector<std::vector<orderedPair> > > Tree_Neighbors_LFR;
				for (size_t j = 0; j <= nLevels; j++) {
					std::vector<std::vector<orderedPair> > level;
					for (size_t k = 0; k < tree[j].size(); k++) {
						std::vector<orderedPair> box;
						level.push_back(box);
					}
					Tree_Neighbors_LFR.push_back(level);
				}
				assign_Tree_Interactions_LFR(Tree_Neighbors_LFR);
				std::vector<orderedPair> leafNodesOld = leafNodes;
				for (size_t l = 0; l < leafNodesOld.size(); l++) {
					makeLevelRestriction_Box(l, leafNodes[l], Tree_Neighbors_LFR, leafNodesOld);
				}

				int cnt = 0;
 				for (size_t l = 0; l < leafNodesOld.size(); l++) {//removing all the non-leaves from leafNodes vector
					if (leafNodesOld[l].x == -1) {
						leafNodes.erase(leafNodes.begin()+l-cnt);
						cnt++;
					}
				}

				for (size_t j = 0; j <= nLevels; j++) {
					for (size_t k = 0; k < Tree_Neighbors_LFR[j].size(); k++) {
						Tree_Neighbors_LFR[j][k].clear();
					}
					Tree_Neighbors_LFR[j].clear();
				}
				Tree_Neighbors_LFR.clear();
			}
		}

		int getIndexOfInLeafNodes(orderedPair op) {
			for (size_t l = 0; l < leafNodes.size(); l++) {
				if (op.x == leafNodes[l].x && op.y == leafNodes[l].y) {
					return l;
				}
			}
		}

		void makeLevelRestriction_Box(int l, orderedPair leaf, std::vector<std::vector<std::vector<orderedPair> > >& Tree_Neighbors_LFR, std::vector<orderedPair>& leafNodesOld) {
			//l is index of leaf in leafNodes
			/*
			cout << "leafNodesOld.size(): " << leafNodesOld.size() << endl;
			for (size_t r = 0; r < leafNodesOld.size(); r++) {
				cout << r << ", "<< leafNodesOld[r].x << ", " << leafNodesOld[r].y << endl;
			}
			*/
			int j = leaf.x;
			int k = leaf.y;
			for (size_t i = 0; i < Tree_Neighbors_LFR[j][k].size(); i++) {
				int nj = Tree_Neighbors_LFR[j][k][i].x;
				int nk = Tree_Neighbors_LFR[j][k][i].y;
				std::vector<orderedPair> newLeaves;

				if (j-nj >= 2 || nj-j >= 2) {
					if (j-nj >= 2) {
						int indexOfBox = getIndexOfInLeafNodes(Tree_Neighbors_LFR[j][k][i]);//index of Tree_Neighbors_LFR[j][k][i]
						if (leafNodesOld[indexOfBox].x != -1) {
							refineBox(indexOfBox, Tree_Neighbors_LFR[j][k][i], leafNodesOld);
						}
					}
					else if(nj-j >= 2) {
						if (leafNodesOld[l].x != -1) {
							refineBox(l, leaf, leafNodesOld);
						}
					}
				}
			}
		}

		void refineBox(int l, orderedPair leaf, std::vector<orderedPair> &leafNodesOld) {
			int j = leaf.x; //level
			int k = leaf.y; //box index
			int b = tree[j][k].boxNumber;
			leafNodesOld[l].x = -1;//making it a non-leaf node
			tree[j][k].isLeaf = false;
			for (size_t c = 0; c < 4; c++) {
				FMM2DBox box;
				box.isLeaf = true;
				box.level = j+1;
				box.boxNumber		=	b*4+c;
				box.parentNumber	=	b;
				box.radius = 0.5*tree[j][k].radius;
				if (c==0) {
					box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
				}
				else if (c==1) {
					box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y-0.5*tree[j][k].radius;
				}
				else if (c==2) {
					box.center.x		=	tree[j][k].center.x+0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
				}
				else {
					box.center.x		=	tree[j][k].center.x-0.5*tree[j][k].radius;
					box.center.y		=	tree[j][k].center.y+0.5*tree[j][k].radius;
				}
				orderedPair op;
				op.x = j+1;
				op.y = tree[j+1].size();
				leafNodes.push_back(op);
				tree[j+1].push_back(box);
				indexTree[j+1].push_back(box.boxNumber);
			}
		}

		void assignLeafCharges() {
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				tree[j][b].multipoles	=	0.5*(Vec::Ones(rank));//+Eigen::VectorXd::Random(rank));
			}
		}

		void assignLeafCharges(Vec &charges) {
			chargesAll = charges;
			int start = 0;
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				tree[j][b].multipoles	=	charges.segment(start, rank);
				start += rank;
			}
		}

		void assignLeafChargeLocations(std::vector<pts2D> &particles_out) {
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				int startIndex = gridPoints.size();
				for (size_t i = 0; i < rank; i++) {
					tree[j][b].chargeLocations.push_back(startIndex+i);
				}
				shift_scale_Nodes(standardChebNodes, gridPoints, tree[j][b].center.x, tree[j][b].center.y, boxRadius[j]);
			}
			particles_X = gridPoints;//object of base class FMM_Matrix
			particles_Y = gridPoints;
			particles_out = gridPoints;
		}

		void collectBoundary() {
			for (size_t t = 0; t < leafNodes.size(); t++) {
				int j = leafNodes[t].x;
				int k = leafNodes[t].y;
				pts2D leftRightBoundaryTemp, bottomTopBoundaryTemp;
				leftRightBoundaryTemp.x = tree[j][k].center.x - tree[j][k].radius;//left
				leftRightBoundaryTemp.y = tree[j][k].center.x + tree[j][k].radius;//right
				bottomTopBoundaryTemp.x = tree[j][k].center.y - tree[j][k].radius;//bottom
				bottomTopBoundaryTemp.y = tree[j][k].center.y + tree[j][k].radius;//top
				leftRightBoundary.push_back(leftRightBoundaryTemp);
				bottomTopBoundary.push_back(bottomTopBoundaryTemp);
			}
		}

		void assignNonLeafChargeLocations() {
			for (int j=nLevels-1; j>1; --j) {
				int J	=	j+1;
				#pragma omp parallel for
				for (int k=0; k<tree[j].size(); ++k) {
					if (!tree[j][k].isLeaf) {
						int b = tree[j][k].boxNumber;
						int KboxNumber;
						int K[4];
						std::vector<int>::iterator indx;
						KboxNumber = 4*b+0;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[0] = indx-indexTree[J].begin();
						KboxNumber = 4*b+1;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[1] = indx-indexTree[J].begin();
						KboxNumber = 4*b+2;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[2] = indx-indexTree[J].begin();
						KboxNumber = 4*b+3;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[3] = indx-indexTree[J].begin();

						for (int c=0; c<4; ++c) {
							//Block containing n elements, starting at position i: vector.segment(i,n)
							for (int i = 0; i < tree[J][K[c]].chargeLocations.size(); i++) {
								tree[j][k].chargeLocations.push_back(tree[J][K[c]].chargeLocations[i]);
							}
						}
					}
				}
			}
		}

		void assignLeafChebNodes() {
			for (size_t k = 0; k < leafNodes.size(); k++) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				shift_scale_Nodes(standardChebNodes, tree[j][b].chebNodes, tree[j][b].center.x, tree[j][b].center.y, boxRadius[j]);
			}
		}

		double errInApproximation(Eigen::VectorXd& trueValue, Eigen::VectorXd& approximateValue) {
			Eigen::VectorXd error(trueValue.size());
			for (int k=0; k<trueValue.size(); ++k) {
				error(k)	=	fabs(trueValue(k)-approximateValue(k));
				//error(k)	=	fabs((trueValue-approximateValue)(k)/trueValue(k));
			}
			return error.maxCoeff();///trueValue.maxCoeff();
		}

		double errInApproximation(Vec& trueValue, Vec& approximateValue) {
			Vec error;
			error = trueValue - approximateValue;
			VectorXd errorAbs = error.cwiseAbs();
			VectorXd trueValueAbs = trueValue.cwiseAbs();
			return errorAbs.maxCoeff();///trueValueAbs.maxCoeff();
		}

		void outputAdaptiveGrid(std::string filename) {
			double xcenter = 0.0;
			double ycenter = 0.0;
			double Lx = L;
			double Ly = L;

			std::ofstream myfile;
			myfile.open(filename.c_str());
			myfile << "\\documentclass{standalone}" << std::endl;
			myfile << "\\usepackage{tikz}" << std::endl;
			myfile << "\\begin{document}" << std::endl;
			myfile << "\\begin{tikzpicture}" << std::endl;
			for (int k=0; k<(int)leafNodes.size(); ++k) {
				int j = leafNodes[k].x;
				int b = leafNodes[k].y;
				myfile << "\\draw (" << tree[j][b].center.x-tree[j][b].radius << ",";
				myfile << tree[j][b].center.y-tree[j][b].radius << ") rectangle (";
				//myfile << tree[j][b].center.y-tree[j][b].radius << ") rectangle node{\\tiny " << b << "} (";
				//myfile << tree[j][b].center.y-tree[j][b].radius << ") rectangle node{\\tiny " << tree[j][b].boxNumber << "} (";
				myfile << tree[j][b].center.x+tree[j][b].radius << ",";
				myfile << tree[j][b].center.y+tree[j][b].radius << ");" << std::endl;
			}
			long double push	=	0.125;
			myfile<< "\\node at (" << xcenter-Lx-push << "," << ycenter-Ly-push << ") {\\tiny$(" << xcenter-Lx << "," << ycenter-Ly << ")$};" << std::endl;
			myfile<< "\\node at (" << xcenter-Lx-push << "," << ycenter+Ly+push << ") {\\tiny$(" << xcenter-Lx << "," << ycenter+Ly << ")$};" << std::endl;
			myfile<< "\\node at (" << xcenter+Lx+push << "," << ycenter-Ly-push << ") {\\tiny$(" << xcenter+Lx << "," << ycenter-Ly << ")$};" << std::endl;
			myfile<< "\\node at (" << xcenter+Lx+push << "," << ycenter+Ly+push << ") {\\tiny$(" << xcenter+Lx << "," << ycenter+Ly << ")$};" << std::endl;
			myfile << "\\end{tikzpicture}" << std::endl;
			myfile << "\\end{document}" << std::endl;
			myfile.close();
		}
		//////////////////////////////////////////////////////////////

	double find_distance(int j1, int k1, int j2, int k2) {
		pts2D r1 = tree[j1][k1].center;
		pts2D r2 = tree[j2][k2].center;
		return	sqrt((r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y));
	}

	void assign_Child_Interaction(int c, int j, int k) {
		// j: level
		// k: parent box
		int parentboxNumber = tree[j][k].boxNumber;
		int boxA = 4*parentboxNumber+c;//child box number; its index=?
		std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxA);
		int boxA_index = indx-indexTree[j+1].begin();
		for (int n=0; n<tree[j][k].neighborNumbers.size(); ++n) {//parents neighbors; so u need its index to access it which is k
			//children of neighbors of parent which are not neighbors to child=IL
			int pnj = tree[j][k].neighborNumbers[n].x;//level
			int pni = tree[j][k].neighborNumbers[n].y;//index
			int pnn = tree[pnj][pni].boxNumber;//boxNumber
			if (j+1 >= level_LFR || (j+1 < level_LFR && tree[j+1][boxA_index].isLeaf)) { //LFR
				if (tree[pnj][pni].isLeaf) {
					if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj][pni].center.x) >= 3.0*boxRadius[j+1] + boxRadius[pnj]-machinePrecision
					 || fabs(tree[j+1][boxA_index].center.y - tree[pnj][pni].center.y) >= 3.0*boxRadius[j+1] + boxRadius[pnj]-machinePrecision) {
						 orderedPair op;
						 op.x = pnj;
						 op.y = pni;
						 tree[j+1][boxA_index].InteractionList.push_back(op);
					}
					else {
						orderedPair op;
						op.x = pnj;
						op.y = pni;
						tree[j+1][boxA_index].neighborNumbers.push_back(op);
					}
				}
				else {
					for (int nc=0; nc<4; ++nc) { //children of parents neighbors
						int boxB = 4*pnn+nc;//its index=?
						std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxB);
						int boxB_index = indx-indexTree[j+1].begin();
						if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj+1][boxB_index].center.x) >= 3.0*boxRadius[j+1] + boxRadius[pnj+1]-machinePrecision
						 || fabs(tree[j+1][boxA_index].center.y - tree[pnj+1][boxB_index].center.y) >= 3.0*boxRadius[j+1] + boxRadius[pnj+1]-machinePrecision) {
							 orderedPair op;
							 op.x = pnj+1;
							 op.y = boxB_index;
							 tree[j+1][boxA_index].InteractionList.push_back(op);
						}
						else {
							if (!tree[j+1][boxA_index].isLeaf) {
								orderedPair op;
								op.x = pnj+1;
								op.y = boxB_index;
								tree[j+1][boxA_index].neighborNumbers.push_back(op);
							}
							else {
								int pj = pnj+1;
								int pk = boxB_index;

								std::vector<orderedPair> futureGeneration;
								getFutureGeneration_LFR(j+1, boxA_index, pj, pk, futureGeneration);
								for (size_t d = 0; d < futureGeneration.size(); d++) {
									tree[j+1][boxA_index].neighborNumbers.push_back(futureGeneration[d]);
								}
							}
						}
					}
				}
			}
			else { //HFR
				if (tree[pnj][pni].isLeaf) {
					if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj][pni].center.x) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[j+1] + boxRadius[pnj]-machinePrecision
					 || fabs(tree[j+1][boxA_index].center.y - tree[pnj][pni].center.y) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[j+1] + boxRadius[pnj]-machinePrecision) {
							double arg = atan2(tree[j+1][boxA_index].center.y-tree[pnj][pni].center.y, tree[j+1][boxA_index].center.x-tree[pnj][pni].center.x);
	 						arg = fmod(arg+2*PI+PI, 2*PI);
	 						int coneNum = int(arg/ConeAperture[j+1]);
							orderedPair op;
							op.x = pnj;
							op.y = pni;
							tree[j+1][boxA_index].ConeTree[coneNum].InteractionList.push_back(op);
					}
					else {
						orderedPair op;
						op.x = pnj;
						op.y = pni;
						tree[j+1][boxA_index].neighborNumbers.push_back(op);
					}
				}
				else {
					for (int nc=0; nc<4; ++nc) { //children of parents neighbors
						int boxB = 4*pnn+nc;//its index=?
						std::vector<int>::iterator indx = std::find(indexTree[j+1].begin(), indexTree[j+1].end(), boxB);
						int boxB_index = indx-indexTree[j+1].begin();
						if ( fabs(tree[j+1][boxA_index].center.x - tree[pnj+1][boxB_index].center.x) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[j+1] + boxRadius[pnj+1]-machinePrecision
						 || fabs(tree[j+1][boxA_index].center.y - tree[pnj+1][boxB_index].center.y) >= kappa*boxRadius[j+1]*boxRadius[j+1] + boxRadius[pnj+1]-machinePrecision) {
							 double arg = atan2(tree[j+1][boxA_index].center.y-tree[pnj+1][boxB_index].center.y, tree[j+1][boxA_index].center.x-tree[pnj+1][boxB_index].center.x);
	 	 						arg = fmod(arg+2*PI+PI, 2*PI);
	 	 						int coneNum = int(arg/ConeAperture[j+1]);
	 							orderedPair op;
	 							op.x = pnj+1;
	 							op.y = boxB_index;
	 							tree[j+1][boxA_index].ConeTree[coneNum].InteractionList.push_back(op);
						}
						else {
							orderedPair op;
							op.x = pnj+1;
							op.y = boxB_index;
							tree[j+1][boxA_index].neighborNumbers.push_back(op);
						}
					}
				}
			}
		}
	}

	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		if (!tree[j][k].isLeaf) {
			#pragma omp parallel for
			for (int c=0; c<4; ++c) {
				assign_Child_Interaction(c,j,k);
			}
		}
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		#pragma omp parallel for
		for (int k=0; k<tree[j].size(); ++k) {
			assign_Box_Interactions(j,k);//k is index number of box in tree[j] vector
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions() {
		//j=0, no neighbors, no IL
		//j=1, no IL, neighbors yes
		//neighbor includes self
		int j = 1;
		for (int c=0; c<4; ++c) {
			for (int n=0; n<4; ++n) {
				orderedPair op;
				op.x = 1;
				op.y = n;
				tree[j][c].neighborNumbers.push_back(op);
			}
		}
		for (j=1; j<=nLevels-1; ++j) {
			assign_Level_Interactions(j);
		}
	}

	void createCones() {
		N = leafNodes.size()*rank;
		// cout << "level_LFR: " << level_LFR << endl;
		// cout << "Number of particles: " << N << endl;
		nCones.push_back(nCones_LFR*pow(2.0, level_LFR-1));
		ConeAperture.push_back(2*PI/nCones[0]);
		for (size_t j = 1; j <= level_LFR-1; j++) {
			nCones.push_back(nCones[j-1]/2.0);
			ConeAperture.push_back(ConeAperture[j-1]*2.0);
		}
		for (size_t j = 0; j <= level_LFR-1; j++) {
			for (int k=0; k<tree[j].size(); ++k) {
				if (!tree[j][k].isLeaf) {
					tree[j][k].ConeTree.clear();////
					for (int c=0; c<nCones[j]; ++c) {
						FMM2DCone cone;
						cone.angle	=	ConeAperture[j]/2.0 + c*ConeAperture[j];
						tree[j][k].ConeTree.push_back(cone);
					}
				}
			}
		}
	}

	void assign_NonLeaf_Charges() {
		for (int j=nLevels-1; j>1; --j) {
			int J	=	j+1;
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (!tree[j][k].isLeaf) {
					int b = tree[j][k].boxNumber;
					int KboxNumber;
					int K[4];
					std::vector<int>::iterator indx;
					KboxNumber = 4*b+0;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[0] = indx-indexTree[J].begin();
					KboxNumber = 4*b+1;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[1] = indx-indexTree[J].begin();
					KboxNumber = 4*b+2;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[2] = indx-indexTree[J].begin();
					KboxNumber = 4*b+3;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[3] = indx-indexTree[J].begin();

					int NumCharges = tree[J][K[0]].chargeLocations.size() + tree[J][K[1]].chargeLocations.size() +tree[J][K[2]].chargeLocations.size() + tree[J][K[3]].chargeLocations.size();
					tree[j][k].multipoles = Vec::Zero(NumCharges);
					int start = 0;
					for (int c=0; c<4; ++c) {
						//Block containing n elements, starting at position i: vector.segment(i,n)
						int NumElem = tree[J][K[c]].chargeLocations.size();
						tree[j][k].multipoles.segment(start, NumElem) = tree[J][K[c]].multipoles;
						start += NumElem;
					}
				}
			}
		}
	}

	void getParticlesFromChildrenLFR_outgoing_col(int j, int k, std::vector<int>& searchNodes) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					// if (tree[J][K].incoming_checkPoints.size() == 0) {
					// 	cout << "problem LFR_outgoing_row: " << j << ", " << k << ", " << J << ", " << K << endl;
					// }
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone].outgoing_chargePoints.end());
					}
				}
			}
		}
	}

	void getParticlesFromChildrenLFR_outgoing_row(int j, int k, std::vector<int>& searchNodes, int j1, int k1) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone].incoming_checkPoints.end());
					}
				}
			}
		}
	}

	void getNodes_LFR() {
		for (int j=nLevels; j>=2; j--) {
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].incoming_chargePoints.size() > 0)
					tree[j][k].incoming_chargePoints.clear();
				if (tree[j][k].incoming_checkPoints.size() > 0)
					tree[j][k].incoming_checkPoints.clear();
				if (tree[j][k].outgoing_chargePoints.size() > 0)
					tree[j][k].outgoing_chargePoints.clear();
				if (tree[j][k].outgoing_checkPoints.size() > 0)
					tree[j][k].outgoing_checkPoints.clear();
				if (tree[j][k].user_checkPoints.size() > 0)
					tree[j][k].user_checkPoints.clear();
			}
		}
		for (int j=nLevels; j>=level_LFR; j--) {
			getNodes_LFR_outgoing_level(j);
			getNodes_LFR_incoming_level(j);
		}
	}

	void getNodes_LFR_outgoing_box(int j, int k, int &n_rows, int &n_cols, int &ComputedRank) {
		int ILcheck = 0;
		if (tree[j][k].active == true) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenLFR_outgoing_col(j, k, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;//indices
			for (int l=0; l<tree[j][k].InteractionList.size(); ++l) {
				int jIL = tree[j][k].InteractionList[l].x;
				int kIL = tree[j][k].InteractionList[l].y;
				std::vector<int> chargeLocations;
				getParticlesFromChildrenLFR_outgoing_row(jIL, kIL, chargeLocations, j, k);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = IL_Nodes.size();
			n_cols = boxA_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = IL_Nodes;
			col_indices = boxA_Nodes;//object of base class G_LowRank
			std::vector<int> row_bases, col_bases;
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].outgoing_ILActive = true;
				Mat dummy;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy, tree[j][k].outgoing_Ar);
				int minN = n_rows;
				if (n_rows > n_cols) {
					minN = n_cols;
				}
				if(ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].outgoing_checkPoints.push_back(IL_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].outgoing_chargePoints.push_back(boxA_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(IL_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(boxA_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].outgoing_Atilde_dec = Atilde.colPivHouseholderQr();
				}
			}
			// if (n_rows == 0) {
			else {
				tree[j][k].outgoing_ILActive = false;
				getParticlesFromChildrenLFR_outgoing_col(j, k, tree[j][k].outgoing_chargePoints);
			}
		}
		else {
			tree[j][k].outgoing_ILActive = false;
		}
	}

	void getNodes_LFR_outgoing_level(int j) { //LFR; box interactions
	  /*
	  computes
	  1. x^{B,o}
	  2. y^{B,o}
	  3. matrix decomposition: K(x^{B,o}, y^{B,o})
	  */
	  //for (int j=nLevels; j>=2; j--) {
	    int rankPerLevel = 0;
	    int n_rows_checkpoint;
	    int n_cols_checkpoint;
	    std::vector<int> boxA_Particles_checkpoint;
	    std::vector<int> IL_Particles_checkpoint;
			int ComputedRank, n_rows, n_cols;
			//#pragma omp parallel for
	    for (int k=0; k<tree[j].size(); ++k) {
				getNodes_LFR_outgoing_box(j, k, n_rows, n_cols, ComputedRank);
				if (rankPerLevel < ComputedRank) {
					rankPerLevel = ComputedRank;
					n_rows_checkpoint = n_rows;
					n_cols_checkpoint = n_cols;
				}
	    }
	  //}
		cout << "O;	j: " << j << "	Nboxes: " << tree[j].size() << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
	}


	void getParticlesFromChildrenLFR_incoming_row(int j, int k, std::vector<int>& searchNodes) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					// if (tree[J][K].incoming_checkPoints.size() == 0) {
					// 	cout << "problem LFR_outgoing_row: " << j << ", " << k << ", " << J << ", " << K << endl;
					// }
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone].incoming_checkPoints.end());
					}
				}
			}
		}
	}

	void getParticlesFromChildrenLFR_incoming_col(int j, int k, std::vector<int>& searchNodes) {
		if (tree[j][k].isLeaf) {
			searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			for (int c = 0; c < 4; c++) {
				int KboxNumber = 4*b+c;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (J >= level_LFR || (J<level_LFR && tree[J][K].isLeaf)) { //LFR
					// if (tree[J][K].incoming_checkPoints.size() == 0) {
					// 	cout << "problem LFR_outgoing_row: " << j << ", " << k << ", " << J << ", " << K << endl;
					// }
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());
				}
				else { //HFR
					for (size_t cone = 0; cone < nCones[J]; cone++) {
						searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone].outgoing_chargePoints.end());
					}
				}
			}
		}
	}

	void getNodes_LFR_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
		int ILcheck = 0;
		if (tree[j][k].active == true) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenLFR_incoming_row(j, k, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;//indices
			for (int l=0; l<tree[j][k].InteractionList.size(); ++l) {
				int jIL = tree[j][k].InteractionList[l].x;
				int kIL = tree[j][k].InteractionList[l].y;
				std::vector<int> chargeLocations;
				getParticlesFromChildrenLFR_incoming_col(jIL, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = boxA_Nodes.size();
			n_cols = IL_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = boxA_Nodes;//object of base class G_LowRank
			col_indices = IL_Nodes;
			std::vector<int> row_bases, col_bases;
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].incoming_ILActive = true;
				Mat dummy1, dummy2;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy1, dummy2);
				int minN = n_rows;
				if (n_rows > n_cols) {
					minN = n_cols;
				}
				if(ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(boxA_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(IL_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].incoming_Atilde_dec = Atilde.colPivHouseholderQr();
				}
			}
			// if (n_cols == 0) {
			else {
				tree[j][k].incoming_ILActive = false;
				getParticlesFromChildrenLFR_incoming_row(j, k, tree[j][k].incoming_checkPoints);
			}
		}
		else {
			tree[j][k].incoming_ILActive = false;
		}
	}

	void getNodes_LFR_incoming_level(int j) { //LFR; box interactions
		/*
		computes
		1. x^{B,o}
		2. y^{B,o}
		3. matrix decomposition: K(x^{B,o}, y^{B,o})
		*/
		//for (int j=nLevels; j>=2; j--) {
			int rankPerLevel = 0;
			int n_rows_checkpoint;
			int n_cols_checkpoint;
			std::vector<int> boxA_Particles_checkpoint;
			std::vector<int> IL_Particles_checkpoint;
			int ComputedRank, n_rows, n_cols;
			//#pragma omp parallel for
	    for (int k=0; k<tree[j].size(); ++k) {
				getNodes_LFR_incoming_box(j, k, n_rows, n_cols, ComputedRank);
				if (rankPerLevel < ComputedRank) {
					rankPerLevel = ComputedRank;
					n_rows_checkpoint = n_rows;
					n_cols_checkpoint = n_cols;
				}
			}
		//}
		cout << "I;	j: " << j << "	Nboxes: " << tree[j].size() << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
	}

	void LFR_M2M_ILActive_True(int j, int k) {
		std::vector<int> source_points;// = tree[j][k].chargeLocations//source points
		Vec source_densities;
		if (tree[j][k].isLeaf) {
			int Veclength = tree[j][k].multipoles.size();
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int l = 0; l < tree[j][k].multipoles.size(); l++) {
				source_densities(g) = tree[j][k].multipoles(l);
				++g;
			}
			source_points.insert(source_points.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			int KboxNumber;
			int K[4];
			std::vector<int>::iterator indx;
			KboxNumber = 4*b+0;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[0] = indx-indexTree[J].begin();
			KboxNumber = 4*b+1;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[1] = indx-indexTree[J].begin();
			KboxNumber = 4*b+2;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[2] = indx-indexTree[J].begin();
			KboxNumber = 4*b+3;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[3] = indx-indexTree[J].begin();

			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					source_densities(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
				source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		tree[j][k].outgoing_potential = tree[j][k].outgoing_Ar*source_densities;//u^{B,o}
		tree[j][k].outgoing_charges = tree[j][k].outgoing_Atilde_dec.solve(tree[j][k].outgoing_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
		// if (tree[j][k].outgoing_charges.norm() > pow(10,10.0)) {
		// 	cout << "j: " << j << "	k: " << k << "	o_C.n: " << tree[j][k].outgoing_charges.norm() << endl;
		// }
	}

	void LFR_M2M_ILActive_False(int j, int k) {
		if (tree[j][k].isLeaf) {
			tree[j][k].outgoing_charges = tree[j][k].multipoles;
		}
		else {
			int J = j+1;
			int b = tree[j][k].boxNumber;
			int KboxNumber;
			int K[4];
			std::vector<int>::iterator indx;
			KboxNumber = 4*b+0;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[0] = indx-indexTree[J].begin();
			KboxNumber = 4*b+1;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[1] = indx-indexTree[J].begin();
			KboxNumber = 4*b+2;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[2] = indx-indexTree[J].begin();
			KboxNumber = 4*b+3;
			indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
			K[3] = indx-indexTree[J].begin();

			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			tree[j][k].outgoing_charges = Vec::Zero(Veclength);
			int g = 0;
			for (int child = 0; child < 4; child++) {
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					tree[j][k].outgoing_charges(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
			}
		}
	}

	void LFR_M2M() {//outgoing operations
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
	  x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		A = K(x^{B,o}, y^{B,o})
		x^{B,o}=tree[j][k].outgoing_checkPoints
		y^{B,o}=tree[j][k].outgoing_chargePoints
		*/
		for (int j=nLevels; j>=2; --j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (tree[j][k].active == true) {
					if (tree[j][k].outgoing_ILActive) {
						if (tree[j][k].outgoing_chargePoints.size() > 0) {
							LFR_M2M_ILActive_True(j, k);
						}
					}
					else {
						LFR_M2M_ILActive_False(j, k);
					}
				}
			}
		}
	}

	void Assemble_LFR_M2L() {
		/*
		u^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o}: tree[j][boxA].ConeTree[coneB].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<=nLevels; ++j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {//BoxA
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (tree[j][k].active == true) {
					if (!tree[j][k].isLeaf) {
						//tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].user_checkPoints.size());
						if (tree[j][k].incoming_ILActive) {
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;
								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].M2L.push_back(getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints));
										//tree[j][k].incoming_potential += R*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									// }
								}
								else {
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									tree[j][k].M2L.push_back(getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints));
									//tree[j][k].incoming_potential += R*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
					else {
						//tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].chargeLocations.size());
						if (tree[j][k].incoming_ILActive) {
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;
								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();

									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].M2L.push_back(getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints));
										//tree[j][k].incoming_potential += R*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									// }
								}
								else {
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									tree[j][k].M2L.push_back(getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].outgoing_chargePoints));
									//tree[j][k].incoming_potential += R*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
				}
			}
		}
	}

	void LFR_M2L() {
		/*
		u^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o}: tree[j][boxA].ConeTree[coneB].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<=nLevels; ++j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {//BoxA
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (tree[j][k].active == true) {
					if (!tree[j][k].isLeaf) {
						tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].user_checkPoints.size());
						if (tree[j][k].incoming_ILActive) {
							#pragma omp parallel for
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;
								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										//Mat R = getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints);
										tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									}
								}
								else {
									int n_rows = tree[j][k].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									//Mat R = getMatrix(tree[j][k].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints);
									tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
					else {
						tree[j][k].incoming_potential	=	Vec::Zero(tree[j][k].chargeLocations.size());
						if (tree[j][k].incoming_ILActive) {
							#pragma omp parallel for
							for (int l = 0; l < tree[j][k].InteractionList.size(); l++) {
								int jIL = tree[j][k].InteractionList[l].x;
								int kIL = tree[j][k].InteractionList[l].y;
								if (jIL<level_LFR && !tree[jIL][kIL].isLeaf) {//HFR
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										//Mat R = getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints);
										tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									}
								}
								else {
									int n_rows = tree[j][k].chargeLocations.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									//Mat R = getMatrix(tree[j][k].chargeLocations, tree[jIL][kIL].outgoing_chargePoints);
									tree[j][k].incoming_potential += tree[j][k].M2L[l]*tree[jIL][kIL].outgoing_charges;
								}
							}
						}
					}
				}
			}
		}
	}

	void Assemble_LFR_L2L_IL_True(int j, int k) {
		//tree[j][k].multipoles//source densities
		//tree[j][k].chargeLocations//source points
		//x^{B,o}=tree[j][k].outgoing_checkPoints
		//kernel evaluation between x^{B,o} and source points
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();

		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			//x^{C,i}: tree[J][4*k+c].incoming_checkPoints
			//f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
			//y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
			Mat R;
			if (!tree[J][K[child]].isLeaf) {
				tree[J][K[child]].L2L = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].incoming_chargePoints);
			}
			else {
				if (j >= level_LFR) {//LFR
					tree[J][K[child]].L2L = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].incoming_chargePoints);
				}
				else {
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[J][K[child]].L2L = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						}
					}
				}
			}
		}
	}

	void LFR_L2L_IL_True(int j, int k) {
		//tree[j][k].multipoles//source densities
		//tree[j][k].chargeLocations//source points
		//x^{B,o}=tree[j][k].outgoing_checkPoints
		//kernel evaluation between x^{B,o} and source points
		tree[j][k].incoming_charges = tree[j][k].incoming_Atilde_dec.solve(tree[j][k].incoming_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();

		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			//x^{C,i}: tree[J][4*k+c].incoming_checkPoints
			//f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
			//y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
			Mat R;
			if (!tree[J][K[child]].isLeaf) {
				tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].incoming_charges;//u^{B,o}
			}
			else {
				if (j >= level_LFR) {//LFR
					tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].incoming_charges;//u^{B,o}
				}
				else {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
				}
			}
		}
	}

	void LFR_L2L_IL_False(int j, int k) {
		int J = j+1;
		int b = tree[j][k].boxNumber;

		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		int offset = 0;

		for (int child = 0; child < 4; child++) {
			if (!tree[J][K[child]].isLeaf) {
				for (size_t i = 0; i < tree[J][K[child]].user_checkPoints.size(); i++) {
					tree[J][K[child]].incoming_potential(i) += tree[j][k].incoming_potential(offset+i);
				}
				offset += tree[J][K[child]].user_checkPoints.size();
			}
			else {
				if (j >= level_LFR) {//LFR
					for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
						tree[J][K[child]].incoming_potential(i) += tree[j][k].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].chargeLocations.size();
				}
				else {
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
							tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
						}
					}
					offset += tree[J][K[child]].chargeLocations.size();
				}
			}
		}
	}

	void LFR_L2L() {//outgoing operations
		for (int j=2; j<nLevels; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (j<level_LFR && !tree[j][k].isLeaf) {
					continue;
				}
				if (!tree[j][k].isLeaf) {
					if (tree[j][k].active == true) {
						if (tree[j][k].incoming_ILActive) {
							if (tree[j][k].incoming_chargePoints.size() > 0) {
								LFR_L2L_IL_True(j, k);
							}
						}
						else {
							LFR_L2L_IL_False(j, k);
						}
					}
				}
			}
		}
	}


		void Assemble_LFR_L2L() {//outgoing operations
			#pragma omp parallel for
			for (int j=2; j<nLevels; ++j) {//parent
				#pragma omp parallel for
				for (int k=0; k<tree[j].size(); ++k) {
					if (j<level_LFR && !tree[j][k].isLeaf) {
						continue;
					}
					if (!tree[j][k].isLeaf) {
						if (tree[j][k].active == true) {
							if (tree[j][k].incoming_ILActive) {
								Assemble_LFR_L2L_IL_True(j, k);
							}
						}
					}
				}
			}
		}

	void evaluate_NearField_Using_Precomputations_H() {
		for (int t=0; t<leafNodes.size(); ++t) {
			int j = leafNodes[t].x;
			int k = leafNodes[t].y;
			if (tree[j][k].active == true) {
				for (size_t nColleagues = 0; nColleagues < 9; nColleagues++) {
					int nj = j;
					int nk = tree[j][k].colleagueNeighbors[nColleagues];
					if (nk != -1) {
						Mat boxOperator = colleagueNeighborInteraction[j-2][nColleagues];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				for (size_t nFine = 0; nFine < 12; nFine++) {
					int nj = j+1;
					int nk = tree[j][k].fineNeighbors[nFine];
					if (nk != -1) {
						Mat boxOperator = fineNeighborInteraction[j-2][nFine];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				for (size_t nSFine = 0; nSFine < 20; nSFine++) {
					int nj = j+1;
					int nk = tree[j][k].separatedFineNeighbors[nSFine];
					if (nk != -1) {
						Mat boxOperator = separatedFineNeighborInteraction[j-2][nSFine];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
				for (size_t nCoarse = 0; nCoarse < 12; nCoarse++) {
					int nj = j-1;
					int nk = tree[j][k].coarseNeighbors[nCoarse];
					if (nk != -1) {
						Mat boxOperator = coarseNeighborInteraction[j-2][nCoarse];
						tree[j][k].incoming_potential += boxOperator*tree[nj][nk].multipoles;
					}
				}
			}
		}
	}

	void evaluate_NearField() {
	//we are always making sure that we have a tree high enough to be in LFR;
	//so, near filed needs to be done only in LFR
		//#pragma omp parallel for
		for (int t=0; t<leafNodes.size(); ++t) {
			int j = leafNodes[t].x;
			int k = leafNodes[t].y;
			if (tree[j][k].active == true) {
				//Neighbor Interaction
				for (int l = 0; l < tree[j][k].neighborNumbers.size(); l++) {//excluding self
					int nj = tree[j][k].neighborNumbers[l].x;
					int nk = tree[j][k].neighborNumbers[l].y;
					int n_rows = tree[j][k].chargeLocations.size();
					int n_cols = tree[nj][nk].chargeLocations.size();
					Mat R = getMatrix(tree[j][k].chargeLocations, tree[nj][nk].chargeLocations);
					tree[j][k].incoming_potential += R*tree[nj][nk].multipoles;//u^{B,o}
				}
			}
		}
	}

	void collectPotential(Vec &potential) {
		potential = VectorXcd::Zero(N);
		int start = 0;
		for (size_t t = 0; t < leafNodes.size(); t++) {
			int j = leafNodes[t].x;
			int k = leafNodes[t].y;
			potential.segment(start, rank) = tree[j][k].incoming_potential;
			start += rank;
		}
	}

void getUserCheckPoints() {
	for (size_t j = nLevels; j >= level_LFR; j--) {//LFR
		#pragma omp parallel for
		for (size_t k = 0; k < tree[j].size(); k++) {
			if (tree[j][k].isLeaf) {
				tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
			}
			else {
				if (tree[j][k].incoming_ILActive) {
					tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].incoming_checkPoints.begin(), tree[j][k].incoming_checkPoints.end());
				}
				else {
					int J = j+1;
					int b = tree[j][k].boxNumber;
					int KboxNumber;
					int K[4];
					std::vector<int>::iterator indx;
					KboxNumber = 4*b+0;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[0] = indx-indexTree[J].begin();
					KboxNumber = 4*b+1;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[1] = indx-indexTree[J].begin();
					KboxNumber = 4*b+2;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[2] = indx-indexTree[J].begin();
					KboxNumber = 4*b+3;
					indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
					K[3] = indx-indexTree[J].begin();
					for (size_t ch = 0; ch < 4; ch++) {
						tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[J][K[ch]].user_checkPoints.begin(), tree[J][K[ch]].user_checkPoints.end());
					}
				}
			}
		}
	}
	for (size_t j = level_LFR-1; j >= 2; j--) {//HFR
		#pragma omp parallel for
		for (size_t k = 0; k < tree[j].size(); k++) {
			if (tree[j][k].isLeaf) {
				tree[j][k].user_checkPoints.insert(tree[j][k].user_checkPoints.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
			}
			else {//HFR
				for (size_t cone = 0; cone < nCones[j]; cone++) {
					if (tree[j][k].ConeTree[cone].incoming_ILActive) {
						tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[j][k].ConeTree[cone].incoming_checkPoints.begin(), tree[j][k].ConeTree[cone].incoming_checkPoints.end());
					}
					else {
						int J = j+1;
						int b = tree[j][k].boxNumber;
						int KboxNumber;
						int K[4];
						std::vector<int>::iterator indx;
						KboxNumber = 4*b+0;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[0] = indx-indexTree[J].begin();
						KboxNumber = 4*b+1;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[1] = indx-indexTree[J].begin();
						KboxNumber = 4*b+2;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[2] = indx-indexTree[J].begin();
						KboxNumber = 4*b+3;
						indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
						K[3] = indx-indexTree[J].begin();
						int cone_child =  cone/2;
						for (size_t ch = 0; ch < 4; ch++) {
							if (tree[J][K[ch]].isLeaf) {
								tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].user_checkPoints.begin(), tree[J][K[ch]].user_checkPoints.end());
							}
							else {
								if (J == level_LFR) {
									tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].user_checkPoints.begin(), tree[J][K[ch]].user_checkPoints.end());
								}
								else {
									tree[j][k].ConeTree[cone].user_checkPoints.insert(tree[j][k].ConeTree[cone].user_checkPoints.end(), tree[J][K[ch]].ConeTree[cone_child].user_checkPoints.begin(), tree[J][K[ch]].ConeTree[cone_child].user_checkPoints.end());
								}
							}
						}
					}
				}
			}
		}
	}
}

void getNodes_HFR() {
	for (int j=level_LFR-1; j>=2; --j) {
		for (int k=0; k<tree[j].size(); ++k) {
			if (!tree[j][k].isLeaf) {
				for (int cone=0; cone<nCones[j]; ++cone) {
					if (tree[j][k].ConeTree[cone].incoming_chargePoints.size() > 0)
						tree[j][k].ConeTree[cone].incoming_chargePoints.clear();
					if (tree[j][k].ConeTree[cone].incoming_checkPoints.size() > 0)
						tree[j][k].ConeTree[cone].incoming_checkPoints.clear();
					if (tree[j][k].ConeTree[cone].outgoing_chargePoints.size() > 0)
						tree[j][k].ConeTree[cone].outgoing_chargePoints.clear();
					if (tree[j][k].ConeTree[cone].outgoing_checkPoints.size() > 0)
						tree[j][k].ConeTree[cone].outgoing_checkPoints.clear();
					if (tree[j][k].ConeTree[cone].user_checkPoints.size() > 0)
						tree[j][k].ConeTree[cone].user_checkPoints.clear();
				}
			}
		}
	}
	for (int j=level_LFR-1; j>=2; --j) {
		getNodes_HFR_outgoing_level(j);
		getNodes_HFR_incoming_level(j);
	}
}


void getParticlesFromChildrenHFR_outgoing_row(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getParticlesFromChildrenHFR_outgoing_col(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getNodes_HFR_outgoing_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
	if (tree[j][k].active == true) {
		for (int cone=0; cone<nCones[j]; ++cone) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenHFR_outgoing_col(j, k, cone, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;
			for (int b = 0; b < tree[j][k].ConeTree[cone].InteractionList.size(); b++) {
					int jB = tree[j][k].ConeTree[cone].InteractionList[b].x;
					int boxB = tree[j][k].ConeTree[cone].InteractionList[b].y;
					double arg = atan2(tree[j][k].center.y-tree[jB][boxB].center.y, tree[j][k].center.x-tree[jB][boxB].center.x);
					double argB = fmod(arg+2*PI, 2*PI);
					int coneB = int(argB/ConeAperture[jB]);//=coneB; direction of nBox wrt k
					std::vector<int> chargeLocations;
					getParticlesFromChildrenHFR_outgoing_row(jB, boxB, coneB, chargeLocations);
					IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = IL_Nodes.size();
			n_cols = boxA_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = IL_Nodes;
			col_indices = boxA_Nodes;
			std::vector<int> row_bases, col_bases;
			// if(j==5 && k==112 && cone==7) {
			// 	cout << "j: " << j << "	k: " << k << "	cone: " << cone << "	n_rows: " << n_rows << "	n_cols: " << n_cols << endl;
			// }
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].ConeTree[cone].outgoing_ILActive = true;
				Mat dummy;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy, tree[j][k].ConeTree[cone].outgoing_Ar);
				// if(j==5 && k==112 && cone==7) {
				// 	cout << "row_bases.size(): " << row_bases.size() << "	col_bases.size(): " << col_bases.size() << "	ComputedRank: " << ComputedRank << endl;
				// }
				if(ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].ConeTree[cone].outgoing_checkPoints.push_back(IL_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].ConeTree[cone].outgoing_chargePoints.push_back(boxA_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(IL_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(boxA_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].ConeTree[cone].outgoing_Atilde_dec = Atilde.colPivHouseholderQr();
				}
			}
			// if (n_rows == 0) {
			else {
				tree[j][k].ConeTree[cone].outgoing_ILActive = false;
				getParticlesFromChildrenHFR_outgoing_col(j, k, cone, tree[j][k].ConeTree[cone].outgoing_chargePoints);
			}
			// if(j==5 && k==112 && cone==7) {
			// 	cout << "j: " << j << "	k: " << k << "	cone: " << cone << "	OIL: " << tree[j][k].ConeTree[cone].outgoing_ILActive << endl;
			// }
		}
	}
}

void getNodes_HFR_outgoing_level(int j) { //HFR; cone interactions
	//for (int j=level_LFR-1; j>=2; --j) {
		int rankPerLevel = 0;
		int n_rows_checkpoint;
		int n_cols_checkpoint;
		int kMax;
		int n_rows, n_cols, ComputedRank;
		std::vector<int> boxA_Particles_checkpoint;
		std::vector<int> IL_Particles_checkpoint;
		for (int k=0; k<tree[j].size(); ++k) {
			if (tree[j][k].isLeaf) {
				getNodes_LFR_outgoing_box(j, k, n_rows, n_cols, ComputedRank);
			}
			else {
				getNodes_HFR_outgoing_box(j, k, n_rows, n_cols, ComputedRank);
			}
			if (rankPerLevel < ComputedRank) {
				rankPerLevel = ComputedRank;
				n_rows_checkpoint = n_rows;
				n_cols_checkpoint = n_cols;
				kMax = k;
			}
		}
	//}
		cout << "O;	j: " << j << "	Nboxes: " << tree[j].size() << "	k: " << kMax << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
	}

void getParticlesFromChildrenHFR_incoming_row(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].incoming_checkPoints.begin(), tree[J][K].incoming_checkPoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.begin(), tree[J][K].ConeTree[cone_child].incoming_checkPoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getParticlesFromChildrenHFR_incoming_col(int j, int k, int cone_parent, std::vector<int>& searchNodes) {
	if (tree[j][k].isLeaf) {
		searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
	}
	else {
		int J = j+1;//child
		int b = tree[j][k].boxNumber;
		if (j == level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			for (int child = 0; child < 4; child++) {
				int KboxNumber = 4*b+child;
				std::vector<int>::iterator indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
				int K = indx-indexTree[J].begin();
				if (tree[J][K].isLeaf) {
					searchNodes.insert(searchNodes.end(), tree[J][K].outgoing_chargePoints.begin(), tree[J][K].outgoing_chargePoints.end());//outgoing_chargePoints
				}
				else {
					searchNodes.insert(searchNodes.end(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.begin(), tree[J][K].ConeTree[cone_child].outgoing_chargePoints.end());//outgoing_chargePoints
				}
			}
		}
	}
}

void getNodes_HFR_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
	if (tree[j][k].active == true) {
		for (int cone=0; cone<nCones[j]; ++cone) {
			std::vector<int> boxA_Nodes;
			getParticlesFromChildrenHFR_incoming_row(j, k, cone, boxA_Nodes);

			//sort( boxA_Nodes.begin(), boxA_Nodes.end() );
			//boxA_Nodes.erase( unique( boxA_Nodes.begin(), boxA_Nodes.end() ), boxA_Nodes.end() );

			std::vector<int> IL_Nodes;
			for (int b = 0; b < tree[j][k].ConeTree[cone].InteractionList.size(); b++) {
					int jB = tree[j][k].ConeTree[cone].InteractionList[b].x;
					int boxB = tree[j][k].ConeTree[cone].InteractionList[b].y;
					double arg = atan2(tree[j][k].center.y-tree[jB][boxB].center.y, tree[j][k].center.x-tree[jB][boxB].center.x);
					double argB = fmod(arg+2*PI, 2*PI);
					int coneB = int(argB/ConeAperture[jB]);//=coneB; direction of nBox wrt k
					std::vector<int> chargeLocations;
					getParticlesFromChildrenHFR_incoming_col(jB, boxB, coneB, chargeLocations);
					IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}

			//sort( IL_Nodes.begin(), IL_Nodes.end() );
			//IL_Nodes.erase( unique( IL_Nodes.begin(), IL_Nodes.end() ), IL_Nodes.end() );

			n_rows = boxA_Nodes.size();
			n_cols = IL_Nodes.size();
			int tol_pow = TOL_POW;
			double tol_ACA = pow(10,-1.0*tol_pow);
			row_indices = boxA_Nodes;
			col_indices = IL_Nodes;
			std::vector<int> row_bases, col_bases;
			if (n_rows != 0 && n_cols != 0) {
				tree[j][k].ConeTree[cone].incoming_ILActive = true;
				Mat dummy1, dummy2;
				ACA_only_nodes(row_bases, col_bases, ComputedRank, tol_ACA, dummy1, dummy2);
				if(ComputedRank > 0) {
					for (int r = 0; r < row_bases.size(); r++) {
						tree[j][k].ConeTree[cone].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
					}
					for (int c = 0; c < col_bases.size(); c++) {
						tree[j][k].ConeTree[cone].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
					}
					std::vector<int> row_indices_local;
					for (size_t r = 0; r < row_bases.size(); r++) {
						row_indices_local.push_back(boxA_Nodes[row_bases[r]]);
					}
					std::vector<int> col_indices_local;
					for (size_t c = 0; c < col_bases.size(); c++) {
						col_indices_local.push_back(IL_Nodes[col_bases[c]]);
					}
					Mat Atilde = getMatrix(row_indices_local, col_indices_local);
					tree[j][k].ConeTree[cone].incoming_Atilde_dec = Atilde.colPivHouseholderQr();
				}
			}
			// if (n_cols == 0) {
			else {
				tree[j][k].ConeTree[cone].incoming_ILActive = false;
				getParticlesFromChildrenHFR_incoming_row(j, k, cone, tree[j][k].ConeTree[cone].incoming_checkPoints);
			}
		}
	}
}

void getNodes_HFR_incoming_level(int j) { //HFR; cone interactions
		int rankPerLevel = 0;
		int n_rows_checkpoint;
		int n_cols_checkpoint;
		int kMax;
		std::vector<int> boxA_Particles_checkpoint;
		std::vector<int> IL_Particles_checkpoint;
		int n_rows, n_cols, ComputedRank;
		for (int k=0; k<tree[j].size(); ++k) {
			if (tree[j][k].isLeaf){
				getNodes_LFR_incoming_box(j, k, n_rows, n_cols, ComputedRank);
			}
			else {
				getNodes_HFR_incoming_box(j, k, n_rows, n_cols, ComputedRank);
			}
			if (rankPerLevel < ComputedRank) {
				rankPerLevel 			= ComputedRank;
				n_rows_checkpoint = n_rows;
				n_cols_checkpoint = n_cols;
				kMax 							= k;
			}
		}
	cout << "I;	j: " << j << "	Nboxes: " << tree[j].size() << "	k: " << kMax << "	rows,cols: " << n_rows_checkpoint << "," << n_cols_checkpoint << "	Crank: " << rankPerLevel << endl;
}

	void Assemble_HFR_M2M_ILActiveTrue(int j, int k, int cone_parent) {
		/*
		A = K(x^{B,o}, y^{B,o})
		x^{B,o}=tree[j][k].outgoing_checkPoints
		y^{B,o}=tree[j][k].outgoing_chargePoints
		*/
		std::vector<int> source_points;// = tree[j][k].chargeLocations//source points
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		if (j==level_LFR-1) {
			for (int child = 0; child < 4; child++) {
				source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			int Veclength = 0;
			for (int child = 0; child < 4; child++) {
				if (tree[J][K[child]].isLeaf) {
					source_points.insert(source_points.end(), tree[J][K[child]].outgoing_chargePoints.begin(), tree[J][K[child]].outgoing_chargePoints.end());//outgoing_chargePoints
				}
				else {
					source_points.insert(source_points.end(), tree[J][K[child]].ConeTree[cone_child].outgoing_chargePoints.begin(), tree[J][K[child]].ConeTree[cone_child].outgoing_chargePoints.end());//outgoing_chargePoints
				}
			}
		}
		int n_rows = tree[j][k].ConeTree[cone_parent].outgoing_checkPoints.size();//outgoing_checkPoints
		int n_cols = source_points.size();
		if (n_rows != 0 && n_cols != 0) {
			tree[j][k].ConeTree[cone_parent].M2M = getMatrix(tree[j][k].ConeTree[cone_parent].outgoing_checkPoints, source_points);
			// if (tree[j][k].ConeTree[cone_parent].M2M.norm() > pow(10, 10.0)) {
			// 	cout << "j: " << j << "	k: " << k << "	M2M.n(): " << tree[j][k].ConeTree[cone_parent].M2M.norm() << endl;
			// }
		}
	}

	void HFR_M2M_ILActiveTrue(int j, int k, int cone_parent) {
		/*
		A = K(x^{B,o}, y^{B,o})
		x^{B,o}=tree[j][k].outgoing_checkPoints
		y^{B,o}=tree[j][k].outgoing_chargePoints
		*/
		std::vector<int> source_points;// = tree[j][k].chargeLocations//source points
		Vec source_densities;
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		if (j==level_LFR-1) {
			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					source_densities(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			int Veclength = 0;
			// #pragma omp parallel for
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					Veclength += tree[J][K[child]].outgoing_charges.size();
					// if(j==4 && k==28 && cone_parent==14 || j==5 && k==112 && cone_parent==7 || j==5 && k==114 && cone_parent==7) {
					// 	cout << "j: " << j << "	k: " << k << "	child: " << child << "J: " << J << "	K[child]: " << K[child] << "	OILA_c: " << tree[J][K[child]].outgoing_ILActive << "	ocharges: " << tree[J][K[child]].outgoing_charges.size() << "	ochargeP: " << tree[J][K[child]].outgoing_chargePoints.size() << endl;
					// }
				}
				else {
					Veclength += tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size();
					// if(j==4 && k==28 && cone_parent==14 || j==5 && k==112 && cone_parent==7 || j==5 && k==114 && cone_parent==7) {
					// 	cout << "j: " << j << "	k: " << k << "	child: " << child << "J: " << J << "	K[child]: " << K[child] << "	cone_child: " << cone_child << "	OILA_c: " << tree[J][K[child]].ConeTree[cone_child].outgoing_ILActive << "	ocharges: " << tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size() << "	ochargeP: " << tree[J][K[child]].ConeTree[cone_child].outgoing_chargePoints.size() << endl;
					// }
				}
			}
			source_densities = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
						source_densities(g) = tree[J][K[child]].outgoing_charges(l);
						++g;
					}
				}
				else {
					for (int l = 0; l < tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size(); l++) {
						source_densities(g) = tree[J][K[child]].ConeTree[cone_child].outgoing_charges(l);
						++g;
					}
				}
			}
		}
		int n_rows = tree[j][k].ConeTree[cone_parent].outgoing_checkPoints.size();//outgoing_checkPoints
		int n_cols = source_densities.size();
		// if(j==4 && k==28 && cone_parent==14) {
		// 	cout << "j: " << j << "	k: " << k << "	M2M: " << tree[j][k].ConeTree[cone_parent].M2M.rows() << ", " << tree[j][k].ConeTree[cone_parent].M2M.cols() << "	source_densities: " << source_densities.size() << endl;
		// 	cout << "ochargep: " << tree[j][k].ConeTree[cone_parent].outgoing_chargePoints.size() << "	ocheckp: " << tree[j][k].ConeTree[cone_parent].outgoing_checkPoints.size() << endl;
		// }
		// if(j==5 && k==112 && cone_parent==7 || j==5 && k==114 && cone_parent==7) {
		// 	cout << "j: " << j << "	k: " << k << "	coneP: " << cone_parent << endl;
		// 	cout << "M2M: " << tree[j][k].ConeTree[cone_parent].M2M.rows() << ", " << tree[j][k].ConeTree[cone_parent].M2M.cols() << "	source_densities: " << source_densities.size() << endl;
		// 	cout << "ochargep: " << tree[j][k].ConeTree[cone_parent].outgoing_chargePoints.size() << "	ocheckp: " << tree[j][k].ConeTree[cone_parent].outgoing_checkPoints.size() << endl;
		// }
		if (n_rows != 0 && n_cols != 0) {
			tree[j][k].ConeTree[cone_parent].outgoing_potential = tree[j][k].ConeTree[cone_parent].M2M*source_densities;//u^{B,o}
			tree[j][k].ConeTree[cone_parent].outgoing_charges = tree[j][k].ConeTree[cone_parent].outgoing_Atilde_dec.solve(tree[j][k].ConeTree[cone_parent].outgoing_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
		}
	}

	void HFR_M2M_ILActiveFalse(int j, int k, int cone_parent) {
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		if (j==level_LFR-1) {
			int Veclength = tree[J][K[0]].outgoing_charges.size()+tree[J][K[1]].outgoing_charges.size()+tree[J][K[2]].outgoing_charges.size()+tree[J][K[3]].outgoing_charges.size();
			tree[j][k].ConeTree[cone_parent].outgoing_charges = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
					tree[j][k].ConeTree[cone_parent].outgoing_charges(g) = tree[J][K[child]].outgoing_charges(l);
					++g;
				}
			}
		}
		else {
			int cone_child =  cone_parent/2; //of j+1 level;//l'
			int Veclength = 0;
			// #pragma omp parallel for
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					Veclength += tree[J][K[child]].outgoing_charges.size();
				}
				else {
					Veclength += tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size();
				}
			}
			tree[j][k].ConeTree[cone_parent].outgoing_charges = Vec::Zero(Veclength);// = tree[j][k].multipoles//source densities
			int g = 0;
			for (int child = 0; child < 4; child++) {//child
				if (tree[J][K[child]].isLeaf) {
					for (int l = 0; l < tree[J][K[child]].outgoing_charges.size(); l++) {
						tree[j][k].ConeTree[cone_parent].outgoing_charges(g) = tree[J][K[child]].outgoing_charges(l);
						++g;
					}
				}
				else {
					for (int l = 0; l < tree[J][K[child]].ConeTree[cone_child].outgoing_charges.size(); l++) {
						tree[j][k].ConeTree[cone_parent].outgoing_charges(g) = tree[J][K[child]].ConeTree[cone_child].outgoing_charges(l);
						++g;
					}
				}
			}
		}
	}

	void HFR_M2M() {//outgoing operations
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=level_LFR-1; j>=2; --j) {
			// #pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true) {
					// #pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						// cout << "j: " << j << "	k: " << k << "	cone_parent: " << cone_parent << "	ILA: " << tree[j][k].ConeTree[cone_parent].ILActive << " ochp: " << tree[j][k].ConeTree[cone_parent].outgoing_chargePoints.size() << endl;
						// cout << "5,112, OIL: " << tree[j][k].ConeTree[cone_parent].outgoing_ILActive << "	IIL: " << tree[j][k].ConeTree[cone_parent].incoming_ILActive << endl;
						if (tree[j][k].ConeTree[cone_parent].outgoing_ILActive) {
							if (tree[j][k].ConeTree[cone_parent].outgoing_chargePoints.size() > 0) {
								HFR_M2M_ILActiveTrue(j, k, cone_parent);
							}
						}
						else {
							HFR_M2M_ILActiveFalse(j, k, cone_parent);
						}
					}
				}
			}
		}
	}

	void Assemble_HFR_M2M() {//outgoing operations
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=level_LFR-1; j>=2; --j) {
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true) {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						if (tree[j][k].ConeTree[cone_parent].outgoing_ILActive) {
							Assemble_HFR_M2M_ILActiveTrue(j, k, cone_parent);
						}
					}
				}
			}
		}
	}

	void Assemble_HFR_M2L() {
		/*
		l:tree[j][boxA].ConeTree[coneA]; l':tree[j][boxB].ConeTree[coneB]
		u^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k < tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active) {
					#pragma omp parallel for
					for (int coneB = 0; coneB < nCones[j]; coneB++) {
						//tree[j][k].ConeTree[coneB].incoming_potential = Vec::Zero(tree[j][k].ConeTree[coneB].user_checkPoints.size());
						if (tree[j][k].ConeTree[coneB].incoming_ILActive) {
							int numIL = tree[j][k].ConeTree[coneB].InteractionList.size();
							for (int l = 0; l < numIL; l++) {
								int jIL = tree[j][k].ConeTree[coneB].InteractionList[l].x;
								int kIL = tree[j][k].ConeTree[coneB].InteractionList[l].y;
								if (tree[jIL][kIL].isLeaf) {//LFR
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].outgoing_charges.size();
									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].ConeTree[coneB].M2L.push_back(getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints));
										//tree[j][k].ConeTree[coneB].incoming_potential += R*tree[jIL][kIL].outgoing_charges;//u^{B,o}
									// }
								}
								else {
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									// if (n_rows != 0 && n_cols != 0) {
										tree[j][k].ConeTree[coneB].M2L.push_back(getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints));
										//tree[j][k].ConeTree[coneB].incoming_potential += R*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									// }
								}
							}
						}
					}
				}
			}
		}
	}

	void HFR_M2L() {
		/*
		l:tree[j][boxA].ConeTree[coneA]; l':tree[j][boxB].ConeTree[coneB]
		u^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_potential
		x^{B,i,l'}: tree[j][boxB].ConeTree[coneB].incoming_checkPoints
		f^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_charges
		y^{A,o,l}: tree[j][boxA].ConeTree[coneA].outgoing_chargePoints
		*/
		#pragma omp parallel for
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k < tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active) {
					#pragma omp parallel for
					for (int coneB = 0; coneB < nCones[j]; coneB++) {
						tree[j][k].ConeTree[coneB].incoming_potential = Vec::Zero(tree[j][k].ConeTree[coneB].user_checkPoints.size());
						if (tree[j][k].ConeTree[coneB].incoming_ILActive) {
							int numIL = tree[j][k].ConeTree[coneB].InteractionList.size();
							#pragma omp parallel for
							for (int l = 0; l < numIL; l++) {
								int jIL = tree[j][k].ConeTree[coneB].InteractionList[l].x;
								int kIL = tree[j][k].ConeTree[coneB].InteractionList[l].y;
								if (tree[jIL][kIL].isLeaf) {//LFR
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										/*Mat R = getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].outgoing_chargePoints);
										Mat Err = tree[j][k].ConeTree[coneB].M2L[l]-R;
										if (Err.norm() != 0.0) {
											cout << "j: " << j << "	k: " << k << "	coneB: " << coneB << "	Err: " << Err.norm() << endl;
										}*/
										tree[j][k].ConeTree[coneB].incoming_potential += tree[j][k].ConeTree[coneB].M2L[l]*tree[jIL][kIL].outgoing_charges;//u^{B,o}
									}
								}
								else {
									double arg = atan2(tree[j][k].center.y-tree[jIL][kIL].center.y, tree[j][k].center.x-tree[jIL][kIL].center.x);
									double argA = fmod(arg+2*PI, 2*PI);
									int coneA = int(argA/ConeAperture[jIL]);//=coneB; direction of nBox wrt k
									int n_rows = tree[j][k].ConeTree[coneB].user_checkPoints.size();
									int n_cols = tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints.size();//outgoing_chargePoints
									int RHS_size = tree[jIL][kIL].ConeTree[coneA].outgoing_charges.size();
									if (n_rows != 0 && n_cols != 0 && RHS_size != 0) {
										/*Mat R = getMatrix(tree[j][k].ConeTree[coneB].user_checkPoints, tree[jIL][kIL].ConeTree[coneA].outgoing_chargePoints);
										Mat Err = tree[j][k].ConeTree[coneB].M2L[l]-R;
										if (Err.norm() != 0.0) {
											cout << "j: " << j << "	k: " << k << "	coneB: " << coneB << "	Err: " << Err.norm() << endl;
										}*/
										tree[j][k].ConeTree[coneB].incoming_potential += tree[j][k].ConeTree[coneB].M2L[l]*tree[jIL][kIL].ConeTree[coneA].outgoing_charges;//u^{B,o}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	void Assemble_HFR_L2L_ILActiveTrue (int j, int k, int cone_parent) {
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			//cout << "J: " << J << "	K: " << K[child] << "	iL: " << tree[J][K[child]].isLeaf << endl;
			if (j==level_LFR-1) {
				if (level_LFR != nLevels) {
					if (!tree[J][K[child]].isLeaf) {
						int n_rows = tree[J][K[child]].user_checkPoints.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
					else {
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
				}
				else {
					int n_rows = tree[J][K[child]].chargeLocations.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
			}
			else {
				int cone_child = cone_parent/2;
				if (!tree[J][K[child]].isLeaf) {
					//cout << "J: " << J << "	K: " << K[child] << "	rc: " << tree[J][K[child]].ConeTree[cone_child].L2L.rows() << ", " << tree[J][K[child]].ConeTree[cone_child].L2L.cols() << "	rhs: " << tree[j][k].ConeTree[cone_parent].incoming_charges.size() << endl;

					int n_rows = tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].ConeTree[cone_child].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						//tree[J][K[child]].ConeTree[cone_child].incoming_potential += tree[J][K[child]].ConeTree[cone_child].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
					//cout << "done" << endl;
				}
				else {
					int n_rows = tree[J][K[child]].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					//cout << "J: " << J << "	K: " << K[child] << "	rc: " << tree[J][K[child]].L2L.rows() << ", " << tree[J][K[child]].L2L.cols() << "	rhs: " << tree[j][k].ConeTree[cone_parent].incoming_charges.size() << ", " << tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size() << endl;
					if (n_rows != 0 && n_cols != 0) {
						tree[j][k].ConeTree[cone_parent].L2L[child] = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						//cout << "R.rc: " << tree[J][K[child]].L2L.rows() << ", " << tree[J][K[child]].L2L.cols() << endl;
						//tree[J][K[child]].incoming_potential += tree[J][K[child]].L2L*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
					//cout << "done" << endl;
				}
			}
		}
	}

	void HFR_L2L_ILActiveTrue (int j, int k, int cone_parent) {
		/*
		x^{C,i}: tree[J][4*k+c].incoming_checkPoints
		f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
		y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
		x^{C,i}: tree[J][4*k+c].incoming_checkPoints
		f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
		y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
		x^{C,i,l'}: tree[J][4*k+c].ConeTree[cone_child].incoming_checkPoints
		f^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_charges
		y^{B,i,l}: tree[j][k].ConeTree[cone_parent].incoming_chargePoints
		*/
		tree[j][k].ConeTree[cone_parent].incoming_charges = tree[j][k].ConeTree[cone_parent].incoming_Atilde_dec.solve(tree[j][k].ConeTree[cone_parent].incoming_potential);//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		#pragma omp parallel for
		for (int child = 0; child < 4; child++) {
			if (j==level_LFR-1) {
				if (level_LFR != nLevels) {
					if (!tree[J][K[child]].isLeaf) {
						int n_rows = tree[J][K[child]].user_checkPoints.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							/*Mat R = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
							if (Err.norm() != 0.0) {
								cout << "J: " << J << "	K[child]: " << K[child] << "	cone_parent: " << cone_parent << "	Err: " << Err.norm() << endl;
							}*/
							tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
					else {
						int n_rows = tree[J][K[child]].chargeLocations.size();
						int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
						if (n_rows != 0 && n_cols != 0) {
							/*Mat R = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
							Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
							if (Err.norm() != 0.0) {
								cout << "J: " << J << "	K[child]: " << K[child] << "	cone_parent: " << cone_parent << "	Err: " << Err.norm() << endl;
							}*/
							tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
						}
					}
				}
				else {
					int n_rows = tree[J][K[child]].chargeLocations.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						/*Mat R = getMatrix(tree[J][K[child]].chargeLocations, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
						if (Err.norm() != 0.0) {
							cout << "J: " << J << "	K[child]: " << K[child] << "	cone_parent: " << cone_parent << "	Err: " << Err.norm() << endl;
						}*/
						tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
			}
			else {
				int cone_child = cone_parent/2;
				if (!tree[J][K[child]].isLeaf) {
					int n_rows = tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						/*Mat R = getMatrix(tree[J][K[child]].ConeTree[cone_child].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
						if (Err.norm() != 0.0) {
							cout << "J: " << J << "	K[child]: " << K[child] << "	cone_child: " << cone_child << "	Err: " << Err.norm() << endl;
						}*/
						tree[J][K[child]].ConeTree[cone_child].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
				else {
					int n_rows = tree[J][K[child]].user_checkPoints.size();
					int n_cols = tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size();
					if (n_rows != 0 && n_cols != 0) {
						/*Mat R = getMatrix(tree[J][K[child]].user_checkPoints, tree[j][k].ConeTree[cone_parent].incoming_chargePoints);
						Mat Err = tree[j][k].ConeTree[cone_parent].L2L[child]-R;
						if (Err.norm() != 0.0) {
							cout << "J: " << J << "	K[child]: " << K[child] << "	cone_child: " << cone_child << "	Err: " << Err.norm() << endl;
						}*/
						tree[J][K[child]].incoming_potential += tree[j][k].ConeTree[cone_parent].L2L[child]*tree[j][k].ConeTree[cone_parent].incoming_charges;//u^{B,o}
					}
				}
			}
		}
	}

	void HFR_L2L_ILActiveFalse(int j, int k, int cone_parent) {
		int J = j+1;
		int b = tree[j][k].boxNumber;
		int KboxNumber;
		int K[4];
		std::vector<int>::iterator indx;
		KboxNumber = 4*b+0;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[0] = indx-indexTree[J].begin();
		KboxNumber = 4*b+1;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[1] = indx-indexTree[J].begin();
		KboxNumber = 4*b+2;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[2] = indx-indexTree[J].begin();
		KboxNumber = 4*b+3;
		indx = std::find(indexTree[J].begin(), indexTree[J].end(), KboxNumber);
		K[3] = indx-indexTree[J].begin();
		int offset = 0;
		for (int child = 0; child < 4; child++) {
			if (j==level_LFR-1) {
				if (level_LFR != nLevels) {
					if (!tree[J][K[child]].isLeaf) {
						for (size_t i = 0; i < tree[J][K[child]].user_checkPoints.size(); i++) {
							tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
						}
						offset += tree[J][K[child]].user_checkPoints.size();
					}
					else {
						for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
							tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
						}
						offset += tree[J][K[child]].chargeLocations.size();
					}
				}
				else {
					for (size_t i = 0; i < tree[J][K[child]].chargeLocations.size(); i++) {
						tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].chargeLocations.size();
				}
			}
			else {
				int cone_child = cone_parent/2;
				if (tree[J][K[child]].isLeaf) {
					for (size_t i = 0; i < tree[J][K[child]].user_checkPoints.size(); i++) {
						tree[J][K[child]].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].user_checkPoints.size();
				}
				else {
					for (size_t i = 0; i < tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size(); i++) {
						tree[J][K[child]].ConeTree[cone_child].incoming_potential(i) += tree[j][k].ConeTree[cone_parent].incoming_potential(offset+i);
					}
					offset += tree[J][K[child]].ConeTree[cone_child].user_checkPoints.size();
				}
			}
		}
	}

	void Assemble_HFR_L2L() {//outgoing operations; finds locals untill level: level_LFR
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true && !tree[j][k].isLeaf) {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						if (tree[j][k].ConeTree[cone_parent].incoming_ILActive) {
							Assemble_HFR_L2L_ILActiveTrue(j,k,cone_parent);
						}
					}
				}
			}
		}
	}

	void HFR_L2L() {//outgoing operations; finds locals untill level: level_LFR
		/*
		tree[j][k].multipoles//source densities
		tree[j][k].chargeLocations//source points
		x^{B,o}=tree[j][k].outgoing_checkPoints
		kernel evaluation between x^{B,o} and source points
		*/
		for (int j=2; j<level_LFR; ++j) {//parent
			#pragma omp parallel for
			for (int k=0; k<tree[j].size(); ++k) {
				if (tree[j][k].isLeaf){
					continue;
				}
				if (tree[j][k].active == true && !tree[j][k].isLeaf) {
					#pragma omp parallel for
					for (int cone_parent = 0; cone_parent < nCones[j]; cone_parent++) {//pick l
						if (tree[j][k].ConeTree[cone_parent].incoming_ILActive) {
							if (tree[j][k].ConeTree[cone_parent].incoming_chargePoints.size() > 0) {
								HFR_L2L_ILActiveTrue(j,k,cone_parent);
							}
						}
						else {
							HFR_L2L_ILActiveFalse(j,k,cone_parent);
						}
					}
				}
			}
		}
	}
//////////////////////////////////////////////////////////////////
};

#endif

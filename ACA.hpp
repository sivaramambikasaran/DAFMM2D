//#ifndef _singularNodes_HPP__
//#define _singularNodes_HPP__
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <set>
#include <algorithm>
#include <bits/stdc++.h>

const std::complex<double> I(0.0, 1.0);
// #define EIGEN_DONT_PARALLELIZE
using namespace std;
using namespace Eigen;
const double PI	=	3.1415926535897932384;

typedef MatrixXcd Mat;
typedef VectorXcd Vec;
typedef std::complex<double> kernel_dtype;

struct pts2D {
	double x,y;
};

class FMM_Matrix {
public:
	std::vector<pts2D> particles_X;
	std::vector<pts2D> particles_Y;

	FMM_Matrix(std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y){
			this->particles_X = particles_X;
			this->particles_Y = particles_Y;
	}

	// virtual kernel_dtype getMatrixEntry(const unsigned i, const unsigned j) {
	// 	cout << "virtual getInteraction" << endl;
	// 	return 0.0;
	// }
	//
	kernel_dtype getMatrixEntry(const unsigned i, const unsigned j);

	virtual kernel_dtype getInteraction(const pts2D r1, const pts2D r2) {
		cout << "virtual getInteraction" << endl;
		return 0.0;
	}

	Vec getRow(const int j, std::vector<int> col_indices) {
		int n_cols = col_indices.size();
		Vec row(n_cols);
    //#pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, col_indices[k]);
    }
    return row;
  }

  Vec getCol(const int k, std::vector<int> row_indices) {
		int n_rows = row_indices.size();
    Vec col(n_rows);
    //#pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
			//cout << "j: " << j << "	row_indices[j]: " << row_indices[j] << "	k: " << k << endl;
      col(j) = this->getMatrixEntry(row_indices[j], k);
    }
    return col;
  }

  Mat getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
		int n_rows = row_indices.size();
		int n_cols = col_indices.size();
    Mat mat(n_rows, n_cols);
    //#pragma omp parallel for
    for (int j=0; j < n_rows; ++j) {
        //#pragma omp parallel for
        for (int k=0; k < n_cols; ++k) {
            mat(j,k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
        }
    }
    return mat;
  }

	Mat getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index) {
	  Mat mat(row_end_index-row_start_index, col_end_index-col_start_index);
	  #pragma omp parallel for
	  for (int j=row_start_index; j < row_end_index; ++j) {
	      #pragma omp parallel for
	      for (int k=col_start_index; k < col_end_index; ++k) {
	          mat(j,k) = this->getMatrixEntry(j, k);
	      }
	  }
	  return mat;
	}

};

class LowRank: public FMM_Matrix {
public:
	double tol_ACA;
	std::vector<int> row_indices;
	std::vector<int> col_indices;
	LowRank(int tol_pow, std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y, std::vector<int> row_indices, std::vector<int> col_indices):
	FMM_Matrix(particles_X, particles_Y) {
			this->tol_ACA = pow(10,-1.0*tol_pow);
			//this->particles_X = particles_X;
			//this->particles_Y = particles_Y;
			this->row_indices = row_indices;
			this->col_indices = col_indices;
	}

	void maxAbsVector(const Vec& v, const std::set<int>& allowed_indices,
																	double max, int& index
																 ) {
			std::set<int>::iterator it;
			index = *allowed_indices.begin();
			max   = abs(v(index));

			for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) {
					if(abs(v(*it)) > max) {
							index   =   *it;
							max     =   abs(v(index));
					}
			}
	}

	void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, double epsilon, Mat &Ac, Mat &Ar) {
		int col_index;
		int row_index;
		int N1 = row_indices.size();
		int N2 = col_indices.size();
		Vec row(N2), col(N1), v(N2), u(N1);
		Vec row_temp, col_temp;
		std::vector<Vec> Uvec;
		std::vector<Vec> Vvec;
		computed_rank = 0;
		double max;
		int min = N1;
		if (N1 > N2) {
			min = N2;
		}
		Ac = Mat(N1, min);
		Ar = Mat(min, N2);

		std::set<int> remaining_row_ind;
		std::set<int> remaining_col_ind;

		for(int l = 0; l < N1; l++) {
				remaining_row_ind.insert(l);
		}
		for(int l= 0; l < N2; l++) {
				remaining_col_ind.insert(l);
		}
			int l_local;
			for(l_local = 0; l_local < N1; l_local++) {
				row_index=l_local;
				row = getRow(row_indices[row_index], col_indices);
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				if(abs(row(int(col_index))) > pow(10,-16)) {
					break;
				}
			}
			if (l_local == N1) {
				Ac = Ac.block(0,0,N1,computed_rank);
				Ar = Ar.block(0,0,computed_rank,N2);
				return;
			}
			v=row/row(int(col_index));
			row_bases.push_back(row_index);
			col_bases.push_back(int(col_index));
			col = getCol(col_indices[col_index], row_indices);
			u	=	col;
			Uvec.push_back(u);
			Vvec.push_back(v);
			Ac.col(computed_rank) = col;
			Ar.row(computed_rank) = row;
			remaining_col_ind.erase(col_index);
			remaining_row_ind.erase(row_index);
			computed_rank = 1;

			double normS	=	0.0;
			double prev_normS = DBL_MAX;
			this->maxAbsVector(col, remaining_row_ind, max, row_index);

			while (abs(prev_normS-normS) >= epsilon*prev_normS && computed_rank < min) {
				row_bases.push_back(int(row_index));
				row = getRow(row_indices[row_index], col_indices);
				row_temp = row;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Vvec[l].size(); s++) {
						row(s)	-=	Uvec[l](row_index)*Vvec[l](s);
					}
				}
				this->maxAbsVector(row, remaining_col_ind, max, col_index);
				col_bases.push_back(int(col_index));
				for (size_t l = 0; l < row.size(); l++) {
					v(l) = row(l)/row(int(col_index));
				}
				col = getCol(col_indices[col_index], row_indices);
				col_temp = col;
				for (int l=0; l<=computed_rank-1; ++l){
					for (size_t s = 0; s < Uvec[l].size(); s++) {
						col(s)	-=	Vvec[l](col_index)*Uvec[l](s);
					}
				}
				u	=	col;
				if(u.norm()< 1e-16 || v.norm()< 1e-16) {
					row_bases.pop_back();
					col_bases.pop_back();
					break;
				}
				Uvec.push_back(u);
				Vvec.push_back(v);
				Ac.col(computed_rank) = col_temp;
				Ar.row(computed_rank) = row_temp;

				++computed_rank;
				remaining_col_ind.erase(col_index);
				remaining_row_ind.erase(row_index);
				if (computed_rank != 2)
					prev_normS = normS;
				normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
				for (int l=0; l<=computed_rank-1; ++l){
					normS	+=	2*abs(Uvec[l].dot(u) * Vvec[l].dot(v));
				}
				normS	=	sqrt(normS);
				this->maxAbsVector(col, remaining_row_ind, max, row_index);
			}

		Ac = Ac.block(0,0,N1,computed_rank);
		Ar = Ar.block(0,0,computed_rank,N2);
	}
};

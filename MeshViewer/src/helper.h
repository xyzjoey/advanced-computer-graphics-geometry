#ifndef HELPER_H
#define HELPER_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

struct Triplets //for sparse matrix initialization
{
	static std::vector<Eigen::Triplet<double>> triplets;
	static void reset(int expectedEntries) { triplets.clear(); triplets.reserve(expectedEntries); }
	static void push(int i, int j, double v_ij) { triplets.push_back(Eigen::Triplet<double>(i, j, v_ij)); }
	static Eigen::SparseMatrix<float> createMatrix(int row, int col)
	{
		Eigen::SparseMatrix<float> A(row, col);
		A.setFromTriplets(triplets.begin(), triplets.end());
		return A;
	}
};

Eigen::SparseMatrix<float> sparseI(int size);

double cotWeight(Vertex* p, Vertex* pNbr, Vertex* pNbrPrev, Vertex* pNbrNext);

float euclidean(float x, float y, float z);

#endif
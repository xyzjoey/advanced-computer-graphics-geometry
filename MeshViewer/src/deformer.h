#ifndef DEFORMER_H
#define DEFORMER_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mesh.h"

// Deform mesh using Laplacian coordinates
class Deformer {
public:
	Deformer();
	~Deformer();

	void setMesh(Mesh* mesh);

	/*====== Programming Assignment 2 ======*/
	// This is the place where the editing techniques take place
	void deform();
	void setHandler(int handlerInd); //
	/*====== Programming Assignment 2 ======*/

private:
	/*====== Programming Assignment 2 ======*/
	// Build left hand side matrix and pre-factorize it
	void buildSystemMat();
	

	int handlerInd; // vertex index of handler // -1 --> none
	//Mesh* smoothMesh;

	Eigen::SparseMatrix<float> mLaplacian;
	Eigen::SparseMatrix<float> mAT;
	Eigen::SparseMatrix<float> mANormalizeT;
	Eigen::VectorXf mb;
	Eigen::VectorXf mDelta;

	//std::vector<double> mSigmaLengths;
	std::vector<float> mRoiDeltaLengths;

	//const bool mCOTWEIGHT = true;
	const bool mLOCALTRANSFORM = false;
	
	std::vector<int> mVtxToHandleMap;

	/*====== Programming Assignment 2 ======*/

	void clear();

	Mesh* mMesh;
	std::vector< Vertex* > mRoiList; // unconstrained vertices
	std::vector< Vertex* > mHandleList; // constrained vertices
	// Solver for pre-factorizing the system matrix of the deformation
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >* mCholeskySolver;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >* mCholeskySolverNormalize;
};

#endif // DEFORMER_H

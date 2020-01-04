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
	void deform(bool rotationInvariant = true);
	void setHandleOffset(int handleFlag, Eigen::Vector3f offset);
	/*====== Programming Assignment 2 ======*/

private:
	/*====== Programming Assignment 2 ======*/
	// Build left hand side matrix and pre-factorize it
	void buildSystemMat();

	Eigen::MatrixXf mRotations; // alhpa, beta, gamma

	Eigen::SparseMatrix<float> mAT;
	Eigen::SparseMatrix<float> mASmoothT;
	Eigen::VectorXf mb;
	Eigen::VectorXf mbSmooth;
	
	/*====== Programming Assignment 2 ======*/

	void clear();

	Mesh* mMesh;
	Mesh* mMeshSmooth;
	std::vector< Vertex* > mRoiList; // unconstrained vertices
	std::vector< Vertex* > mHandleList; // constrained vertices
	// Solver for pre-factorizing the system matrix of the deformation
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >* mCholeskySolver;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >* mCholeskySolverSmooth;
};

#endif // DEFORMER_H

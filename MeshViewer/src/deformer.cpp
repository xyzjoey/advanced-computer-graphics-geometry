#include "deformer.h"
#include "helper.h"
#include <iostream>

Deformer::Deformer() : mMesh(nullptr),
                       mCholeskySolver(nullptr) {
}

Deformer::~Deformer() {
	clear();
}

void Deformer::clear() {
	if (mCholeskySolver) {
		delete mCholeskySolver;
	}
	mCholeskySolver = nullptr;
	mRoiList.clear();
	mHandleList.clear();
}

void Deformer::setMesh(Mesh* mesh) {
	mMesh = mesh;
	clear();
	// Record the handle vertices
	for (Vertex* vert : mMesh->vertices()) {
		if (vert->flag() > 0 || vert->isBoundary()) {
			mHandleList.push_back(vert);
		} else {
			mRoiList.push_back(vert);
		}
	}
	// Build system matrix for deformation
	buildSystemMat();
}

void Deformer::setHandler(int handlerInd) {
	this->handlerInd = handlerInd;
	std::cout << handlerInd << std::endl;
}

void Deformer::buildSystemMat() {
	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Build the matrix of the linear system for 
	/* deformation and do factorization, in order
	/* to reuse and speed up in Deformer::deform().
	/* Handle vertices are maked by Vertex::flag() > 0
	/* Movements of the specified handle are already
	/* recorded in Vertex::position()
	/**********************************************/

	std::vector<Vertex*> vertices = mMesh->vertices();
	int vtxNum = vertices.size();
	int roiNum = mRoiList.size();
	int handleNum = mHandleList.size();

	std::vector<int> vtxToHandleMap(vtxNum); //vtxToRoiMap[i] == k: ith vtx is kth roi, if not roi --> -1
	std::vector<std::vector<int>> nbrInds(vtxNum); //neighbor indices
	Eigen::VectorXf pos(vtxNum * 3);
	Eigen::VectorXf delta(vtxNum * 3);
	std::vector<float> roiDeltaLengths(roiNum * 3);
	Eigen::SparseMatrix<float> laplacian, T;
	Eigen::SparseMatrix<float> A, ANormalize;
	Eigen::VectorXf b;
	Eigen::SparseMatrix<double> systemMat, systemMatNormalize;

	// set std::vector<int> vtxToRoiMap
	for (int i = 0; i < vtxNum; ++i) {
		vtxToHandleMap[i] = -1; // initialize with -1
	}
	for (int i = 0; i < handleNum; ++i) {
		int vtxInd = mHandleList[i]->index();
		vtxToHandleMap[vtxInd] = i;
	}

	// set std::vector<std::vector<int>> nbrInds
	for (int i = 0; i < vtxNum; ++i) {
		OneRingVertex ring(vertices[i]);
		Vertex* neighbor = nullptr;
		while (neighbor = ring.nextVertex()) nbrInds[i].push_back(neighbor->index());
	}

	// set Eigen::VectorXf pos
	for (int i = 0; i < vtxNum; ++i) {
		pos[i*3 + 0] = vertices[i]->position()[0]; //x
		pos[i*3 + 1] = vertices[i]->position()[1]; //y
		pos[i*3 + 2] = vertices[i]->position()[2]; //z
	}

	// set Eigen::SparseMatrix<float> laplacian (with constraints)
	Triplets::reset(vtxNum * 3 * 10);
	for (int i = 0; i < vtxNum; ++i) {
		// set one for self
		Triplets::push(i*3 + 0, i*3 + 0, 1); //x
		Triplets::push(i*3 + 1, i*3 + 1, 1); //y
		Triplets::push(i*3 + 2, i*3 + 2, 1); //z

		// uniform weights for neighbors
		if (vtxToHandleMap[i] == -1) { // if is unconstrained
			double weight = -1.0 / nbrInds[i].size();
			for (int j = 0; j < nbrInds[i].size(); ++j) {
				int nbrInd = nbrInds[i][j];
				Triplets::push(i * 3 + 0, nbrInd * 3 + 0, weight); //x
				Triplets::push(i * 3 + 1, nbrInd * 3 + 1, weight); //y
				Triplets::push(i * 3 + 2, nbrInd * 3 + 2, weight); //z
			}
		}
	}
	laplacian = Triplets::createMatrix(vtxNum * 3, vtxNum * 3);

	// set Eigen::VectorXf delta
	delta = laplacian * pos;

	// delta length for normalization
	for (int i = 0; i < roiNum; ++i) {
		int vtxInd = mRoiList[i]->index();
		roiDeltaLengths[i] = euclidean(delta[vtxInd*3 + 0], delta[vtxInd*3 + 1], delta[vtxInd*3 + 2]);
	}

	// local transformation
	Triplets::reset(vtxNum * 9 * 10);
	for (int i = 0; i < roiNum; ++i) { //skip handle vertices
		int vtxInd = mRoiList[i]->index();
		int nbrNum = nbrInds[vtxInd].size();

		// set C
		Eigen::MatrixXf C((nbrNum + 1)*3, 7);
		for (int j = 0; j < nbrNum + 1; ++j) {
			int vtxIndCurr;
			if (j == 0) vtxIndCurr = vtxInd; // j==0 --> self
			else vtxIndCurr = nbrInds[vtxInd][j-1]; // else --> (j-1)th neighbor

			Vertex* vtxCurr = vertices[vtxIndCurr];

			const float x = vtxCurr->position()[0];
			const float y = vtxCurr->position()[1];
			const float z = vtxCurr->position()[2];

			C(j*3+0,0) = x; C(j*3+0,1) = 0; C(j*3+0,2) = z; C(j*3+0,3) = -y; C(j*3+0,4) = 1; C(j*3+0,5) = 0; C(j*3+0,6) = 0;
			C(j*3+1,0) = y; C(j*3+1,1) = -z; C(j*3+1,2) = 0; C(j*3+1,3) = x; C(j*3+1,4) = 0; C(j*3+1,5) = 1; C(j*3+1,6) = 0;
			C(j*3+2,0) = z; C(j*3+2,1) = y; C(j*3+2,2) = -x; C(j*3+2,3) = 0; C(j*3+2,4) = 0; C(j*3+2,5) = 0; C(j*3+2,6) = 1;
		}

		const float dx = delta[vtxInd*3 + 0];
		const float dy = delta[vtxInd*3 + 1];
		const float dz = delta[vtxInd*3 + 2];

		// set D
		Eigen::MatrixXf D(3, 7);
		D(0,0) = dx; D(0,1) = 0; D(0,2) = dz; D(0,3) = -dy; D(0,4) = 1; D(0,5) = 0; D(0,6) = 0;
		D(1,0) = dy; D(1,1) = -dz; D(1,2) = 0; D(1,3) = dx; D(1,4) = 0; D(1,5) = 1; D(1,6) = 0;
		D(2,0) = dz; D(2,1) = dy; D(2,2) = -dx; D(2,3) = 0; D(2,4) = 0; D(2,5) = 0; D(2,6) = 1;

		// compute M = D * ((C.T * C)^-1 * C.T) // (3 x (nbrNum + 1)*3)
		Eigen::MatrixXf M = D * (C.transpose() * C).inverse() * C.transpose();

		// set T
		for (int j = 0; j < nbrNum + 1; ++j) {
			int vtxIndCurr;
			if (j == 0) vtxIndCurr = vtxInd; // j==0 --> self
			else vtxIndCurr = nbrInds[vtxInd][j - 1]; // else --> (j-1)th neighbor

			Triplets::push(vtxInd*3 + 0, vtxIndCurr*3 + 0, M(0,j*3+0)); Triplets::push(vtxInd*3 + 0, vtxIndCurr*3 + 1, M(0,j*3+1)); Triplets::push(vtxInd*3 + 0, vtxIndCurr*3 + 2, M(0,j*3+2));
			Triplets::push(vtxInd*3 + 1, vtxIndCurr*3 + 0, M(1,j*3+0)); Triplets::push(vtxInd*3 + 1, vtxIndCurr*3 + 1, M(1,j*3+1)); Triplets::push(vtxInd*3 + 1, vtxIndCurr*3 + 2, M(1,j*3+2));
			Triplets::push(vtxInd*3 + 2, vtxIndCurr*3 + 0, M(2,j*3+0)); Triplets::push(vtxInd*3 + 2, vtxIndCurr*3 + 1, M(2,j*3+1)); Triplets::push(vtxInd*3 + 2, vtxIndCurr*3 + 2, M(2,j*3+2));
		}
	}
	T = Triplets::createMatrix(vtxNum * 3, vtxNum * 3);

	// set A
	//A = laplacian;
	A = laplacian - T;
	ANormalize = laplacian;

	// set b
	if (mLOCALTRANSFORM) b = Eigen::VectorXf::Zero(vtxNum * 3);
	else b = delta;

	// store members
	mAT = A.transpose();
	mANormalizeT = ANormalize.transpose();
	mb = b;
	mLaplacian = laplacian;
	mDelta = delta;
	mRoiDeltaLengths = roiDeltaLengths;
	mVtxToHandleMap = vtxToHandleMap;

	// set systemMat
	systemMat = (mAT * A).cast<double>();
	systemMatNormalize = (mANormalizeT * ANormalize).cast<double>();

	/*====== Programming Assignment 2 ======*/

	 //Please refer to the following link for the usage of sparse linear system solvers in Eigen
	 //https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	// Do factorization
	if (systemMat.nonZeros() > 0) { // compute systemMat
		mCholeskySolver = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
		mCholeskySolver->compute(systemMat);
		if (mCholeskySolver->info() != Eigen::Success) {
			// Decomposition failed
			std::cout << "Sparse decomposition failed\n";
		}
		else {
			std::cout << "Sparse decomposition succeeded\n";
		}
	}
	if (systemMatNormalize.nonZeros() > 0) { // compute systemMatNormalize
		mCholeskySolverNormalize = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
		mCholeskySolverNormalize->compute(systemMatNormalize);
		if (mCholeskySolverNormalize->info() != Eigen::Success) {
			// Decomposition failed
			std::cout << "Sparse decomposition failed\n";
		}
		else {
			std::cout << "Sparse decomposition succeeded\n";
		}
	}
}

void Deformer::deform() {
	if (mCholeskySolver == nullptr) {
		return;
	}

	/*====== Programming Assignment 2 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* This is the place where the editing techniques 
	/* take place.
	/* Solve for the new vertex positions after the 
	/* specified handles move using the factorized
	/* matrix from Deformer::buildSystemMat(), i.e.,
	/* mCholeskySolver defined in deformer.h
	/**********************************************/

	// Please refer to the following link for the usage of sparse linear system solvers in Eigen
	// https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	std::vector<Vertex*> vertices = mMesh->vertices();
	int vtxNum = vertices.size();
	int roiNum = mRoiList.size();
	int handleNum = mHandleList.size();

	Eigen::VectorXd x; //solution

	// set b
	for (int i = 0; i < handleNum; ++i) { // set handle to new positions
		int vtxInd = mHandleList[i]->index();
		mb[vtxInd * 3 + 0] = mHandleList[i]->position()[0]; // x
		mb[vtxInd * 3 + 1] = mHandleList[i]->position()[1]; // y
		mb[vtxInd * 3 + 2] = mHandleList[i]->position()[2]; // z
	}

	if (mLOCALTRANSFORM) {
		// set b
		for (int i = 0; i < roiNum; ++i) { // set roi to zero
			int vtxInd = mRoiList[i]->index();
			mb[vtxInd * 3 + 0] = 0; // x
			mb[vtxInd * 3 + 1] = 0; // y
			mb[vtxInd * 3 + 2] = 0; // z
		}

		// solve
		x = mCholeskySolver->solve((mAT * mb).cast<double>());

		// compute delta and scale --> set b
		Eigen::VectorXf deltaNew = mLaplacian * x.cast<float>();
		for (int i = 0; i < roiNum; ++i) {
			int vtxInd = mRoiList[i]->index();

			float deltaLengthOld = mRoiDeltaLengths[i];
			float deltaLengthNew = euclidean(deltaNew[vtxInd * 3 + 0], deltaNew[vtxInd * 3 + 1], deltaNew[vtxInd * 3 + 2]);
			float scale = deltaLengthOld / deltaLengthNew;

			mb[vtxInd * 3 + 0] = scale * deltaNew[vtxInd * 3 + 0]; // x
			mb[vtxInd * 3 + 1] = scale * deltaNew[vtxInd * 3 + 1]; // y
			mb[vtxInd * 3 + 2] = scale * deltaNew[vtxInd * 3 + 2]; // z
		}
	}

	// solve for normalization
	x = mCholeskySolverNormalize->solve((mANormalizeT * mb).cast<double>());

	// set new positions for roi
	for (int i = 0; i < roiNum; ++i) {
		int vtxInd = mRoiList[i]->index();
		vertices[vtxInd]->setPosition(Eigen::Vector3f(x[vtxInd*3 + 0], // x
													  x[vtxInd*3 + 1], // y
													  x[vtxInd*3 + 2])); // z
	}

	/*====== Programming Assignment 2 ======*/
}

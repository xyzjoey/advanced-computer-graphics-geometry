#include "deformer.h"
#include "helper.h"
#include <iostream>

Deformer::Deformer() : mMesh(nullptr), mMeshSmooth(nullptr),
                       mCholeskySolver(nullptr), mCholeskySolverSmooth(nullptr) {
}

Deformer::~Deformer() {
	clear();
}

void Deformer::clear() {
	if (mCholeskySolver) {
		delete mCholeskySolver;
	}
	mCholeskySolver = nullptr;
	if (mCholeskySolverSmooth) {
		delete mCholeskySolverSmooth;
	}
	mCholeskySolverSmooth = nullptr;

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

void Deformer::setHandleOffset(int handleFlag, Eigen::Vector3f offset) {
	// move handle position of mMeshSmooth
	int handleNum = mHandleList.size();
	for (int i = 0; i < handleNum; ++i) {
		if (mHandleList[i]->flag() != handleFlag) continue;

		int vtxInd = mHandleList[i]->index();
		Eigen::Vector3f p = mMeshSmooth->vertices()[vtxInd]->position();
		mMeshSmooth->vertices()[vtxInd]->setPosition(p + offset);
	}
}

Eigen::SparseMatrix<float> computeLaplacian(const std::vector<Vertex*>& vertices, const std::vector<int> roiInds, const std::vector<std::vector<int>>& nbrInds, bool useCotWeight) {

	//helper
	auto sum = [](std::vector<double> v) { double sum = 0.0; for (auto& n : v) sum += n; return sum; };
	//helper
	auto cotWeightByInds = [vertices, nbrInds](int i, int j)
	{
		int neighborNum = nbrInds[i].size();
		Vertex* p = vertices[i];
		Vertex* pNbr = vertices[nbrInds[i][j]];
		Vertex* pNbrPrev = vertices[nbrInds[i][(j - 1 + neighborNum) % neighborNum]]; // get previous neighbor
		Vertex* pNbrNext = vertices[nbrInds[i][(j + 1) % neighborNum]]; // get next neighbor

		return cotWeight(p, pNbr, pNbrPrev, pNbrNext);
	};

	int vtxNum = vertices.size();
	int roiNum = roiInds.size();

	Triplets::reset(vtxNum * 3 * 10);

	// self vertices
	for (int i = 0; i < vtxNum * 3; ++i) Triplets::push(i, i, 1); //all one

	// neighbor vertices
	for (int i = 0; i < roiNum; ++i) { // do this only for roi region
		int vtxInd = roiInds[i];
		
		const std::vector<int>& nbrIndsCurr = nbrInds[vtxInd];
		int nbrNum = nbrIndsCurr.size();

		// compute weight for each neighbor
		std::vector<double> weights;
		for (int j = 0; j < nbrNum; ++j) {
			if (useCotWeight) weights.push_back(cotWeightByInds(vtxInd, j));
			else weights.push_back(1.0);
		}

		// set weight
		for (int j = 0; j < nbrNum; ++j) {
			double weight = -weights[j] / sum(weights);
			int nbrInd = nbrIndsCurr[j];
			Triplets::push(vtxInd * 3 + 0, nbrInd * 3 + 0, weight); //x
			Triplets::push(vtxInd * 3 + 1, nbrInd * 3 + 1, weight); //y
			Triplets::push(vtxInd * 3 + 2, nbrInd * 3 + 2, weight); //z
		}
	}

	return Triplets::createMatrix(vtxNum * 3, vtxNum * 3);
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

	// get smoothed mesh
	if (mMeshSmooth != nullptr) delete mMeshSmooth;
	mMeshSmooth = new Mesh(mMesh); // clone mesh
	mMeshSmooth->implicitUmbrellaSmooth(true); // do smoothing

	// get const / member variable
	const std::vector<Vertex*>& vertices = mMesh->vertices();
	const std::vector<Vertex*>& verticesSmooth = mMeshSmooth->vertices();
	const int vtxNum = vertices.size();
	const int roiNum = mRoiList.size();
	const int handleNum = mHandleList.size();

	// initialize
	Eigen::MatrixXf rotations(vtxNum, 3); // alpha, beta, gamma
	std::vector<int> roiInds(roiNum);
	std::vector<std::vector<int>> nbrInds(vtxNum); //neighbor indices

	Eigen::VectorXf pos(vtxNum * 3); // positions
	Eigen::VectorXf posSmooth(vtxNum * 3); // positions of smoothed mesh
	Eigen::VectorXf delta(vtxNum * 3);
	Eigen::VectorXf deltaSmooth(vtxNum * 3);
	Eigen::SparseMatrix<float> laplacian; // vtxNum*3 x vtxNum*3
	Eigen::SparseMatrix<float> laplacianSmooth; // vtxNum*3 x vtxNum*3

	Eigen::SparseMatrix<float> A;
	Eigen::SparseMatrix<float> ASmooth;
	Eigen::VectorXf b;
	Eigen::VectorXf bSmooth;
	Eigen::SparseMatrix<double> systemMat;
	Eigen::SparseMatrix<double> systemMatSmooth;

	// extract rotations from smoothed mesh
	Eigen::Vector3f normal, tangent, binormal, edge;
	for (int i = 0; i < vtxNum; ++i) {
		// normal
		normal = verticesSmooth[i]->normal();
		// tangent
		edge = firstNeighbor(verticesSmooth[i])->position() - verticesSmooth[i]->position();
		tangent = projectPlane(edge, normal).normalized();
		// binormal
		binormal = normal.cross(tangent).normalized();

		// d = original vertex - smoothed vertex
		Eigen::Vector3f d = vertices[i]->position() - verticesSmooth[i]->position();

		// rotation components of d
		rotations(i, 0) = d.dot(normal); // alpha
		rotations(i, 1) = d.dot(tangent); // beta
		rotations(i, 2) = d.dot(binormal); // gamma
	}

	// get roi indices
	for (int i = 0; i < roiNum; ++i) {
		roiInds[i] = mRoiList[i]->index();
	}

	// get neighbor indices
	for (int i = 0; i < vtxNum; ++i) {
		OneRingVertex ring(vertices[i]);
		Vertex* neighbor = nullptr;
		while (neighbor = ring.nextVertex()) nbrInds[i].push_back(neighbor->index());
	}

	// get positions
	for (int i = 0; i < vtxNum; ++i) {
		pos[i * 3 + 0] = vertices[i]->position()[0]; //x
		pos[i * 3 + 1] = vertices[i]->position()[1]; //y
		pos[i * 3 + 2] = vertices[i]->position()[2]; //z
	}
	for (int i = 0; i < vtxNum; ++i) {
		posSmooth[i * 3 + 0] = verticesSmooth[i]->position()[0]; //x
		posSmooth[i * 3 + 1] = verticesSmooth[i]->position()[1]; //y
		posSmooth[i * 3 + 2] = verticesSmooth[i]->position()[2]; //z
	}

	// compute laplacian (with constraints)	
	laplacian = computeLaplacian(vertices, roiInds, nbrInds, true);
	laplacianSmooth = computeLaplacian(verticesSmooth, roiInds, nbrInds, true);

	// compute delta
	delta = laplacian * pos;
	deltaSmooth = laplacianSmooth * posSmooth;

	// set A
	A = laplacian;
	ASmooth = laplacianSmooth;

	// set b
	b = delta;
	bSmooth = deltaSmooth;

	// store members for reuse
	mAT = A.transpose();
	mASmoothT = ASmooth.transpose();
	mb = b;
	mbSmooth = bSmooth;
	mRotations = rotations;

	// set systemMat
	systemMat = (mAT * A).cast<double>();
	systemMatSmooth = (mASmoothT * ASmooth).cast<double>();

	/*====== Programming Assignment 2 ======*/

	 //Please refer to the following link for the usage of sparse linear system solvers in Eigen
	 //https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html

	if (systemMat.nonZeros() > 0) { // compute systemMatNormalize
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
	if (systemMatSmooth.nonZeros() > 0) { // compute systemMatNormalize
		mCholeskySolverSmooth = new Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >();
		mCholeskySolverSmooth->compute(systemMatSmooth);
		if (mCholeskySolverSmooth->info() != Eigen::Success) {
			// Decomposition failed
			std::cout << "Sparse decomposition failed\n";
		}
		else {
			std::cout << "Sparse decomposition succeeded\n";
		}
	}
}

void Deformer::deform(bool rotationInvariant) {
	if (mCholeskySolver == nullptr || mCholeskySolverSmooth == nullptr) {
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

	const std::vector<Vertex*>& vertices = mMesh->vertices();
	const std::vector<Vertex*>& verticesSmooth = mMeshSmooth->vertices();
	const int vtxNum = vertices.size();
	const int roiNum = mRoiList.size();
	const int handleNum = mHandleList.size();

	Eigen::VectorXd x; // to solve

	if (rotationInvariant) {
		// set b
		for (int i = 0; i < handleNum; ++i) { // set handle new positions (of smoothed mesh)
			int vtxInd = mHandleList[i]->index();
			mbSmooth[vtxInd * 3 + 0] = verticesSmooth[vtxInd]->position()[0]; // x
			mbSmooth[vtxInd * 3 + 1] = verticesSmooth[vtxInd]->position()[1]; // y
			mbSmooth[vtxInd * 3 + 2] = verticesSmooth[vtxInd]->position()[2]; // z
		}

		// solve smoothed mesh positions
		x = mCholeskySolverSmooth->solve((mASmoothT * mbSmooth).cast<double>());

		// set new positions to smoothed mesh
		for (int i = 0; i < roiNum; ++i) {
			int vtxInd = mRoiList[i]->index();
			verticesSmooth[vtxInd]->setPosition(Eigen::Vector3f(x[vtxInd * 3 + 0], // x
																x[vtxInd * 3 + 1], // y
																x[vtxInd * 3 + 2])); // z
		}

		// recompute normals
		mMeshSmooth->computeVertexNormals();

		// apply rotations
		Eigen::Vector3f normal, tangent, binormal, edge;
		for (int i = 0; i < roiNum; ++i) {
			int vtxInd = mRoiList[i]->index();

			// normal
			normal = verticesSmooth[vtxInd]->normal();
			// tangent
			edge = firstNeighbor(verticesSmooth[vtxInd])->position() - verticesSmooth[vtxInd]->position();
			tangent = projectPlane(edge, normal).normalized();
			// binormal
			binormal = normal.cross(tangent).normalized();

			// d
			Eigen::Vector3f d = mRotations(vtxInd, 0) * normal + // alpha * normal
				mRotations(vtxInd, 1) * tangent + // beta * tangent
				mRotations(vtxInd, 2) * binormal; // gamma * binormal

			// set final positions
			vertices[vtxInd]->setPosition(verticesSmooth[vtxInd]->position() + d);
		}

	} else {
	
		// set b
		for (int i = 0; i < handleNum; ++i) { // set handle new positions
			int vtxInd = mHandleList[i]->index();
			mb[vtxInd * 3 + 0] = mHandleList[i]->position()[0]; // x
			mb[vtxInd * 3 + 1] = mHandleList[i]->position()[1]; // y
			mb[vtxInd * 3 + 2] = mHandleList[i]->position()[2]; // z
		}

		// solve mesh positions
		x = mCholeskySolver->solve((mAT * mb).cast<double>());
	
		// set positions
		for (int i = 0; i < roiNum; ++i) {
			int vtxInd = mRoiList[i]->index();
			vertices[vtxInd]->setPosition(Eigen::Vector3f(x[vtxInd * 3 + 0], // x
														x[vtxInd * 3 + 1], // y
														x[vtxInd * 3 + 2])); // z
		}
	}
	/*====== Programming Assignment 2 ======*/
}

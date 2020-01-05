#include "mesh.h"
#include "helper.h"
#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <queue>

HEdge::HEdge(bool b) {
	mBoundary = b;

	mTwin = nullptr;
	mPrev = nullptr;
	mNext = nullptr;

	mStart = nullptr;
	mFace = nullptr;

	mFlag = false;
	mValid = true;
}

HEdge* HEdge::twin() const {
	return mTwin;
}

HEdge* HEdge::setTwin(HEdge* e) {
	mTwin = e;
	return mTwin;
}

HEdge* HEdge::prev() const {
	return mPrev;
}

HEdge* HEdge::setPrev(HEdge* e) {
	mPrev = e;
	return mPrev;
}

HEdge* HEdge::next() const {
	return mNext;
}

HEdge* HEdge::setNext(HEdge* e) {
	mNext = e;
	return mNext;
}

Vertex* HEdge::start() const {
	return mStart;
}

Vertex* HEdge::setStart(Vertex* v) {
	mStart = v;
	return mStart;
}

Vertex* HEdge::end() const {
	return mNext->start();
}

Face* HEdge::leftFace() const {
	return mFace;
}

Face* HEdge::setFace(Face* f) {
	mFace = f;
	return mFace;
}

bool HEdge::flag() const {
	return mFlag;
}

bool HEdge::setFlag(bool b) {
	mFlag = b;
	return mFlag;
}

bool HEdge::isBoundary() const {
	return mBoundary;
}

bool HEdge::isValid() const {
	return mValid;
}

bool HEdge::setValid(bool b) {
	mValid = b;
	return mValid;
}

OneRingHEdge::OneRingHEdge(const Vertex* v) {
	if (v == nullptr) {
		mStart = nullptr;
		mNext = nullptr;
	} else {
		mStart = v->halfEdge();
		mNext = v->halfEdge();
	}
}

HEdge* OneRingHEdge::nextHEdge() {
	HEdge* ret = mNext;
	if (mNext != nullptr && mNext->prev()->twin() != mStart) {
		mNext = mNext->prev()->twin();
	} else {
		mNext = nullptr;
	}
	return ret;
}

OneRingVertex::OneRingVertex(const Vertex* v): ring(v) {
}

Vertex* OneRingVertex::nextVertex() {
	HEdge* he = ring.nextHEdge();
	return he != nullptr ? he->end() : nullptr;
}

Vertex::Vertex() : mHEdge(nullptr), mFlag(0) {
	mPosition = Eigen::Vector3f::Zero();
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}

Vertex::Vertex(const Eigen::Vector3f& v): mPosition(v), mHEdge(nullptr), mFlag(0) {
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}

Vertex::Vertex(float x, float y, float z): mHEdge(nullptr), mFlag(0) {
	mPosition = Eigen::Vector3f(x, y, z);
	mColor = VCOLOR_BLUE;
	mNormal = Eigen::Vector3f::Zero();
}


const Eigen::Vector3f& Vertex::position() const {
	return mPosition;
}

const Eigen::Vector3f& Vertex::setPosition(const Eigen::Vector3f& p) {
	mPosition = p;
	return mPosition;
}

const Eigen::Vector3f& Vertex::normal() const {
	return mNormal;
}

const Eigen::Vector3f& Vertex::setNormal(const Eigen::Vector3f& n) {
	mNormal = n;
	return mNormal;
}

const Eigen::Vector3f& Vertex::color() const {
	return mColor;
}

const Eigen::Vector3f& Vertex::setColor(const Eigen::Vector3f& c) {
	mColor = c;
	return mColor;
}

HEdge* Vertex::halfEdge() const {
	return mHEdge;
}

HEdge* Vertex::setHalfEdge(HEdge* he) {
	mHEdge = he;
	return mHEdge;
}

int Vertex::index() const {
	return mIndex;
}

int Vertex::setIndex(int i) {
	mIndex = i;
	return mIndex;
}

int Vertex::flag() const {
	return mFlag;
}

int Vertex::setFlag(int f) {
	mFlag = f;
	return mFlag;
}

bool Vertex::isValid() const {
	return mValid;
}

bool Vertex::setValid(bool b) {
	mValid = b;
	return mValid;
}

bool Vertex::isBoundary() const {
	OneRingHEdge ring(this);
	HEdge* curr = nullptr;
	while (curr = ring.nextHEdge()) {
		if (curr->isBoundary()) {
			return true;
		}
	}
	return false;
}

int Vertex::valence() const {
	int count = 0;
	OneRingVertex ring(this);
	Vertex* curr = nullptr;
	while (curr = ring.nextVertex()) {
		++count;
	}
	return count;
}

Face::Face() : mHEdge(nullptr), mValid(true) {
}

HEdge* Face::halfEdge() const {
	return mHEdge;
}

HEdge* Face::setHalfEdge(HEdge* he) {
	mHEdge = he;
	return mHEdge;
}

bool Face::isBoundary() const {
	HEdge* curr = mHEdge;
	do {
		if (curr->twin()->isBoundary()) {
			return true;
		}
		curr = curr->next();
	} while (curr != mHEdge);
	return false;
}

bool Face::isValid() const {
	return mValid;
}

bool Face::setValid(bool b) {
	mValid = b;
	return mValid;
}

Mesh::Mesh() {
	mVertexPosFlag = true;
	mVertexNormalFlag = true;
	mVertexColorFlag = true;
}

Mesh::~Mesh() {
	clear();
}

Mesh::Mesh(const Mesh* other) {
	if (other == nullptr) return;

	const Eigen::MatrixX3i& F = other->mFaceMat;

	int vtxNum = other->mVertexList.size();
	int faceNum = other->mFaceList.size();

	// vertices
	for (int i = 0; i < vtxNum; ++i) {
		mVertexList.push_back(new Vertex(other->mVertexList[i]->position()));
		mVertexList[i]->setNormal(other->mVertexList[i]->normal());
	}

	// faces
	for (int i = 0; i < faceNum; ++i) {
		addFace(F(i, 0), F(i, 1), F(i, 2));
	}

	std::vector< HEdge* > hedgeList;
	for (int i = 0; i < mBHEdgeList.size(); ++i) {
		if (mBHEdgeList[i]->start()) {
			hedgeList.push_back(mBHEdgeList[i]);
		}
		// TODO
	}
	mBHEdgeList = hedgeList;

	for (int i = 0; i < mVertexList.size(); ++i) {
		mVertexList[i]->adjHEdges.clear();
		mVertexList[i]->setIndex(i);
		mVertexList[i]->setFlag(0);
	}
}

const std::vector< HEdge* >& Mesh::edges() const {
	return mHEdgeList;
}

const std::vector< HEdge* >& Mesh::boundaryEdges() const {
	return mBHEdgeList;
}

const std::vector< Vertex* >& Mesh::vertices() const {
	return mVertexList;
}

const std::vector< Face* >& Mesh::faces() const {
	return mFaceList;
}


bool Mesh::isVertexPosDirty() const {
	return mVertexPosFlag;
}

void Mesh::setVertexPosDirty(bool b) {
	mVertexPosFlag = b;
}

bool Mesh::isVertexNormalDirty() const {
	return mVertexNormalFlag;
}

void Mesh::setVertexNormalDirty(bool b) {
	mVertexNormalFlag = b;
}

bool Mesh::isVertexColorDirty() const {
	return mVertexColorFlag;
}

void Mesh::setVertexColorDirty(bool b) {
	mVertexColorFlag = b;
}

bool Mesh::loadMeshFile(const std::string filename) {
	// Use libigl to parse the mesh file
	bool iglFlag = igl::read_triangle_mesh(filename, mVertexMat, mFaceMat);
	if (iglFlag) {
		clear();

		// Construct the half-edge data structure.
		int numVertices = mVertexMat.rows();
		int numFaces = mFaceMat.rows();

		// Fill in the vertex list
		for (int vidx = 0; vidx < numVertices; ++vidx) {
			mVertexList.push_back(new Vertex(mVertexMat(vidx, 0),
			                                 mVertexMat(vidx, 1),
			                                 mVertexMat(vidx, 2)));
		}
		// Fill in the face list
		for (int fidx = 0; fidx < numFaces; ++fidx) {
			addFace(mFaceMat(fidx, 0), mFaceMat(fidx, 1), mFaceMat(fidx, 2));
		}

		std::vector< HEdge* > hedgeList;
		for (int i = 0; i < mBHEdgeList.size(); ++i) {
			if (mBHEdgeList[i]->start()) {
				hedgeList.push_back(mBHEdgeList[i]);
			}
			// TODO
		}
		mBHEdgeList = hedgeList;

		for (int i = 0; i < mVertexList.size(); ++i) {
			mVertexList[i]->adjHEdges.clear();
			mVertexList[i]->setIndex(i);
			mVertexList[i]->setFlag(0);
		}
	} else {
		std::cout << __FUNCTION__ << ": mesh file loading failed!\n";
	}
	return iglFlag;
}

static void _setPrevNext(HEdge* e1, HEdge* e2) {
	e1->setNext(e2);
	e2->setPrev(e1);
}

static void _setTwin(HEdge* e1, HEdge* e2) {
	e1->setTwin(e2);
	e2->setTwin(e1);
}

static void _setFace(Face* f, HEdge* e) {
	f->setHalfEdge(e);
	e->setFace(f);
}

void Mesh::addFace(int v1, int v2, int v3) {
	Face* face = new Face();

	HEdge* hedge[3];
	HEdge* bhedge[3]; // Boundary half-edges
	Vertex* vert[3];

	for (int i = 0; i < 3; ++i) {
		hedge[i] = new HEdge();
		bhedge[i] = new HEdge(true);
	}
	vert[0] = mVertexList[v1];
	vert[1] = mVertexList[v2];
	vert[2] = mVertexList[v3];

	// Connect prev-next pointers
	for (int i = 0; i < 3; ++i) {
		_setPrevNext(hedge[i], hedge[(i + 1) % 3]);
		_setPrevNext(bhedge[i], bhedge[(i + 1) % 3]);
	}

	// Connect twin pointers
	_setTwin(hedge[0], bhedge[0]);
	_setTwin(hedge[1], bhedge[2]);
	_setTwin(hedge[2], bhedge[1]);

	// Connect start pointers for bhedge
	bhedge[0]->setStart(vert[1]);
	bhedge[1]->setStart(vert[0]);
	bhedge[2]->setStart(vert[2]);
	for (int i = 0; i < 3; ++i) {
		hedge[i]->setStart(vert[i]);
	}

	// Connect start pointers
	// Connect face-hedge pointers
	for (int i = 0; i < 3; ++i) {
		vert[i]->setHalfEdge(hedge[i]);
		vert[i]->adjHEdges.push_back(hedge[i]);
		_setFace(face, hedge[i]);
	}
	vert[0]->adjHEdges.push_back(bhedge[1]);
	vert[1]->adjHEdges.push_back(bhedge[0]);
	vert[2]->adjHEdges.push_back(bhedge[2]);

	// Merge boundary if needed
	for (int i = 0; i < 3; ++i) {
		Vertex* start = bhedge[i]->start();
		Vertex* end = bhedge[i]->end();

		for (int j = 0; j < end->adjHEdges.size(); ++j) {
			HEdge* curr = end->adjHEdges[j];
			if (curr->isBoundary() && curr->end() == start) {
				_setPrevNext(bhedge[i]->prev(), curr->next());
				_setPrevNext(curr->prev(), bhedge[i]->next());
				_setTwin(bhedge[i]->twin(), curr->twin());
				bhedge[i]->setStart(nullptr); // Mark as unused
				curr->setStart(nullptr); // Mark as unused
				break;
			}
		}
	}

	// Finally add hedges and faces to list
	for (int i = 0; i < 3; ++i) {
		mHEdgeList.push_back(hedge[i]);
		mBHEdgeList.push_back(bhedge[i]);
	}
	mFaceList.push_back(face);
}

Eigen::Vector3f Mesh::initBboxMin() const {
	return (mVertexMat.colwise().minCoeff()).transpose();
}

Eigen::Vector3f Mesh::initBboxMax() const {
	return (mVertexMat.colwise().maxCoeff()).transpose();
}

void Mesh::groupingVertexFlags() {
	// Init to 255
	for (Vertex* vert : mVertexList) {
		if (vert->flag() != 0) {
			vert->setFlag(255);
		}
	}
	// Group handles
	int id = 0;
	std::vector< Vertex* > tmpList;
	for (Vertex* vert : mVertexList) {
		if (vert->flag() == 255) {
			++id;
			vert->setFlag(id);

			// Do search
			tmpList.push_back(vert);
			while (!tmpList.empty()) {
				Vertex* v = tmpList.back();
				tmpList.pop_back();

				OneRingVertex orv = OneRingVertex(v);
				while (Vertex* v2 = orv.nextVertex()) {
					if (v2->flag() == 255) {
						v2->setFlag(id);
						tmpList.push_back(v2);
					}
				}
			}
		}
	}
}

void Mesh::clear() {
	for (int i = 0; i < mHEdgeList.size(); ++i) {
		delete mHEdgeList[i];
	}
	for (int i = 0; i < mBHEdgeList.size(); ++i) {
		delete mBHEdgeList[i];
	}
	for (int i = 0; i < mVertexList.size(); ++i) {
		delete mVertexList[i];
	}
	for (int i = 0; i < mFaceList.size(); ++i) {
		delete mFaceList[i];
	}

	mHEdgeList.clear();
	mBHEdgeList.clear();
	mVertexList.clear();
	mFaceList.clear();
}

std::vector< int > Mesh::collectMeshStats() {
	int V = 0; // # of vertices
	int E = 0; // # of half-edges
	int F = 0; // # of faces
	int B = 0; // # of boundary loops
	int C = 0; // # of connected components
	int G = 0; // # of genus

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Collect mesh information as listed above.
	/**********************************************/

	V = mVertexList.size();
	E = mHEdgeList.size();
	F = mFaceList.size();
	B = countBoundaryLoops();
	C = countConnectedComponents();
	G = C - (V - E / 2 + F + B) / 2;

	/*====== Programming Assignment 0 ======*/

	std::vector< int > stats;
	stats.push_back(V);
	stats.push_back(E);
	stats.push_back(F);
	stats.push_back(B);
	stats.push_back(C);
	stats.push_back(G);
	return stats;
}

int Mesh::countBoundaryLoops() {
	int count = 0;

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Helper function for Mesh::collectMeshStats()
	/**********************************************/

	//const, helper
	enum { INVISITED, VISITED };
	auto nextBHedge = [](HEdge* edge)
	{
		OneRingHEdge ring(edge->end());
		HEdge* curr = ring.nextHEdge();

		while (!curr->isBoundary() || curr == edge->twin()) curr = ring.nextHEdge();
		return curr;
	};

	//initialize
	for (int i = 0; i < mBHEdgeList.size(); ++i) mBHEdgeList[i]->setFlag(INVISITED);

	//check each boundary edge
	for (int i = 0; i < mBHEdgeList.size(); ++i) {
		if (mBHEdgeList[i]->flag() == VISITED) continue;

		HEdge* startEdge = mBHEdgeList[i];
		HEdge* currEdge = startEdge;

		//walk through a boundary
		do {
			currEdge->setFlag(VISITED);
			currEdge->twin()->setFlag(VISITED);
			currEdge = nextBHedge(currEdge);
		} while (currEdge && currEdge != startEdge);

		++count;
	}

	/*====== Programming Assignment 0 ======*/

	return count;
}

int Mesh::countConnectedComponents() {
	int count = 0;

	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Helper function for Mesh::collectMeshStats()
	/* Count the number of connected components of
	/* the mesh. (Hint: use a stack)
	/**********************************************/

	//const, helper
	enum { INVISITED, VISITED };
	auto setFlag = [](Vertex* vtx, bool flag) { vtx->halfEdge()->setFlag(flag); };
	auto isVisited = [](Vertex* vtx) { return vtx->halfEdge()->flag() == VISITED; };
	int i = 0;
	auto nextInvisited = [&, this]() -> Vertex*
	{
		for (; i < mVertexList.size(); ++i) if (!isVisited(mVertexList[i])) return mVertexList[i];
		return nullptr;
	};

	//initialize
	for (int i = 0; i < mVertexList.size(); ++i) setFlag(mVertexList[i], INVISITED);

	Vertex* start = nullptr;
	std::queue<Vertex*> q;

	//count each connected component
	while (start = nextInvisited()) {

		//BFS
		q.push(start);
		while (!q.empty()) {
			//pop from queue
			Vertex* curr = q.front();
			q.pop();

			if (isVisited(curr)) continue;

			//push neighboring vertices to queue
			Vertex* neighbor = nullptr;
			OneRingVertex ring(curr);
			while (neighbor = ring.nextVertex()) if (!isVisited(neighbor)) q.push(neighbor);

			setFlag(curr, VISITED);
		}

		++count;
	}

	/*====== Programming Assignment 0 ======*/

	return count;
}

void Mesh::computeVertexNormals() {
	/*====== Programming Assignment 0 ======*/

	/**********************************************/
	/*          Insert your code here.            */
	/**********************************************/
	/*
	/* Compute per-vertex normal using neighboring
	/* facet information. (Hint: remember using a 
	/* weighting scheme. Plus, do you notice any
	/* disadvantages of your weighting scheme?)
	/**********************************************/

	/*====== Programming Assignment 0 ======*/
	
	//helper
	auto vector = [](Vertex* start, Vertex* end) { return end->position() - start->position(); };

	//for each vertex
	for (int i = 0; i < mVertexList.size(); ++i) {
		//initialize
		Eigen::Vector3f total(0, 0, 0);

		Vertex* vtx = mVertexList[i];
		OneRingVertex ring(vtx);

		Vertex* start = ring.nextVertex();
		Vertex* vtx1 = nullptr;
		Vertex* vtx2 = start;

		//traverse neighbor vertices and accumulate cross product
		while (vtx2) {
			vtx1 = vtx2;
			vtx2 = ring.nextVertex();
			total = total + vector(vtx, vtx1).cross(vector(vtx, vtx2 ? vtx2 : start));
		}

		vtx->setNormal(total.normalized());
	}


	// Notify mesh shaders
	setVertexNormalDirty(true);
}


void Mesh::umbrellaSmooth(bool cotangentWeights) {
	/*====== Programming Assignment 1 ======*/

	const float LAMBDA = 1;
	const int vtxNum = mVertexList.size();

	std::vector<std::vector<int>> neighborInds(vtxNum);
	Eigen::MatrixXf P(vtxNum, 3);
	Eigen::MatrixXf Pnew;
	Eigen::SparseMatrix<float> L;

	//Prepare neighborInds (indices of neighboring vertices)
	for (int i = 0; i < vtxNum; ++i) {
		OneRingVertex ring(mVertexList[i]);
		Vertex* neighbor = nullptr;
		while (neighbor = ring.nextVertex()) neighborInds[i].push_back(neighbor->index());
	}

	//Prepare P (position matrix)
	for (int i = 0; i < vtxNum; ++i) {
		P(i, 0) = mVertexList[i]->position()[0];
		P(i, 1) = mVertexList[i]->position()[1];
		P(i, 2) = mVertexList[i]->position()[2];
	}

	if (cotangentWeights) {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 1: Implement the cotangent weighting
		/* scheme for explicit mesh smoothing.
		/*
		/* Hint:
		/* It is advised to double type to store the
		/* weights to avoid numerical issues.
		/**********************************************/

		//helper
		auto sum = [](std::vector<double> v) { double sum = 0.0; for (auto& n : v) sum += n; return sum; };
		auto weight = [this, neighborInds](int i, int j)
		{
			int neighborNum = neighborInds[i].size();
			Vertex* p = mVertexList[i];
			Vertex* pNbr = mVertexList[neighborInds[i][j]];
			Vertex* pNbrPrev = mVertexList[neighborInds[i][(j - 1 + neighborNum) % neighborNum]];
			Vertex* pNbrNext = mVertexList[neighborInds[i][(j + 1) % neighborNum]];

			return cotWeight(p, pNbr, pNbrPrev, pNbrNext);
		};

		//Prepare L
		Triplets::reset(vtxNum * 10);
		for (int i = 0; i < vtxNum; ++i) {
			Triplets::push(i, i, -1);
			std::vector<double> weights;
			for (int j = 0; j < neighborInds[i].size(); ++j) weights.push_back(weight(i, j));
			for (int j = 0; j < neighborInds[i].size(); ++j) Triplets::push(i, neighborInds[i][j], weights[j] / sum(weights));
		}
		L = Triplets::createMatrix(vtxNum, vtxNum);

	}
	else {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 2: Implement the uniform weighting
		/* scheme for explicit mesh smoothing.
		/**********************************************/

		//Prepare L
		Triplets::reset(vtxNum * 10);
		for (int i = 0; i < vtxNum; ++i) {
			Triplets::push(i, i, -1);
			double weight = 1.0 / neighborInds[i].size();
			for (int j = 0; j < neighborInds[i].size(); ++j) Triplets::push(i, neighborInds[i][j], weight);
		}
		L = Triplets::createMatrix(vtxNum, vtxNum);
	}

	//Compute Pnew
	Pnew = (sparseI(vtxNum) + LAMBDA * L) * P;
	//Assign Pnew to vertices
	for (int i = 0; i < vtxNum; ++i) mVertexList[i]->setPosition(Eigen::Vector3f(Pnew(i, 0), Pnew(i, 1), Pnew(i, 2)));


	/*====== Programming Assignment 1 ======*/

	computeVertexNormals();
	// Notify mesh shaders
	setVertexPosDirty(true);
	setVertexNormalDirty(true);
}

void Mesh::implicitUmbrellaSmooth(bool cotangentWeights) {
	/*====== Programming Assignment 1 ======*/

	/* A sparse linear system Ax=b solver using the conjugate gradient method. */
	auto fnConjugateGradient = [](const Eigen::SparseMatrix< float >& A,
	                              const Eigen::VectorXf& b,
	                              int maxIterations,
	                              float errorTolerance,
	                              Eigen::VectorXf& x)
	{
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Params:
		/*  A: 
		/*  b: 
		/*  maxIterations:	Max number of iterations
		/*  errorTolerance: Error tolerance for the early stopping condition
		/*  x:				Stores the final solution, but should be initialized. 
		/**********************************************/
		/*
		/* Step 1: Implement the biconjugate gradient
		/* method.
		/* Hint: https://en.wikipedia.org/wiki/Biconjugate_gradient_method
		/**********************************************/

		Eigen::VectorXf x1, x2, r1, r2, r1New, r2New, p1, p2;
		float alpha, beta;

		x1 = Eigen::VectorXf::Ones(b.size());
		r1 = b - A * x1;
		r2 = r1;
		p1 = r1;
		p2 = r2;

		double bNorm = b.norm();

		for (int i = 0; i < maxIterations; ++i) {
			alpha = (r2.transpose() * r1 / (p2.transpose() * A * p1))(0, 0);
			x1 = x1 + alpha * p1;
			x2 = x2 + alpha * p2;
			r1New = r1 - alpha * A * p1;
			r2New = r2 - alpha * A.transpose() * p2;
			beta = (r2New.transpose() * r1New / (r2.transpose() * r1))(0, 0);
			p1 = r1New + beta * p1;
			p2 = r2New + beta * p2;

			r1 = r1New;
			r2 = r2New;

			if (r1.norm() / bNorm <= errorTolerance) break;
		}

		x = x1;
	};
	
	/* IMPORTANT:
	/* Please refer to the following link about the sparse matrix construction in Eigen. */
	/* http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title3 */

	// constant
	const float LAMBDA = 10;
	const int MAX_ITER = 20;
	const float ERROR_TOL = 1e-6;
	const int vtxNum = mVertexList.size();

	// initialize
	std::vector<std::vector<int>> neighborInds(vtxNum);
	Eigen::VectorXf Px(vtxNum), Py(vtxNum), Pz(vtxNum);
	Eigen::VectorXf PxNew, PyNew, PzNew;
	Eigen::SparseMatrix<float> L;

	// get neighborInds (indices of neighboring vertices)
	for (int i = 0; i < vtxNum; ++i) {
		OneRingVertex ring(mVertexList[i]);
		Vertex* neighbor = nullptr;
		while (neighbor = ring.nextVertex()) neighborInds[i].push_back(neighbor->index());
	}

	// set Px Py Pz (position matrix)
	for (int i = 0; i < vtxNum; ++i) Px(i) = mVertexList[i]->position()[0];
	for (int i = 0; i < vtxNum; ++i) Py(i) = mVertexList[i]->position()[1];
	for (int i = 0; i < vtxNum; ++i) Pz(i) = mVertexList[i]->position()[2];

	// compute L
	if (cotangentWeights) {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 2: Implement the cotangent weighting 
		/* scheme for implicit mesh smoothing. Use
		/* the above fnConjugateGradient for solving
		/* sparse linear systems.
		/*
		/* Hint:
		/* It is advised to double type to store the
		/* weights to avoid numerical issues.
		/**********************************************/

		//helper
		auto sum = [](std::vector<double> v) { double sum = 0.0; for (auto& n : v) sum += n; return sum; };
		auto weight = [this, neighborInds](int i, int j)
		{
			int neighborNum = neighborInds[i].size();
			Vertex* p = mVertexList[i];
			Vertex* pNbr = mVertexList[neighborInds[i][j]];
			Vertex* pNbrPrev = mVertexList[neighborInds[i][(j - 1 + neighborNum) % neighborNum]];
			Vertex* pNbrNext = mVertexList[neighborInds[i][(j + 1) % neighborNum]];

			return cotWeight(p, pNbr, pNbrPrev, pNbrNext);
		};

		// set L
		Triplets::reset(vtxNum * 10);
		for (int i = 0; i < vtxNum; ++i) {
			Triplets::push(i, i, -1);
			std::vector<double> weights;
			for (int j = 0; j < neighborInds[i].size(); ++j) weights.push_back(weight(i, j));
			for (int j = 0; j < neighborInds[i].size(); ++j) Triplets::push(i, neighborInds[i][j], weights[j] / sum(weights));
		}
		L = Triplets::createMatrix(vtxNum, vtxNum);

	} else {
		/**********************************************/
		/*          Insert your code here.            */
		/**********************************************/
		/*
		/* Step 3: Implement the uniform weighting 
		/* scheme for implicit mesh smoothing. Use
		/* the above fnConjugateGradient for solving
		/* sparse linear systems.
		/**********************************************/

		// set L
		Triplets::reset(vtxNum * 10);
		for (int i = 0; i < vtxNum; ++i) {
			Triplets::push(i, i, -1);
			double weight = 1.0 / neighborInds[i].size();
			for (int j = 0; j < neighborInds[i].size(); ++j) Triplets::push(i, neighborInds[i][j], weight);
		}
		L = Triplets::createMatrix(vtxNum, vtxNum);
	}

	// solve PxNew, PyNew, PzNew
	Eigen::SparseMatrix<float> A = sparseI(vtxNum) - LAMBDA * L;
	fnConjugateGradient(A, Px, MAX_ITER, ERROR_TOL, PxNew);
	fnConjugateGradient(A, Py, MAX_ITER, ERROR_TOL, PyNew);
	fnConjugateGradient(A, Pz, MAX_ITER, ERROR_TOL, PzNew);

	// assign PxNew, PyNew, PzNew to vertices
	for (int i = 0; i < vtxNum; ++i) mVertexList[i]->setPosition(Eigen::Vector3f(PxNew[i], PyNew[i], PzNew[i]));

	/*====== Programming Assignment 1 ======*/

	computeVertexNormals();
	// Notify mesh shaders
	setVertexPosDirty(true);
	setVertexNormalDirty(true);
}

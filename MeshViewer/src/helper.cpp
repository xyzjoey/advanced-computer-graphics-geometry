#include "helper.h"

std::vector<Eigen::Triplet<double>> Triplets::triplets; // define static triplets

Eigen::SparseMatrix<float> sparseI(int size)
{
	Triplets::reset(size);
	for (int i = 0; i < size; ++i) Triplets::push(i, i, 1);
	return Triplets::createMatrix(size, size);
}

double cotWeight(Vertex* p, Vertex* pNbr, Vertex* pNbrPrev, Vertex* pNbrNext)
{
	//Nbr: neighbor

	//const, helper
	const double MAX_WEIGHT = 1 / tan(1e-6);
	auto angle = [](Eigen::Vector3f v1, Eigen::Vector3f v2) { return acos(v1.dot(v2) / (v1.norm() * v2.norm())); };
	auto vector = [](Vertex* start, Vertex* end) { return end->position() - start->position(); };

	//Compute weight
	double angle1 = angle(vector(pNbrPrev, p), vector(pNbrPrev, pNbr));
	double angle2 = angle(vector(pNbrNext, p), vector(pNbrNext, pNbr));
	double weight = 1 / tan(angle1) + 1 / tan(angle2);

	//Handle invalid value
	if (isnan(weight)) weight = 0.0;
	if (weight > MAX_WEIGHT) weight = MAX_WEIGHT;
	return abs(weight);
}

Vertex* firstNeighbor(Vertex* v) {
	OneRingVertex ring(v);
	return ring.nextVertex();
}

double euclidean(double x, double y, double z) {
	return sqrtf(x*x + y*y + z*z);
}

Eigen::Vector3f project(Eigen::Vector3f v, Eigen::Vector3f u) {
	// u is an unit vector
	// project v onto u
	return v.dot(u) * u;
}

Eigen::Vector3f projectPlane(Eigen::Vector3f v, Eigen::Vector3f normal) {
	// project v to the plane with normal
	return v - project(v, normal);
}
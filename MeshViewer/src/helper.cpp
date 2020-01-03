#include "mesh.h"
#include "helper.h"

std::vector<Eigen::Triplet<double>> Triplets::triplets;

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

	//Compute
	double angle1 = angle(vector(pNbrPrev, p), vector(pNbrPrev, pNbr));
	double angle2 = angle(vector(pNbrNext, p), vector(pNbrNext, pNbr));
	double weight = 1 / tan(angle1) + 1 / tan(angle2);

	//Handle invalid value
	if (isnan(weight)) weight = 0.0;
	if (weight > MAX_WEIGHT) weight = MAX_WEIGHT;
	return abs(weight);
}

float euclidean(float x, float y, float z) {
	return sqrtf(x*x + y*y + z*z);
}
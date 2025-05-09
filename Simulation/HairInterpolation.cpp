#include "Simulation/HairInterpolation.h"

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;
using namespace HairInterpolation;

void RenderedHair::timeStep() {
	// Draw simulation model
	SimulationModel *model = Simulation::getCurrent()->getModel();
	const ParticleData &pd = model->getParticles();
	const OrientationData &od = model->getOrientations();

	// COMPUTE BARYCENTRIC WEIGHTS FOR EACH GUIDE HAIR'S ROOT
	// FOR POSITIONING THE RENDERED HAIR

	// interpolated rest points of rendered hair
	vector<Vector3r> interpRestPoints = {root};
	// barycentric weights of each guide hair for positioning rendered hairs
	// weights[g] = weight of guide hair g
	vector<Real> rootWeights = {0, 0, 0};
	
	// store roots of all guide hairs
	vector<Vector3r> roots = {};
	for (unsigned int i = 0; i < model->getLineModels().size(); i++)
	{
		LineModel *lineModel = model->getLineModels()[i];
		
		const unsigned int indexOffset = lineModel->getIndexOffset();
		const unsigned int rootIndex = lineModel->getEdges()[0].m_vert[0] + indexOffset;
		const Vector3r &root = pd.getPosition(rootIndex);
		roots.push_back(root);
	}

	// calculate denominator in weights
	const Real lambdaDenom = (roots[0] - roots[2]).cross(roots[1] - roots[2]).norm();

	// calculate barycentric weights
	for (unsigned int i = 0; i < model->getLineModels().size(); i++)
	{
		LineModel *lineModel = model->getLineModels()[i];
		const Vector3r prevRoot = roots[(i + 3 - 1) % 3];
		const Vector3r currRoot = roots[(i + 3 + 0) % 3];
		const Vector3r nextRoot = roots[(i + 3 + 1) % 3];
		const Real lambda = (root - prevRoot).cross(nextRoot - prevRoot).norm() / lambdaDenom;
		rootWeights[i] = lambda;
	}

	// simulated forces on each rendered hair vertex
	// forceSims[i] = simulated force for rendered hair vertex i
	vector<Vector3r> forceSims(currState->positions.size());
	for (unsigned int i = 0; i < model->getLineModels().size(); i++)
	{
		LineModel *lineModel = model->getLineModels()[i];
		for(unsigned int e=0; e<lineModel->getEdges().size(); e++)
		{
			const unsigned int indexOffset = lineModel->getIndexOffset();
			const unsigned int indexOffsetQuaternions = lineModel->getIndexOffsetQuaternions();
			const unsigned int i1 = lineModel->getEdges()[e].m_vert[0] + indexOffset;
			const unsigned int i2 = lineModel->getEdges()[e].m_vert[1] + indexOffset;
			const unsigned int iq = lineModel->getEdges()[e].m_quat + indexOffsetQuaternions;
			const Vector3r &v1 = pd.getPosition(i1);
			const Vector3r &v2 = pd.getPosition(i2);
			const Quaternionr &q = od.getQuaternion(iq);

			// d_j,3 in equation 4
			Vector3r tangentVec = computeSegmentTangentVector(q);
			// l_j in equation 4
			Real length = (pd.getPosition0(i2) - pd.getPosition0(i1)).norm();
			// k^SS in equation 4
			Real stiffness = model->getRodStretchingStiffness();
			const Vector3r force = -stiffness / length * ((v2 - v1) / length - tangentVec);
			if(i == 0) {
				forceSims[e] = rootWeights[i] * force;
			} else {
				forceSims[e] += rootWeights[i] * force;
			}
		}
	}
	
	updateRenderedHair(forceSims);
}

// compute rest positions and quaternions for each
// vertex/segment of the rendered hair rooted at root
void RenderedHair::init() {
	SimulationModel *model = Simulation::getCurrent()->getModel();
	const ParticleData &pd = model->getParticles();
	const OrientationData &od = model->getOrientations();
	
	// get number of vertices in rendered hair
	// should just be length of shortest guide hair
	const int minNumEdges = std::min(
		model->getLineModels()[0]->getEdges().size(),
		std::min(
			model->getLineModels()[1]->getEdges().size(),
			model->getLineModels()[2]->getEdges().size()
		)
	);

	// COMPUTE BARYCENTRIC WEIGHTS FOR EACH GUIDE HAIR'S ROOT
	// FOR POSITIONING THE RENDERED HAIR

	// interpolated rest positions/orientations of rendered hair
	vector<Vector3r> restPositions(minNumEdges + 1);
	vector<Quaternionr> restOrientations(minNumEdges);
	restPositions[0] = root;
	// barycentric weights of each guide hair for positioning rendered hairs
	// rootWeights[g] = weight of guide hair g
	vector<Real> rootWeights = {0, 0, 0};
	
	// grab each root
	vector<Vector3r> roots = {};
	for (unsigned int i = 0; i < model->getLineModels().size(); i++)
	{
		LineModel *lineModel = model->getLineModels()[i];
		
		const unsigned int indexOffset = lineModel->getIndexOffset();
		const unsigned int rootIndex = lineModel->getEdges()[0].m_vert[0] + indexOffset;
		const Vector3r &root = pd.getPosition(rootIndex);
		roots.push_back(root);
	}

	// calculate denominator in weights
	const Real lambdaDenom = (roots[0] - roots[2]).cross(roots[1] - roots[2]).norm();

	// calculate BC weight from each root
	for (unsigned int i = 0; i < model->getLineModels().size(); i++)
	{
		LineModel *lineModel = model->getLineModels()[i];
		const Vector3r prevRoot = roots[(i + 3 - 1) % 3];
		const Vector3r currRoot = roots[(i + 3 + 0) % 3];
		const Vector3r nextRoot = roots[(i + 3 + 1) % 3];
		const Real lambda = (root - prevRoot).cross(nextRoot - prevRoot).norm() / lambdaDenom;
		rootWeights[i] = lambda;
	}

	// compute interpolated rest positions
	for(unsigned int e=0; e<minNumEdges; e++) {
		Vector3r newPoint(0, 0, 0);
		for (unsigned int i = 0; i < model->getLineModels().size(); i++) {
			LineModel *lineModel = model->getLineModels()[i];
			const unsigned int indexOffset = lineModel->getIndexOffset();
			// root vertex is not interpolated, so only look at v2 of each edge
			const unsigned int i2 = lineModel->getEdges()[e].m_vert[1] + indexOffset;
			const Vector3r &v2 = pd.getPosition(i2);

			// ADD TO INTERPOLATED POINTS
			// root vertex is not interpolated, so only look at v2 of each edge
			newPoint += rootWeights[i] * v2;
		}
		restPositions[e+1] = newPoint;
	}

	// compute rest quaternions based on interpolated positions
	Vector3r from(0, 0, 1);
	for(int i=0; i<minNumEdges; i++)	
	{
		Vector3r to = (restPositions[i + 1] - restPositions[i]).normalized();
		Quaternionr dq = Quaternionr::FromTwoVectors(from, to);
		if(i == 0) restOrientations[i] = dq;
		else restOrientations[i] = dq * restOrientations[i - 1];
		from = to;
	}

	restState->positions = restPositions;
	restState->orientations = restOrientations;
	currState->positions = restPositions;
	currState->orientations = restOrientations;
}

// ================================================================================
// ================================ RECONSTRUCTION ================================
// ================================================================================

// computes the quaternion w^0_j = (q^0_j)^-1 * q^0_(j+1)
// which represents the scaled darboux vector of edge j at rest
// essentially the angular velocity of edge j
Quaternionr RenderedHair::computeDarbouxVector(int j) {
	return restState->orientations[j].conjugate() * restState->orientations[j+1];
}

// EQUATION 3
// computes the phi value (either +1 or -1)
// for the j-th segment of the rendered hair
int RenderedHair::computePhiScalar(int j) {
	// q_j in the inequality
	Quaternionr& q_j = currState->orientations[j];
	// q_(j+1) in the inequality
	Quaternionr& q_jp1 = currState->orientations[j+1];
	// q_j^0 in the inequality
	// NOTE: THIS IS THE DARBOUX VECTOR QUATERNION, NOT THE J-TH REST QUATERNION
	Quaternionr darbouxVec = computeDarbouxVector(j);
	Real leftHandSide = addQuats(q_j.conjugate() * q_jp1, scaleQuat(darbouxVec, -1)).squaredNorm();
	Real rightHandSide = addQuats(q_j.conjugate() * q_jp1, darbouxVec).squaredNorm();
	return leftHandSide <= rightHandSide ? +1 : -1;
}

// EQUATION 14
// computes the b value for the i-th segment of the rendered hair
Quaternionr RenderedHair::computeBQuat(int i) {
	SimulationModel *model = Simulation::getCurrent()->getModel();
	// k^bt in the equation (bend and twisting stiffness)
	Real kBT = getBendTwistStiffness();
	// phi_(i-1) in the equation
	int phi_im1 = computePhiScalar(i - 1);
	// phi_i in the equation
	int phi_i = computePhiScalar(i);
	// q_(i-1) in the equation
	Quaternionr q_im1 = currState->orientations[i - 1];
	// q_(i-1)^0 in the equation
	// Quaternionr q_im1Rest = renderedHairOrientationsRest[i - 1];
	Quaternionr w_im1 = computeDarbouxVector(i - 1);
	// q_i^0 in the equation
	// Quaternionr q_iRest = renderedHairOrientationsRest[i];
	Quaternionr w_i = computeDarbouxVector(i);
	// q_(i+1) in the equation
	Quaternionr q_ip1 = currState->orientations[i + 1];

	// b = -k^BT(
	//		phi_(i-1) * q_(i-1) * w^0_(i-1)
	//		+ phi_i * q_(i+1) * wbar^0_i
	// )
	return scaleQuat(addQuats(scaleQuat(q_im1 * w_im1, phi_im1), scaleQuat(q_ip1 * w_i.conjugate(), phi_i)), -kBT);
}

// EQUATION 17
// computes the lambda value for the i-th segment of the rendered hair
// using the provided simulated net force
Real RenderedHair::computeLambdaScalar(int i, Vector3r netForce) {
	Real l_i = (restState->positions[i + 1] - restState->positions[i]).norm();
	Quaternionr b_i = computeBQuat(i);

	return (2 * l_i * netForce).norm() + b_i.norm();
}

// THEOREM 4.1
// computes inverted matrix (A_i - lambda_i * I_4)^-1 for segment i
Matrix4d RenderedHair::computeInvertedMatrix(int i, Vector3r netForce) {
	Real l_i = (restState->positions[i + 1] - restState->positions[i]).norm();
	Vector3r v = -2 * l_i * netForce;
	Real lambda = computeLambdaScalar(i, netForce);

	Real frontConstant = 1 / (v.squaredNorm() - lambda * lambda);
	Matrix4d matrix({
		{lambda + v.z(),	0,				-v.x(),			 v.y()			},
		{0,					lambda + v.z(),	-v.y(),			-v.x()			},
		{-v.x(),			-v.y(),			lambda - v.z(),	0				},
		{ v.y(),			-v.x(),			0,				lambda - v.z()	}
	});
	return frontConstant * matrix.transpose();
}

// EQUATION 17
// computes the new orientation of segment i
Quaternionr RenderedHair::computeNewOrientation(int i, Vector3r netForce) {
	Matrix4d invertedMatrix = computeInvertedMatrix(i, netForce);
	Quaternionr b = computeBQuat(i);
	Vector4r b_vec(b.w(), b.x(), b.y(), b.z());
	
	Vector4r vectorOutput = invertedMatrix * b_vec;
	Quaternionr quatOutput(vectorOutput.x(), vectorOutput.y(), vectorOutput.z(), vectorOutput.w());
	quatOutput.normalize();
	return quatOutput;
}

// calculate the tangent axis vector for some quaternion orientation
Vector3r RenderedHair::computeSegmentTangentVector(Quaternionr q) {
	Quaternionr axis(0, 0, 0, 1);
	Quaternionr product = q * axis * q.conjugate();
	Vector3r vec_product(product.x(), product.y(), product.z());
	return vec_product;
}


// EQUATION 12
// recursively calculate each new vertex position for the rendered hair
// netForces should be 0-indexed such that netForces[0] is the net force
// on the root vertex (thus, zero force) and netForces[1] is the net force
// on the first actually free vertex
void RenderedHair::updateRenderedHair(vector<Vector3r> netForces) {
	SimulationModel *model = Simulation::getCurrent()->getModel();

	vector<Vector3r> newPositions(currState->positions.size());
	vector<Quaternionr> newOrientations(currState->orientations.size());
	// set root to start
	newPositions[0] = root;

	// k^ss in the equation
	Real kSS = getStretchShearStiffness();
	// current position for recursion
	Vector3r curPos = root;
	for(int m = 0; m < netForces.size() - 1; m++) {
		Real l_m = (restState->positions[m + 1] - restState->positions[m]).norm();
		Vector3r fSIM_m = netForces[m];
		Quaternionr q_m = computeNewOrientation(m, fSIM_m);
		Vector3r d_m3 = computeSegmentTangentVector(q_m);
		curPos += l_m * ((fSIM_m * l_m) / kSS + d_m3);
		newPositions[m+1] = curPos;
		newOrientations[m] = q_m;
	}
	currState->positions = newPositions;
	currState->orientations = newOrientations;
}

// returns the magnitude of the stretching and shearing stiffness vector
Real HairInterpolation::getStretchShearStiffness() {
	SimulationModel *model = Simulation::getCurrent()->getModel();
	// VECTOR MAGNITUDE
	return pow(
		pow(model->getRodStretchingStiffness(), 2) +
		pow(model->getRodShearingStiffnessX(), 2) +
		pow(model->getRodShearingStiffnessY(), 2),
		0.5);
	
	// // ARITHMETIC MEAN
	// return (
	// 	model->getRodStretchingStiffness() +
	// 	model->getRodShearingStiffnessX() +
	// 	model->getRodShearingStiffnessY()
	// 	) / 3.0;

	// // GEOMETRIC MEAN
	// return pow(
	// 	model->getRodStretchingStiffness() *
	// 	model->getRodShearingStiffnessX() *
	// 	model->getRodShearingStiffnessY(),
	// 	1.0 / 3.0);
}

// returns the magnitude of the bending and twisting stiffness vector
Real HairInterpolation::getBendTwistStiffness() {
	SimulationModel *model = Simulation::getCurrent()->getModel();
	// VECTOR MAGNITUDE
	return pow(
		pow(model->getRodTwistingStiffness(), 2) +
		pow(model->getRodBendingStiffnessX(), 2) +
		pow(model->getRodBendingStiffnessY(), 2),
		0.5);

	// // ARITHMETIC MEAN
	// return (
	// 	model->getRodTwistingStiffness() +
	// 	model->getRodBendingStiffnessX() +
	// 	model->getRodBendingStiffnessY()
	// 	) / 3.0;

	// // GEOMETRIC MEAN
	// return pow(
	// 	model->getRodTwistingStiffness() *
	// 	model->getRodBendingStiffnessX() *
	// 	model->getRodBendingStiffnessY(),
	// 	1.0 / 3.0);
}

// adds two quaternions together component-wise
Quaternionr HairInterpolation::addQuats(Quaternionr q1, Quaternionr q2) {
	return Quaternionr(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}

// scales each component of q by s
Quaternionr HairInterpolation::scaleQuat(Quaternionr q, Real s) {
	return Quaternionr(q.w() * s, q.x() * s, q.y() * s, q.z() * s);
}

// returns string representation of vector v as "(v.x, v.y, v.z)"
std::string HairInterpolation::vecToString(Vector3r v) {
	return "(" + std::to_string(v.x()) + ", " + std::to_string(v.y()) + ", " + std::to_string(v.z()) + ")";
}
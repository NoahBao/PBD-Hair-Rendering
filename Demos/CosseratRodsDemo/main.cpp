#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Demos/Visualization/Selection.h"
#include "Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "Simulation/SimulationModel.h"
#include "Simulation/TimeStepController.h"
#include <iostream>
#include "Demos/Visualization/Visualization.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"
#include "Utils/FileSystem.h"
#include "Demos/Common/DemoBase.h"
#include "Simulation/Simulation.h"
#include "../Common/imguiParameters.h"
#include "Simulation/HairInterpolation.h"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
	#define new DEBUG_NEW 
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;
using namespace HairInterpolation;

void timeStep ();
void buildModel ();
void createHelix(const Vector3r &position, const Matrix3r &orientation, Real radius, Real height, Real totalAngle, int nPoints);
void render ();
void reset();

// list of all rendered hairs in between guide hairs
vector<RenderedHair*> renderedHairs;

DemoBase *base;
const unsigned int nParticles = 50;
const Real helixRadius = 0.5;
const Real helixHeight = -5.0;
const Real helixTotalAngle = static_cast<Real>(2.0*M_PI);
const Matrix3r helixOrientation = AngleAxisr(-static_cast<Real>(-0.50 * M_PI), Vector3r(1,0,0)).toRotationMatrix();
bool drawFrames = false;


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS

	base = new DemoBase();
	base->init(argc, argv, "Cosserat rods demo");

	SimulationModel *model = new SimulationModel();
	model->init();
	Simulation::getCurrent()->setModel(model);

	buildModel();

	base->createParameterGUI();

	// add additional parameter just for this demo
	imguiParameters::imguiBoolParameter* param = new imguiParameters::imguiBoolParameter();
	param->description = "Draw frames";
	param->label = "Draw frames";
	param->readOnly = false;
	param->getFct = []() -> bool { return drawFrames; };
	param->setFct = [](bool b) -> void { drawFrames = b; };
	imguiParameters::addParam("Visualization", "Elastic rods", param);


	// OpenGL
	MiniGL::setClientIdleFunc (timeStep);		
	MiniGL::addKeyFunc('r', reset);
	MiniGL::setClientSceneFunc(render);			
	MiniGL::setViewport (40.0f, 0.1f, 500.0f, Vector3r (5.0, 10.0, 30.0), Vector3r (5.0, 0.0, 0.0));
	MiniGL::mainLoop ();	

	Utilities::Timing::printAverageTimes();
	Utilities::Timing::printTimeSums();

	delete Simulation::getCurrent();
	delete base;
	delete model;

	return 0;
}


void reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Simulation::getCurrent()->reset();
	base->getSelectedParticles().clear();

	Simulation::getCurrent()->getModel()->cleanup();
	buildModel();
}


void timeStep ()
{
	const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
	if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
		base->setValue(DemoBase::PAUSE, true);

	if (base->getValue<bool>(DemoBase::PAUSE))
		return;

	// Simulation code
	SimulationModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
	for (unsigned int i = 0; i < numSteps; i++)
	{
		START_TIMING("SimStep");
		Simulation::getCurrent()->getTimeStep()->step(*model);
		STOP_TIMING_AVG;

		base->step();
	}

	// interpolate and reconstruct all rendered hairs
	for(RenderedHair* rh : renderedHairs) {
		rh->timeStep();
	}
}

void buildModel ()
{
	TimeManager::getCurrent ()->setTimeStepSize (static_cast<Real>(0.005));
	
	Vector3r ghRoot1(0, 0, 5);
	Vector3r ghRoot2(-2.5 * std::pow(3, 0.5), 0, -2.5);
	Vector3r ghRoot3(2.5 * std::pow(3, 0.5), 0, -2.5);
	// // IDENTICAL GUIDE HAIRS
	// createHelix(ghRoot1, helixOrientation, helixRadius, helixHeight, helixTotalAngle, nParticles);
	// createHelix(ghRoot2, helixOrientation, helixRadius, helixHeight, helixTotalAngle, nParticles);
	// createHelix(ghRoot3, helixOrientation, helixRadius, helixHeight, helixTotalAngle, nParticles);
	// NON-IDENTICAL GUIDE HAIRS
	createHelix(ghRoot1, helixOrientation, helixRadius, helixHeight, helixTotalAngle, nParticles);
	createHelix(ghRoot2, helixOrientation, helixRadius * 2, helixHeight, helixTotalAngle, nParticles);
	createHelix(ghRoot3, helixOrientation, helixRadius, helixHeight * 2, helixTotalAngle / 2, nParticles);

	
	// The root we specify ends up being the center of the helix,
	// but that's not where the first point lies.
	// The first point lies a distance (helixRadius) in the positive x-direction
	ghRoot1 += Vector3r(helixRadius, 0, 0);
	ghRoot2 += Vector3r(helixRadius * 2, 0, 0);
	ghRoot3 += Vector3r(helixRadius, 0, 0);

	renderedHairs = {};
	// put rendered hairs in exact same spots as guide hairs to see how they line up
	renderedHairs.push_back(new RenderedHair(ghRoot1));
	renderedHairs.push_back(new RenderedHair(ghRoot2));
	renderedHairs.push_back(new RenderedHair(ghRoot3));

	const Real rhStepSize = base->getValue<Real>(DemoBase::INTERPOLATED_HAIR_STEP_SIZE);

	// spawn in rendered hairs that are (rhStepSize) away from each other)
	Real slope1 = (ghRoot1.z() - ghRoot2.z()) / (ghRoot1.x() - ghRoot2.x());
	Real slope2 = (ghRoot3.z() - ghRoot1.z()) / (ghRoot3.x() - ghRoot1.x());
	const Real xStart = rhStepSize * std::ceil(ghRoot2.x() / rhStepSize);
	const Real xEnd = rhStepSize * std::floor(ghRoot3.x() / rhStepSize);
	for(Real x = xStart; x <= xEnd; x += rhStepSize) {
		Real zStart = ghRoot2.z();
		Real zEnd = x <= ghRoot1.x()
			? slope1 * (x - ghRoot2.x()) + ghRoot2.z()
			: slope2 * (x - ghRoot3.x()) + ghRoot3.z();
		for(Real z = zStart; z <= zEnd; z += rhStepSize) {
			renderedHairs.push_back(new RenderedHair(Vector3r(x, 0, z)));
		}
	}
}

void renderLineModels()
{
	// Draw simulation model
	SimulationModel *model = Simulation::getCurrent()->getModel();
	const ParticleData &pd = model->getParticles();
	const OrientationData &od = model->getOrientations();
	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	float green[4] = { 0.0f, 0.8f, 0.0f, 1 };
	float blue[4] = { 0.0f, 0.0f, 0.8f, 1 };

	// position-based interpolated points of each rendered hair
	// pIP[r] = list of vertex positions for rendered hair r
	vector<vector<Vector3r>> positionInterpolatedPoints(renderedHairs.size());
	for(int rhInd = 0; rhInd < positionInterpolatedPoints.size(); rhInd++) {
		positionInterpolatedPoints[rhInd] = vector<Vector3r>(renderedHairs[rhInd]->currState->positions.size());
		positionInterpolatedPoints[rhInd][0] = renderedHairs[rhInd]->root;
	}
	
	// store roots of all guide hairs
	vector<Vector3r> guideRoots(model->getLineModels().size());
	for (unsigned int ghInd = 0; ghInd < model->getLineModels().size(); ghInd++)
	{
		LineModel *lineModel = model->getLineModels()[ghInd];
		
		const unsigned int indexOffset = lineModel->getIndexOffset();
		const unsigned int rootIndex = lineModel->getEdges()[0].m_vert[0] + indexOffset;
		const Vector3r &root = pd.getPosition(rootIndex);
		guideRoots[ghInd] = root;
	}

	// barycentric weights of each guide hair for positioning rendered hairs
	// weights[r][g] = weight of guide hair g on rendered hair r
	vector<vector<Real>> guideWeights(renderedHairs.size());
	for(unsigned int i = 0; i < guideWeights.size(); i++) {
		guideWeights[i] = vector<Real>(model->getLineModels().size());
	}

	// calculate denominator in weights
	const Real lambdaDenom = (guideRoots[0] - guideRoots[2]).cross(guideRoots[1] - guideRoots[2]).norm();

	// calculate barycentric weights of each guide hair for each rendered hair
	for(unsigned int rhInd = 0; rhInd < guideWeights.size(); rhInd++) {
		for (unsigned int ghInd = 0; ghInd < model->getLineModels().size(); ghInd++)
		{
			LineModel *lineModel = model->getLineModels()[ghInd];
			const Vector3r prevRoot = guideRoots[(ghInd + 3 - 1) % 3];
			const Vector3r currRoot = guideRoots[(ghInd + 3 + 0) % 3];
			const Vector3r nextRoot = guideRoots[(ghInd + 3 + 1) % 3];
			const Real lambda = (renderedHairs[rhInd]->root - prevRoot).cross(nextRoot - prevRoot).norm() / lambdaDenom;
			guideWeights[rhInd][ghInd] = lambda;
		}
	}
	
	// render all guide hairs while also interpolating each
	// rendered hair (if user wants to do position-based interpolation)
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
			
			MiniGL::drawSphere(v1, 0.07f, blue);
			if( e == lineModel->getEdges().size() -1 ) MiniGL::drawSphere(v2, 0.07f, blue);
			if(drawFrames) MiniGL::drawCylinder(v1, v2, blue, 0.01f);
			else MiniGL::drawCylinder(v1, v2, blue, 0.07f);
			
			// position-based interpolation
			for(unsigned int rhInd = 0; rhInd < positionInterpolatedPoints.size(); rhInd++) {
				if(i == 0) {
					positionInterpolatedPoints[rhInd][e+1] = Vector3r(0, 0, 0);
				}
				positionInterpolatedPoints[rhInd][e+1] += guideWeights[rhInd][i] * v2;
			}

			//draw coordinate frame at the center of the edges
			if(drawFrames)
			{
				Vector3r vm = 0.5 * (v1 + v2);
				Real scale = static_cast<Real>(0.15);
				Vector3r d1 = q._transformVector(Vector3r(1, 0, 0)) * scale;
				Vector3r d2 = q._transformVector(Vector3r(0, 1, 0)) * scale;
				Vector3r d3 = q._transformVector(Vector3r(0, 0, 1)) * scale;
				MiniGL::drawCylinder(vm, vm + d1, red,   0.01f);
				MiniGL::drawCylinder(vm, vm + d2, green, 0.01f);
				MiniGL::drawCylinder(vm, vm + d3, blue,  0.01f);
			}
		}
	}

	if(base->getValue<bool>(DemoBase::INTERPOLATION_METHOD)) {
		// FORCE-BASED INTERPOLATION DRAWING
		for(RenderedHair* rh : renderedHairs) {
			for(unsigned int i = 0; i < rh->currState->positions.size() - 1; ++i) {
				const Vector3r v1 = rh->currState->positions[i];
				const Vector3r v2 = rh->currState->positions[i + 1];
				MiniGL::drawCylinder(v1, v2, blue, 0.01f);
			}
		}
	} else {
		// POSITION-BASED INTERPOLATION DRAWING
		for(unsigned int rhInd = 0; rhInd < positionInterpolatedPoints.size(); rhInd++) {
			for(unsigned int vInd = 0; vInd < positionInterpolatedPoints[rhInd].size() - 1; vInd++) {
				const Vector3r v1 = positionInterpolatedPoints[rhInd][vInd];
				const Vector3r v2 = positionInterpolatedPoints[rhInd][vInd + 1];
				MiniGL::drawCylinder(v1, v2, blue, 0.01f);
			}
		}
	}

}

void render ()
{
	base->render();

	renderLineModels();
}

/** Create a particle model and orientation model
*/
void createHelix(const Vector3r &position, const Matrix3r &orientation, Real radius, Real height, Real totalAngle, int nPoints)
{
	int nQuaternions = nPoints - 1;
	vector<Vector3r> points(nPoints);
	vector<Quaternionr> quaternions(nQuaternions);
	
	//init particles
	for (int i = 0; i<nPoints; i++)   
	{
		points[i].x() = radius * std::cos(totalAngle / ((Real)nPoints) * (Real)i);
		points[i].y() = radius * std::sin(totalAngle / ((Real)nPoints) * (Real)i);
		points[i].z() = height / ((Real)nPoints) * (Real)i;

		points[i] = orientation * points[i] + position;
	}

	//init quaternions
	Vector3r from(0, 0, 1);
	for(int i=0; i<nQuaternions; i++)	
	{
		Vector3r to = (points[i + 1] - points[i]).normalized();
		Quaternionr dq = Quaternionr::FromTwoVectors(from, to);
		if(i == 0) quaternions[i] = dq;
		else quaternions[i] = dq * quaternions[i - 1];
		from = to;
	}

	vector<unsigned int> indices(2 * nPoints - 1);
	vector<unsigned int> indicesQuaternions(nQuaternions);

	for(int i=0; i < nPoints -1; i++)
	{
		indices[2 * i] = i;
		indices[2 * i + 1] = i + 1;
	}

	for (int i = 0; i < nQuaternions; i++)
	{
		indicesQuaternions[i] = i;
	}

	SimulationModel *model = Simulation::getCurrent()->getModel();
	model->addLineModel(nPoints, nQuaternions, &points[0], &quaternions[0], &indices[0], &indicesQuaternions[0]);

	ParticleData &pd = model->getParticles();
	const int nPointsTotal = pd.getNumberOfParticles();
	for (int i = nPointsTotal - 1; i > nPointsTotal - nPoints; i--)
	{
		pd.setMass(i, 1.0);
	}

	// Set mass of points to zero => make it static
	pd.setMass(nPointsTotal - nPoints, 0.0);

	OrientationData &od = model->getOrientations();
	const unsigned int nQuaternionsTotal = od.getNumberOfQuaternions();
	for(unsigned int i = nQuaternionsTotal - 1; i > nQuaternionsTotal - nQuaternions; i--)
	{
		od.setMass(i, 1.0);
	}
	
	// Set mass of quaternions to zero => make it static
	od.setMass(nQuaternionsTotal - nQuaternions, 0.0);

	// init constraints
	const size_t rodNumber = model->getLineModels().size() - 1;
	const unsigned int offset = model->getLineModels()[rodNumber]->getIndexOffset();
	const unsigned int offsetQuaternions = model->getLineModels()[rodNumber]->getIndexOffsetQuaternions();
	const size_t nEdges = model->getLineModels()[rodNumber]->getEdges().size();
	const LineModel::Edges &edges = model->getLineModels()[rodNumber]->getEdges();
		
	//stretchShear constraints
	for(unsigned int i=0; i < nEdges; i++)
	{
		const unsigned int v1 = edges[i].m_vert[0] + offset;
		const unsigned int v2 = edges[i].m_vert[1] + offset;
		const unsigned int q1 = edges[i].m_quat + offsetQuaternions;
		model->addStretchShearConstraint(v1, v2, q1, model->getRodStretchingStiffness(), model->getRodShearingStiffnessX(), model->getRodShearingStiffnessY());
	}

	//bendTwist constraints
	for(unsigned int i=0; i < nEdges - 1; i++)
	{
		const unsigned int q1 = edges[i].m_quat + offsetQuaternions;
		const unsigned int q2 = edges[i + 1].m_quat + offsetQuaternions;
		model->addBendTwistConstraint(q1, q2, model->getRodTwistingStiffness(), model->getRodBendingStiffnessX(), model->getRodBendingStiffnessY());
	}
	
// 	LOG_INFO << "Number of particles: " << nPoints;
// 	LOG_INFO << "Number of quaternions: " << nQuaternions;
}

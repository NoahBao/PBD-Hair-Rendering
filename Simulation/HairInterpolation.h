#ifndef _HairInterpolation_H
#define _HairInterpolation_H

#include "Common/Common.h"
#include "Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "Simulation/SimulationModel.h"
#include "Simulation/Simulation.h"
#include <iostream>

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;

namespace HairInterpolation {

    class RenderedHair {
        public:
            struct HairState {
                // list of positions for each vertex
                vector<Vector3r> positions;
                // list of orientations for each edge
                vector<Quaternionr> orientations;
            };

            // root of rendered hair
            const Vector3r root;
            // original positions and orientations of rendered hair vertices and edges
            HairState *restState;
            // current positions and orientations of rendered hair vertices and edges
            HairState *currState;

            RenderedHair(const Vector3r r) : root(r), restState(new HairState()), currState(new HairState()) {
                init();
            }

            void timeStep();
            void init();
            Quaternionr computeDarbouxVector(int j);
            int computePhiScalar(int j);
            Quaternionr computeBQuat(int i);
            Real computeLambdaScalar(int i, Vector3r netForce);
            Matrix4d computeInvertedMatrix(int i, Vector3r netForce);
            Quaternionr computeNewOrientation(int i, Vector3r netForce);
            Vector3r computeSegmentTangentVector(Quaternionr q);
            void updateRenderedHair(vector<Vector3r> netForces);
    };
    Real getStretchShearStiffness();
    Real getBendTwistStiffness();
    Quaternionr addQuats(Quaternionr q1, Quaternionr q2);
    Quaternionr scaleQuat(Quaternionr q, Real s);
    std::string vecToString(Vector3r v);
}

#endif
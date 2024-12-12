// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora
// =============================================================================
//
// Demo code about
// - using the ChParticleEmitter to create a cluster of random shapes
// - applying custom force field to particles
// - using Irrlicht to display objects
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/assets/ChTexture.h"
#include <fstream> // Include for file output
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

using namespace chrono;
using namespace chrono::particlefactory;
using namespace chrono::irrlicht;

//     A callback executed at each particle creation can be attached to the emitter.
//     For example, we need that new particles will be bound to Irrlicht visualization:
class MyCreatorForAll : public ChRandomShapeCreator::AddBodyCallback {
public:
    virtual void OnAddBody(std::shared_ptr<ChBody> mbody,
        ChCoordsys<> mcoords,
        ChRandomShapeCreator& mcreator) override {
        // Bind visual model to the visual system
        mbody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
        vis->BindItem(mbody);

        // Bind the collision model to the collision system
        if (mbody->GetCollisionModel())
            coll->Add(mbody->GetCollisionModel());

        // Dsable gyroscopic forces for increased integrator stability
        mbody->SetUseGyroTorque(false);
    }
    ChVisualSystem* vis;
    ChCollisionSystem* coll;
};

int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    // Create a Chrono physical system
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Simulation parameters
    int velcase = 2;
    double initialSpin = 0.0003;
    double maxT = 10.0;
    double time_step = 0.5;
    double spinupT = 0;
    double deltaSpin = 0;
    int SaveResults = 3;
    int SaveT = 240;
    int SaveVideo = 240;
    double density = 2000;
    double restituion = 0.0;
    double adhesion = 0.0;
    double adhesionMult = 0.0;
    double Kn = 200000;
    double Kt = 200000;
    double Gn = 40;
    double Gt = 20;
    // Aggregate characteristics
    double characteristic_lengths[3] = { 102,100,109 };
    double geometrical_slenderness_ratio = 0.921;
    double total_volume = 6.090e+05;
    double porosity = 0.462;
    double bulk_density = 1200;


    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Particle emitter");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox(GetChronoDataFile("skybox/sky_up.jpg"));
    vis->AddTypicalLights();
    vis->AddCamera(ChVector3d(0, 100, -100));

    // Create a rigid body
    auto sphere_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    sphere_mat->SetFriction(0.6f);
    auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    particle_mat->SetFriction(0.6f);
   // sphere_mat->SetAdhesion(adhesion);
    //sphere_mat->SetRestitution(restituion);
    //sphere_mat->SetAdhesionMultDMT(adhesionMult);
    //sphere_mat->SetKn(Kn);
   // sphere_mat->SetKt(Kt);
   // sphere_mat->SetGn(Gn);
   // sphere_mat->SetGt(Gt);
    /*
    auto sphereBody = chrono_types::make_shared<ChBodyEasySphere>(10,          // radius size
        density,         // density
        true,         // visualization?
        true,         // collision?
        sphere_mat);  // contact material
    sphereBody->SetPos(ChVector3d(1, 1, 0));
    sphereBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    sys.Add(sphereBody);
    */
    std::vector<chrono::ChVector3d> vertices;
    vertices.push_back(ChVector3d(10, 10, 10)); // Add more vertices to form a shape
    
    // Create the convex hull body using the vector of vertices
    auto sphereBody = chrono_types::make_shared<chrono::ChBodyEasyConvexHull>(
        vertices,   // vector of vertices
        density,        // density
        true,       // visualization?
        true,       // collision?
        sphere_mat  // contact material
    );

    sphereBody->SetPos(ChVector3d(1, 1, 0));
    sphereBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    sys.Add(sphereBody);
 
    // Create an emitter:
    ChParticleEmitter emitter;

    emitter.ParticlesPerSecond() = 20000;

    emitter.SetUseParticleReservoir(true);
    emitter.ParticleReservoirAmount() = 200;

    // Our ChParticleEmitter object, among the main settings, it requires
    // that you give him four 'randomizer' objects: one is in charge of
    // generating random shapes, one is in charge of generating
    // random positions, one for random alignements, and one for random velocities.
    // In the following we need to instance such objects. (There are many ready-to-use
    // randomizer objects already available in chrono, but note that you could also
    // inherit your own class from these randomizers if the choice is not enough).

    // ---Initialize the randomizer for POSITIONS: random points in a large cube
    auto emitter_positions = chrono_types::make_shared<ChRandomParticlePositionOnGeometry>();
    emitter_positions->SetGeometry(chrono_types::make_shared<ChBox>(60,60,60), ChFrame<>());
    emitter.SetParticlePositioner(emitter_positions);

    // ---Initialize the randomizer for ALIGNMENTS
    auto emitter_rotations = chrono_types::make_shared<ChRandomParticleAlignmentUniform>();
    emitter.SetParticleAligner(emitter_rotations);

    // ---Initialize the randomizer for VELOCITIES, with statistical distribution
    auto mvelo = chrono_types::make_shared<ChRandomParticleVelocityAnyDirection>();
    mvelo->SetModulusDistribution(chrono_types::make_shared<ChUniformDistribution>(0.0, 0.1));
    emitter.SetParticleVelocity(mvelo);

    // ---Initialize the randomizer for ANGULAR VELOCITIES, with statistical distribution
    auto mangvelo = chrono_types::make_shared<ChRandomParticleVelocityAnyDirection>();
    mangvelo->SetModulusDistribution(chrono_types::make_shared<ChUniformDistribution>(initialSpin, initialSpin));
    emitter.SetParticleAngularVelocity(mangvelo);

    // ---Initialize the randomizer for CREATED SHAPES, with statistical distribution

    /*
    // Create a ChRandomShapeCreator object (ex. here for sphere particles)
    auto mcreator_spheres = chrono_types::make_shared<ChRandomShapeCreatorSpheres>();
    mcreator_spheres->SetDiameterDistribution(chrono_types::make_shared<ChUniformDistribution>(2.1, 2.1));
    mcreator_spheres->SetDensityDistribution(chrono_types::make_shared<ChConstantDistribution>(1800));
    emitter.SetParticleCreator(mcreator_spheres);
    */

    // ..as an alternative: create odd shapes with convex hulls, like faceted fragments:
    auto mcreator_hulls = chrono_types::make_shared<ChRandomShapeCreatorConvexHulls>();
    mcreator_hulls->SetNpoints(16);
    mcreator_hulls->SetChordDistribution(chrono_types::make_shared<ChZhangDistribution>(10, 10));
    mcreator_hulls->SetDensityDistribution(chrono_types::make_shared<ChConstantDistribution>(density));

    emitter.SetParticleCreator(mcreator_hulls);
    
    // --- Optional: what to do by default on ALL newly created particles?
    //     A callback executed at each particle creation can be attached to the emitter.
    //     For example, we need that new particles will be bound to Irrlicht visualization:

    // a- define a class that implement your custom OnAddBody method (see top of source file)
    // b- create the callback object...
    auto mcreation_callback = chrono_types::make_shared<MyCreatorForAll>();
    // c- set callback own data that he might need...
    mcreation_callback->vis = vis.get();
    mcreation_callback->coll = sys.GetCollisionSystem().get();
    // d- attach the callback to the emitter!
    emitter.RegisterAddBodyCallback(mcreation_callback);

    // Bind all existing visual shapes to the visualization system
    vis->AttachSystem(&sys);

    // Modify some setting of the physical system for the simulation, if you want
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(40);

    // Turn off default -9.8 downward gravity
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));

    // Write the documents

    std::ofstream gravity("Grav Force");

    gravity << "Gravity Forces for each body  at each time step.\n";

    std::ofstream barycenterD("Barycenter Distances");

    barycenterD << "Barycenter distance in Caardinate each body each time stepn" << std::endl;

    std::ofstream contactforce("Contact Forces");

    contactforce << "Resultant Contact Forces in Cartesian each body each Time step"<< std::endl;

    std::ofstream contactpoint("Contact Points");

    contactpoint << "Contact Point for each pair in cartesian in Time Step" << std::endl;


    // Simulation loop
    double timestep = 0.01;
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        // Create particle flow
        emitter.EmitParticles(sys, timestep);

        // Apply custom forcefield (brute force approach..)
        // A) reset 'user forces accumulators':
        for (auto body : sys.GetBodies()) {
            body->EmptyAccumulators();
        }

        ChVector3<> barycenter(0, 0, 0);
        unsigned int body_index = 0;
        for (auto body : sys.GetBodies()) {
            ChVector3<> distance_to_barycenter = (body->GetPos() - barycenter);
            barycenterD << sys.GetChTime() << "\t" << body_index << "\t" << distance_to_barycenter.x() << "\t" << distance_to_barycenter.y() << "\t" << distance_to_barycenter.z() << "\t" << std::endl;
            ChVector3<> contact_force = body->GetContactForce();
            contactforce << sys.GetChTime() << "\t" << body_index << "\t" << contact_force.x() << "\t" << contact_force.y() << "\t" << contact_force.z() << std::endl;
            body_index++;
        }

        // B) store user computed force:
        // double G_constant = 6.674e-11; // gravitational constant
        double G_constant = 6.674e-9;  // gravitational constant - HACK to speed up simulation
        for (unsigned int i = 0; i < sys.GetBodies().size(); i++) {
            auto abodyA = sys.GetBodies()[i];
            for (unsigned int j = i + 1; j < sys.GetBodies().size(); j++) {
                auto abodyB = sys.GetBodies()[j];
                ChVector3d D_attract = abodyB->GetPos() - abodyA->GetPos();
                double r_attract = D_attract.Length();
                double f_attract = G_constant * (abodyA->GetMass() * abodyB->GetMass()) / (std::pow(r_attract, 2));
                ChVector3d F_attract = (D_attract / r_attract) * f_attract;

                abodyA->AccumulateForce(F_attract, abodyA->GetPos(), false);
                abodyB->AccumulateForce(-F_attract, abodyB->GetPos(), false);
            }
           // gravity << sys.GetChTime() << "\t" << i << "\t" << abodyA->GetAccumulatedForce().x() << "\t" << abodyA->GetAccumulatedForce().y() << "\t" << abodyA->GetAccumulatedForce().z() << std::endl;

        }

        // Perform the integration timestep
        sys.DoStepDynamics(timestep);
    }

    gravity.close();
    barycenterD.close();
    contactforce.close();
    contactpoint.close();

    return 0;
}

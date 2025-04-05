//////////////////////////////////// 
// Large Aggregate convex core
////////////////////////////////////
/*
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChContactContainer.h"
#include <fstream> // Include for file output
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "random"
#include "chrono/serialization//ChArchiveJSON.h"


using namespace chrono;
using namespace chrono::particlefactory;
using namespace chrono::irrlicht;
using namespace chrono::utils;

//     A callback executed at each particle creation can be attached to the emitter.
//     For example, we need that new particles will be bound to Irrlicht visualization:

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
    double density = 2800;
    double restitution = 0.0;
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

    auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    particle_mat->SetFriction(0.6f);
    particle_mat->SetRestitution(restitution);

    std::srand(std::time(nullptr)); // Seed the random number generator
    std::vector<chrono::ChVector3d> vertices;
    std::ofstream verticesfile("Vertices.txt");

    const unsigned int Nbodies = 7999;
    int npoints = 14;
    double mchord = 5;
    double msizeratioYZ = 1.0;
    double msizeratioZ = 1.0;
    double x; double y; double z;

    for (size_t i = 0; i < Nbodies; i++) {
        if (i == 0)
        {
            mchord = 20;
        }
        else
        {
            mchord = 5;
        }

        //	look why they are not added to the system just created
        std::vector<ChVector3d> points;
        points.resize(npoints);
        double hsizex = 0;
        double hsizey = 0;
        double hsizez = 0;
        for (int ip = 0; ip < npoints; ++ip) {
            points[ip].x() = ChRandom::Get();
            points[ip].y() = ChRandom::Get();
            points[ip].z() = ChRandom::Get();
            if (fabs(points[ip].x()) > hsizex)
                hsizex = fabs(points[ip].x());
            if (fabs(points[ip].y()) > hsizey)
                hsizey = fabs(points[ip].y());
            if (fabs(points[ip].z()) > hsizez)
                hsizez = fabs(points[ip].z());
        }
        for (int ip = 0; ip < npoints; ++ip) {
            points[ip].x() *= 0.5 * mchord / hsizex;
            points[ip].y() *= msizeratioYZ * (0.5 * mchord / hsizey);
            points[ip].z() *= msizeratioYZ * (0.5 * mchord / hsizez) * msizeratioZ;
            verticesfile << "[" << points[ip].x() << "," << points[ip].y() << "," << points[ip].z() << "]" << "\t";
        }
        verticesfile << std::endl;

        auto mbody = chrono_types::make_shared<ChBodyEasyConvexHull>(points, density, true, true, particle_mat);
        mbody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));

        // Create the placement box

        if (i == 0)
        {
            x = 0;
            y = 0;
            z = 0;
        }
        else
        {
            x = 200 * (-1.0 + 2.0 * (rand() / (double)RAND_MAX));
            y = 200 * (-1.0 + 2.0 * (rand() / (double)RAND_MAX));
            z = 200 * (-1.0 + 2.0 * (rand() / (double)RAND_MAX));
        }

 
        mbody->SetPos(ChVector3d(x, y, z));
        mbody->EnableCollision(true);
        sys.Add(mbody);
    }
    verticesfile.close();

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Particle emitter with core convex");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox(GetChronoDataFile("skybox/sky_up.jpg"));
    vis->AddTypicalLights();
    vis->AddCamera(ChVector3d(20, 10, -10));
    vis->AttachSystem(&sys);

    // Modify some setting of the physical system for the simulation, if you want
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(500);

    // Turn off default -9.8 downward gravity
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));

    // Write the documents
    std::ofstream pos("Position_Var.txt");
    std::ofstream vel("Velocity_Var.txt");
    std::ofstream mass("Mass_Var.txt");
    std::ofstream quat("Orientation_Var.txt");

    // Simulation loop
    double timestep = 1;
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        // Create particle flow
    //	emitter.EmitParticles(sys, timestep);

        // Apply custom forcefield (brute force approach..)
        // A) reset 'user forces accumulators':
        unsigned int iter = 0;
        for (auto body : sys.GetBodies()) {
            body->EmptyAccumulators();
            iter++;
        }

        // B) store user computed force:g
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
        }

        // Perform the integration timestep
        sys.DoStepDynamics(timestep);

    }

    unsigned int iterafter = 0;
    for (auto body : sys.GetBodies()) {
        body->EmptyAccumulators();
        pos << iterafter << "\t" << body->GetPos().x() << "\t" << body->GetPos().y() << "\t" << body->GetPos().z() << std::endl;
        vel << iterafter << "\t" << body->GetPosDt().x() << "\t" << body->GetPosDt().y() << "\t" << body->GetPosDt().z() << std::endl;
        mass << iterafter << "\t" << body->GetMass() << std::endl;
        quat << iterafter << "\t" << body->GetRot().e0() << "\t" << body->GetRot().e1() << "\t" << body->GetRot().e2() << "\t" << body->GetRot().e3() << std::endl;
        iterafter++;
    }
    mass.close();
    pos.close();
    vel.close();
    quat.close();
    return 0;
}*/



#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChContactContainer.h"
#include <fstream> // Include for file 
#include <filesystem>
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/utils/ChUtilsCreators.h"
#include <random>


using namespace chrono;
using namespace chrono::particlefactory;
using namespace chrono::irrlicht;
using namespace chrono::utils;

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

        // Disable gyroscopic forces for increased integrator stability
        mbody->SetUseGyroTorque(false);
    }
    ChVisualSystem* vis;
    ChCollisionSystem* coll;
};
class ContactBodyReporter : public chrono::ChContactContainer::ReportContactCallback {
public:
    virtual bool OnReportContact(
        const ChVector3<>& pA,            // contact point on object A
        const ChVector3<>& pB,            // contact point on object B
        const ChMatrix33<>& plane_coord, // contact plane coords (normal, U, V)
        const double& distance,          // penetration distance
        const double& eff_radius,        // effective radius of curvature
        const ChVector3<>& react_forces,  // forces in contact plane
        const ChVector3<>& react_torques, // torques in contact plane
        chrono::ChContactable* objA,     // contactable object A
        chrono::ChContactable* objB      // contactable object B
    ) override {
        // Attempt to cast contactable objects to ChBody
        auto bodyA = dynamic_cast<chrono::ChBody*>(objA);
        auto bodyB = dynamic_cast<chrono::ChBody*>(objB);

        if (bodyA && bodyB) {
            std::cout << "Contact between Body A and Body B:" << std::endl;
            std::cout << " - Body A ID: " << bodyA->GetIdentifier() << std::endl;
            std::cout << " - Body B ID: " << bodyB->GetIdentifier() << std::endl;
        }
        else {
            std::cout << "Contact involves non-body objects." << std::endl;
        }
        return true; // Continue reporting contacts
    }
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
    double density = 2800;
    double restitution = 0.0;
    double adhesion = 0.6;
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
    vis->SetWindowTitle("Particle emitter convex fixed particles");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddTypicalLights();
    vis->AddCamera(ChVector3d(0, 25, -25));

    std::filesystem::create_directories("Vertices");
    auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    particle_mat->SetFriction(0.6f);
    particle_mat->SetRestitution(restitution);

    std::vector<ChVector3d> vertices;
    std::ofstream verticesCore("Vertices/VerticesCore.txt");
    for (int i = 0; i < 10; i++)
    {
        double x =  15 * (((double)rand() / (RAND_MAX)) - 1);
        double y = (15 * (((double)rand() / (RAND_MAX)) - 1));
        double z = (0.5 * 15 * (((double)rand() / (RAND_MAX)) - 1));
        vertices.push_back(ChVector3d(x, y, z ));
        std::cout << x << " " << y << " " << z << std::endl;
        verticesCore << x << "\t" << y << "\t" << z << std::endl;
    }
    verticesCore.close();
    
    auto sphere = chrono_types::make_shared<ChBodyEasyConvexHull>(vertices,
        density,
        true,
        true,
        particle_mat);
    sphere->SetPos(ChVector3d(0, 0, 0));
    sphere->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    sys.Add(sphere);

    // Create an emitter:
    ChParticleEmitter emitter;

    emitter.ParticlesPerSecond() = 20000;

    emitter.SetUseParticleReservoir(true);
    emitter.ParticleReservoirAmount() = 9999;


    // Our ChParticleEmitter object, among the main settings, it requires
    // that you give him four 'randomizer' objects: one is in charge of
    // generating random shapes, one is in charge of generating
    // random positions, one for random alignements, and one for random velocities.
    // In the following we need to instance such objects. (There are many ready-to-use
    // randomizer objects already available in chrono, but note that you could also
    // inherit your own class from these randomizers if the choice is not enough).

    // ---Initialize the randomizer for POSITIONS: random points in a large cube
    auto emitter_positions = chrono_types::make_shared<ChRandomParticlePositionOnGeometry>();
    //emitter_positions->SetGeometry(chrono_types::make_shared<ChSphere>(500000), ChFrame<>());
    emitter_positions->SetGeometry(chrono_types::make_shared<ChBox>(500, 500, 500), ChFrame<>());
    emitter.SetParticlePositioner(emitter_positions);

    // ---Initialize the randomizer for ALIGNMENTS
    auto emitter_rotations = chrono_types::make_shared<ChRandomParticleAlignmentUniform>();
    emitter.SetParticleAligner(emitter_rotations);

    // ---Initialize the randomizer for VELOCITIES, with statistical distribution
    auto mvelo = chrono_types::make_shared<ChRandomParticleVelocityAnyDirection>();
    mvelo->SetModulusDistribution(chrono_types::make_shared<ChUniformDistribution>(0.0, 0.0000));
    emitter.SetParticleVelocity(mvelo);

    // ---Initialize the randomizer for ANGULAR VELOCITIES, with statistical distribution
    auto mangvelo = chrono_types::make_shared<ChRandomParticleVelocityAnyDirection>();
    mangvelo->SetModulusDistribution(chrono_types::make_shared<ChUniformDistribution>(initialSpin, initialSpin));
    emitter.SetParticleAngularVelocity(mangvelo);

    //Idea: To speed up the process, create first a bulk of particles that reach equilibrium, over this add many other particles that are not in 
    //equilibrium and make the simulation like this. As the core is already in equilibrium , it will act as one and attract the others around it,
    //until all are in equilibrium. Check if it does make sens to have all the particles in contact in the initial part, or it should be some in contact
    //most of them not, as it is g = 0


        // ---Initialize the randomizer for CREATED SHAPES, with statistical distribution


        // Create a ChRandomShapeCreator object (ex. here for sphere particles)
    auto mcreator_spheres = chrono_types::make_shared<ChRandomShapeCreatorConvexHulls>();
    mcreator_spheres->SetChordDistribution(chrono_types::make_shared<ChConstantDistribution>(5));
    mcreator_spheres->SetNpoints(12);
    mcreator_spheres->SetDensityDistribution(chrono_types::make_shared<ChConstantDistribution>(density));
    mcreator_spheres->SetAddCollisionShape(true);
    emitter.SetParticleCreator(mcreator_spheres);

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
    sys.GetSolver()->AsIterative()->SetMaxIterations(500);

    // Turn off default -9.8 downward gravity
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));

    // Write the documents
    std::ofstream pos("Position_Fixed.txt");
    std::ofstream vel("Velocity_Fixed.txt");
    std::ofstream mass("Mass_Fixed.txt");
    std::ofstream quat("Orientation_Fixed.txt");

    // Simulation loop
    double timestep = 1;
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        // Create particle flow
        emitter.EmitParticles(sys, timestep);

        // Apply custom forcefield (brute force approach..)
        // A) reset 'user forces accumulators':
        unsigned int iter = 0;
        for (auto body : sys.GetBodies()) {
            body->EmptyAccumulators();
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
        }

        // Perform the integration timestep
        sys.DoStepDynamics(timestep);

    }
    // Take the values of the particles at the end when the body is aggregated
    //int volume;
    unsigned int iterafter = 0;
    for (auto body : sys.GetBodies()) {
        body->EmptyAccumulators();
        //  volume = density * body->GetMass();
        pos << iterafter << "\t" << body->GetPos() << std::endl;
        vel << iterafter << "\t" << body->GetPosDt() << std::endl;
        mass << iterafter << "\t" << body->GetMass() << std::endl;
        quat << iterafter << "\t" << body->GetRot().e0() << "\t" << body->GetRot().e1() << "\t" << body->GetRot().e2() << "\t" << body->GetRot().e3() << std::endl;
        iterafter++;
    }
    mass.close();
    pos.close();
    vel.close();
    quat.close();
    return 0;
}

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
// Authors: Josh Diyn
// =============================================================================
//
// Demonstration of using heighfield RigidTerrain patches around a vehicle.
// As an example, a perlin noise heightmap is generated as an array (could import
// from a file if desired..
// The Heightfield patch type uses a y=0 bottom approach, and can be used 
// with Unity/Unreal terrain types.
// Note: Visual mesh generation is the primary performance bottleneck so 
// RigidTerrain patches are limited to 256x256 resolution for the visual mesh
//
// =============================================================================

#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChInteractiveDriver.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_models/vehicle/hmmwv/HMMWV.h"

#include "chrono/physics/ChBodyEasy.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_vehicle/visualization/ChVehicleVisualSystemIrrlicht.h"

using namespace chrono::irrlicht;
using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::hmmwv;

// Main function
int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2025 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    double step_size = 1e-3;
    double tire_step_size = 1e-3;

    // Height map params
    double terrainWidth = 100;  // eg. scale across 100m x 100m
    double terrainHeight = 100;
    int heightMapNx = 1024,
        heightMapNy = 1024;  // note: irrlicht visual mesh for this patch is capped to 512 x 512 to prevent slowdown caused by rendering the mesh
    double heightAmp = 4.5;  // Scale the noise-based height map to this height (total)
        
    // Create vehicle
    HMMWV_Full hmmwv;
    hmmwv.SetContactMethod(ChContactMethod::NSC);
    hmmwv.SetChassisFixed(false);
    hmmwv.SetInitPosition(ChCoordsys<>(ChVector3d(3, 0, 2.75), QUNIT));
    hmmwv.SetEngineType(EngineModelType::SIMPLE);
    hmmwv.SetTransmissionType(TransmissionModelType::AUTOMATIC_SIMPLE_MAP);
    hmmwv.SetDriveType(DrivelineTypeWV::AWD);
    hmmwv.SetBrakeType(BrakeType::SHAFTS);
    hmmwv.SetTireType(TireModelType::TMEASY);
    hmmwv.SetTireStepSize(tire_step_size);
    hmmwv.SetChassisCollisionType(CollisionType::PRIMITIVES);
    
    hmmwv.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    hmmwv.Initialize();

    hmmwv.SetChassisVisualizationType(VisualizationType::MESH);
    hmmwv.SetSuspensionVisualizationType(VisualizationType::PRIMITIVES);
    hmmwv.SetSteeringVisualizationType(VisualizationType::PRIMITIVES);
    hmmwv.SetWheelVisualizationType(VisualizationType::MESH);
    hmmwv.SetTireVisualizationType(VisualizationType::MESH);

    auto sys = hmmwv.GetSystem();

    // Create the shared contact material for all patches
    auto patch_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    patch_mat->SetFriction(0.9f);
    patch_mat->SetRestitution(0.01f);

    // add a box for testing
    auto my_obstacle = chrono_types::make_shared<ChBodyEasyBox>(1, 0.5, 1, 1000, true, true, patch_mat);
    sys->Add(my_obstacle);
    my_obstacle->SetPos(ChVector3d(3, -3, 0.75));
    my_obstacle->SetRot(QuatFromAngleX(30));
    my_obstacle->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/cubetexture_wood.png"));

    // add a sphere for testing
    auto my_obstacle2 = chrono_types::make_shared<ChBodyEasySphere>(1, 100, true, true, patch_mat);
    sys->Add(my_obstacle2);
    my_obstacle2->SetPos(ChVector3d(3, -3, 3));
    my_obstacle2->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/spheretexture.png"));


    sys->SetNumThreads(3,4,8);
    //sys->SetSolverType(ChSolver::Type::BARZILAIBORWEIN); // slightly faster
    //auto solver = sys->GetSolver();
    //solver->AsIterative()->EnableWarmStart(true);
    //solver->AsIterative()->SetMaxIterations(100);  // 50 iterations is enough for this demo


    // Create the terrain
    RigidTerrain terrain(sys);
    
    // time the build - note the visual is what takes the longest. Try turning off the visual for a test of the collision shape speed
    auto start = std::chrono::high_resolution_clock::now();

    // -----------------------------------------------------------------------------
    // replace the old Perlin-heightfield call with the new j = 0 image overload
    // -----------------------------------------------------------------------------
    std::string heightmap_file = GetChronoDataFile("vehicle/terrain/height_maps/test64.bmp"); // other terrain works but terrain3 flat spots are causing issues
    double hMin = 0;     // black in image is mapped to this value
    double hMax = 4.5;   // map white up to hmax

    auto dynPatch = terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector3d(0, 0, -3), QuatFromAngleX(0)),  // patch at world origin
                                     heightmap_file,                                       // grayscale height-map
                                     terrainWidth, terrainHeight,                          // physical X/Y extents (m)
                                     hMin, hMax,                                           // height range (m)
                                     0.001f);                                                // build visual mesh


    dynPatch->SetColor(ChColor(0.5f, 0.7f, 0.6f));
    // Initialize terrain
    terrain.Initialize();
    // timer ouput
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Heightfield Terrain build time: " << std::chrono::duration<double, std::micro>(end - start).count() << " µs\n";

    // driver
    double render_step = 1.0 / 50;
    ChInteractiveDriver driver(hmmwv.GetVehicle());
    driver.SetSteeringDelta(render_step / 1.0);
    driver.SetThrottleDelta(render_step / 1.0);
    driver.SetBrakingDelta(render_step / 0.3);
    driver.Initialize();

    // Create visualisation sys - Irrlicht
    // TODO - shift to VSG for faster vis of high res terrain
    ChVisualSystem::Type vis_type = ChVisualSystem::Type::IRRLICHT;

    std::shared_ptr<ChVehicleVisualSystem> vis;
    auto v = chrono_types::make_shared<ChVehicleVisualSystemIrrlicht>();
    v->SetWindowTitle("Streaming Terrain Demo");
    v->SetChaseCamera(ChVector3d(0, 0, 0.75), 6.0, 0.75);
    v->Initialize();
    v->AddLightDirectional();
    v->AddSkyBox();
    v->AddLogo();
    v->AttachVehicle(&hmmwv.GetVehicle());
    v->AttachDriver(&driver);
    vis = v;
    hmmwv.GetVehicle().EnableRealtime(true);
    ///////////////////////////////////////////////////

    // Sim loop
    while (vis->Run()) {
        double time = sys->GetChTime();

        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        DriverInputs inputs = driver.GetInputs();
        driver.Synchronize(time);
        terrain.Synchronize(time);
        hmmwv.Synchronize(time, inputs, terrain);
        vis->Synchronize(time, inputs);

        driver.Advance(step_size);
        terrain.Advance(step_size);
        hmmwv.Advance(step_size);

        vis->Advance(step_size);
        
    }

    return 0;
}

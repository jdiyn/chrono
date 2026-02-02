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
#include <numeric>
#include <limits>

#include "chrono/input_output/ChUtilsInputOutput.h"
#include "chrono_vehicle/ChVehicleDataPath.h"
#include "chrono_vehicle/driver/ChInteractiveDriver.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono/solver/ChIterativeSolverVI.h"


#include "chrono/physics/ChBodyEasy.h"

#include "chrono/assets/ChVisualSystem.h"

#ifdef CHRONO_IRRLICHT
    #include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
    #include "chrono_vehicle/visualization/ChVehicleVisualSystemIrrlicht.h"
    using namespace chrono::irrlicht;
#endif
#ifdef CHRONO_VSG
    #include "chrono_vehicle/visualization/ChVehicleVisualSystemVSG.h"
   // #include "chrono_vsg/ChVisualSystemVSG.h"
using namespace chrono::vsg3d;
#endif

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
    int heightMapNx = 2056,
        heightMapNy = 2056;  // note: irrlicht visual mesh for this patch is capped to 512 x 512 to prevent slowdown caused by rendering the mesh
    double heightAmp = 4.5;  // Scale the noise-based height map to this height (total)
        
    // Create vehicle
    HMMWV_Reduced hmmwv;
    hmmwv.SetContactMethod(ChContactMethod::NSC);
    hmmwv.SetChassisFixed(false);
    hmmwv.SetInitPosition(ChCoordsys<>(ChVector3d(3, 0,1), QUNIT));
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
    patch_mat->SetFriction(0.5);
    //  patch_mat->SetRollingFriction(0.0001); // this causes problems currently with the heightfield patch response - i htink becasue NSC has issues with rolling friction - deeper chrono issue
    patch_mat->SetRestitution(0.01f);


    // add a box for testing
    auto my_obstacle = chrono_types::make_shared<ChBodyEasyBox>(1, 0.5, 1, 1000, true, true, patch_mat);
    sys->Add(my_obstacle);
    my_obstacle->SetPos(ChVector3d(5, -3.5, 0.75));
    my_obstacle->SetRot(QuatFromAngleX(30));
    my_obstacle->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/cubetexture_wood.png"));

    // add a sphere for testing
    auto sphereObject = chrono_types::make_shared<ChBodyEasySphere>(1, 1000, true, true, patch_mat);
    sys->Add(sphereObject);
    auto sphereStart = ChVector3d(6, -3, 4);
    sphereObject->SetPos(sphereStart);
    sphereObject->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/spheretexture.png"));

    // Set initial speed: rolling in horizontal direction
    double initial_angspeed = 5;
    double initial_linspeed = initial_angspeed;
    //sphereObject->SetAngVelParent(ChVector3d(0, initial_angspeed, 0));
    //sphereObject->SetPosDt(ChVector3d(initial_linspeed, 0, 0));


    sys->SetNumThreads(12,3,2);

    
    // Solver and integrator settings
  //  sys->SetNumThreads(std::max(6, ChOMP::GetNumProcs()-6), 6, 6);

    sys->SetSolverType(ChSolver::Type::PSOR);

    auto solver = std::static_pointer_cast<chrono::ChIterativeSolverVI>(sys->GetSolver());
    solver->SetMaxIterations(250);
    solver->SetOmega(0.8);
    solver->SetSharpnessLambda(1.0);
    solver->SetTolerance(1e-8);

    // Create the terrain
    RigidTerrain terrain(sys);
    
    // time the build - note the visual is what takes the longest. Try turning off the visual for a test of the collision shape speed
    auto start = std::chrono::high_resolution_clock::now();

    std::string heightmap_file = GetChronoDataFile("vehicle/terrain/height_maps/test64.bmp"); // other terrain works but terrain3 flat spots are causing issues
    double hMin = 0;     // black in image is mapped to this value
    double hMax = 4.5;   // map white up to hmax

    //auto dynPatch = terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector3d(4, -5, -1), QuatFromAngleX(0)),  // patch at world origin
    //                                 heightmap_file,                                       // grayscale height-map
    //                                 terrainWidth, terrainHeight,                          // physical X/Y extents (m)
    //                                 hMin, hMax,                                           // height range (m)
    //                                 0.001f);                                              // build visual mesh
    //

    //dynPatch->SetColor(ChColor(0.5f, 0.7f, 0.6f));

    //auto boxpatch = terrain.AddPatch(patch_mat, ChCoordsys<>(), 100, 100, 1, false);
    //boxpatch->SetColor(ChColor(0.5f, 0.7f, 0.6f));
    
    //// Flat heightfield (all heights = 0)
    //std::vector<double> flat_heights(heightMapNx * heightMapNy, 0.0);

    //auto flatPatch = terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector3d(0, 0, -3), QuatFromAngleX(0)),
    //                                  flat_heights,   // heights array
    //                                  heightMapNx,    // grid_nx
    //                                  heightMapNy,    // grid_ny
    //                                  terrainWidth,   // dimX
    //                                  terrainHeight,  // dimY
    //                                  0.001f,         // sweep_sphere_radius
    //                                  true            // visualize
    //);

    //flatPatch->SetColor(ChColor(0.5f, 0.7f, 0.6f));


    // Generate a �bowl� that is 5 m deep at the center, 0 m at the rim
    std::vector<double> bowl_heights(heightMapNx * heightMapNy);
    double halfX = terrainWidth / 2.0;
    double halfY = terrainHeight / 2.0;
    double R = std::min(halfX, halfY);  // horizontal radius (50 m)
    double dx = terrainWidth / (heightMapNx - 1);
    double dy = terrainHeight / (heightMapNy - 1);
    double desiredDepth = 5.0;  // maximum depth in metres

    for (int j = 0; j < heightMapNy; ++j) {
        double y = j * dy - halfY;
        for (int i = 0; i < heightMapNx; ++i) {
            double x = i * dx - halfX;
            double r2 = x * x + y * y;
            double h = 0.0;
            if (r2 <= R * R) {
                // spherical?cap style: -desiredDepth at center, 0 at rim
                h = -desiredDepth * std::sqrt(1.0 - r2 / (R * R));
            }
            bowl_heights[j * heightMapNx + i] = h;
        }
    }

    //auto bowlPatch = terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector3d(0, 0, 4.5), QUNIT), bowl_heights, heightMapNx,
    //                                  heightMapNy, terrainWidth, terrainHeight, 0.001f, true);
    //bowlPatch->SetColor(ChColor(0.5f, 0.7f, 0.6f));

    // =========================================================================
    // PERLIN NOISE HEIGHTFIELD PATCH
    // =========================================================================
    // Generate procedural terrain using Perlin noise for natural-looking hills
    // This demonstrates programmatic heightfield generation without external files
    
    // Perlin noise parameters
    int perlinNx = 2056;
    int perlinNy = 2056;
    double perlinWidth = 80.0;   // 80m x 80m patch
    double perlinLength = 80.0;
    double perlinAmplitude = 6.0;  // Max height variation in meters
    int perlinOctaves = 4;         // Number of noise layers
    double perlinPersistence = 0.5; // Amplitude falloff per octave
    double perlinFrequency = 0.02;  // Base frequency (lower = smoother)
    
    // Simple Perlin-like noise implementation using gradient vectors
    // Permutation table for pseudo-random gradient selection
    std::vector<int> perm(512);
    {
        std::vector<int> p(256);
        std::iota(p.begin(), p.end(), 0);
        std::mt19937 rng(42);  // Fixed seed for reproducibility
        std::shuffle(p.begin(), p.end(), rng);
        for (int i = 0; i < 256; ++i) {
            perm[i] = perm[i + 256] = p[i];
        }
    }
    
    // Gradient vectors (8 directions)
    const double gradients[8][2] = {
        {1, 0}, {-1, 0}, {0, 1}, {0, -1},
        {0.7071, 0.7071}, {-0.7071, 0.7071}, {0.7071, -0.7071}, {-0.7071, -0.7071}
    };
    
    // Fade function for smooth interpolation
    auto fade = [](double t) { return t * t * t * (t * (t * 6 - 15) + 10); };
    
    // Linear interpolation
    auto lerp = [](double a, double b, double t) { return a + t * (b - a); };
    
    // Gradient function
    auto grad = [&](int hash, double x, double y) {
        int h = hash & 7;
        return gradients[h][0] * x + gradients[h][1] * y;
    };
    
    // 2D Perlin noise function
    auto perlinNoise2D = [&](double x, double y) -> double {
        int X = static_cast<int>(std::floor(x)) & 255;
        int Y = static_cast<int>(std::floor(y)) & 255;
        x -= std::floor(x);
        y -= std::floor(y);
        double u = fade(x);
        double v = fade(y);
        int aa = perm[perm[X] + Y];
        int ab = perm[perm[X] + Y + 1];
        int ba = perm[perm[X + 1] + Y];
        int bb = perm[perm[X + 1] + Y + 1];
        return lerp(lerp(grad(aa, x, y), grad(ba, x - 1, y), u),
                    lerp(grad(ab, x, y - 1), grad(bb, x - 1, y - 1), u), v);
    };
    
    // Fractal Brownian Motion (fBm) - layered Perlin noise
    auto fbm = [&](double x, double y, int octaves, double persistence, double frequency) -> double {
        double total = 0.0;
        double amplitude = 1.0;
        double maxValue = 0.0;
        for (int i = 0; i < octaves; ++i) {
            total += perlinNoise2D(x * frequency, y * frequency) * amplitude;
            maxValue += amplitude;
            amplitude *= persistence;
            frequency *= 2.0;
        }
        return total / maxValue;  // Normalize to [-1, 1]
    };
    
    // Generate Perlin noise heightfield
    std::vector<double> perlin_heights(perlinNx * perlinNy);
    double perlinDx = perlinWidth / (perlinNx - 1);
    double perlinDy = perlinLength / (perlinNy - 1);
    
    // First pass: generate raw noise and track actual min/max
    double noiseMin = std::numeric_limits<double>::max();
    double noiseMax = std::numeric_limits<double>::lowest();
    
    for (int j = 0; j < perlinNy; ++j) {
        for (int i = 0; i < perlinNx; ++i) {
            double worldX = i * perlinDx;
            double worldY = j * perlinDy;
            double noise = fbm(worldX, worldY, perlinOctaves, perlinPersistence, perlinFrequency);
            perlin_heights[j * perlinNx + i] = noise;
            noiseMin = std::min(noiseMin, noise);
            noiseMax = std::max(noiseMax, noise);
        }
    }
    
    // Second pass: rescale to [0, perlinAmplitude] using actual range
    // This produces LOCAL heights (BASE-relative, min=0), matching Unreal/Unity terrain export format
    double noiseRange = noiseMax - noiseMin;
    if (noiseRange > 1e-9) {
        for (auto& h : perlin_heights) {
            h = ((h - noiseMin) / noiseRange) * perlinAmplitude;
        }
    }
    
    std::cout << "Perlin noise range: [" << noiseMin << ", " << noiseMax << "] -> scaled to [0, " << perlinAmplitude << "]" << std::endl;
    
    // Add the Perlin noise patch using LOCAL heights (heightsAreLocal = true)
    // This is how you'd import terrain from Unreal/Unity where heights are already BASE-relative
    auto perlinPatch = terrain.AddPatch(patch_mat, 
                                         ChCoordsys<>(ChVector3d(0, 0, -6), QUNIT),  // Position offset (BASE sits at z=-6)
                                         perlin_heights,
                                         perlinNx, perlinNy,
                                         perlinWidth, perlinLength,
                                         0.001f,   // sweep sphere radius
                                         true,     // build visual mesh
                                         true);    // heightsAreLocal=true: heights are already BASE-relative [0, amplitude]
    perlinPatch->SetColor(ChColor(0.6f, 0.55f, 0.4f));  // Sandy/dirt color
    
    // Example of ABSOLUTE heights (heightsAreLocal = false, default):
    // If you had real-world elevation data like GPS altitudes (e.g., 500m to 550m),
    // you'd pass heightsAreLocal=false and the system would automatically shift
    // minHeight to local z=0. Uncomment below to see:
    //
    // std::vector<double> absolute_heights(64 * 64);
    // for (int j = 0; j < 64; ++j) {
    //     for (int i = 0; i < 64; ++i) {
    //         // Simulate real elevation data (500m base + up to 50m variation)
    //         double noise = fbm(i * 0.1, j * 0.1, 4, 0.5, 0.5);
    //         absolute_heights[j * 64 + i] = 500.0 + (noise + 1.0) * 25.0;  // 500m to 550m
    //     }
    // }
    // auto absolutePatch = terrain.AddPatch(patch_mat,
    //                                       ChCoordsys<>(ChVector3d(50, 0, 0), QUNIT),  // BASE at z=0
    //                                       absolute_heights, 64, 64, 20.0, 20.0,
    //                                       0.001f, true, false);  // heightsAreLocal=false (default)



    // Initialize terrain
    terrain.Initialize();
    // timer ouput
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Heightfield Terrain build time: " << std::chrono::duration<double, std::micro>(end - start).count() << " �s\n";

    // driver
    double render_step = 1.0 / 50;
    ChInteractiveDriver driver(hmmwv.GetVehicle());
    driver.SetSteeringDelta(render_step / 1.0);
    driver.SetThrottleDelta(render_step / 1.0);
    driver.SetBrakingDelta(render_step / 0.3);
    driver.Initialize();

    // Create visualisation sys - Irrlicht
    // TODO - shift to VSG for faster vis of high res terrain


    std::shared_ptr<ChVehicleVisualSystem> vis;

#ifdef CHRONO_IRRLICHT
    ChVisualSystem::Type vis_type = ChVisualSystem::Type::IRRLICHT;
    auto v = chrono_types::make_shared<ChVehicleVisualSystemIrrlicht>();
    v->SetWindowTitle("Streaming Terrain Demo");
    v->SetChaseCamera(ChVector3d(0, 0, 0.75), 6.0, 0.75);
    v->Initialize();
    v->AddLightDirectional();
    v->AddSkyBox();
    v->AddLogo();
    v->AttachVehicle(&hmmwv.GetVehicle());
    v->AttachDriver(&driver);
#elif CHRONO_VSG
    ChVisualSystem::Type vis_type = ChVisualSystem::Type::VSG;
    auto v = chrono_types::make_shared<ChVehicleVisualSystemVSG>();
    v->SetWindowTitle("Streaming Terrain Demo");
    v->SetChaseCamera(ChVector3d(0, 0, 0.75), 6.0, 0.75);
    v->Initialize();
    v->AttachVehicle(&hmmwv.GetVehicle());
    v->AttachDriver(&driver);
#endif

    vis = v;
    hmmwv.GetVehicle().EnableRealtime(false);
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


    std::cout << "Sphere Starting point: "
              << sphereStart << "\nFinal resting place of sphere after rolling: " << sphereObject->GetPos() << std::endl;

    return 0;
}

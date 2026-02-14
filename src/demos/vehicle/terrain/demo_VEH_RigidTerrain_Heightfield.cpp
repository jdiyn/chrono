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
#include <iomanip>
#include <map>
#include <set>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <limits>

// ============================================================================
// STRESS TEST KNOBS (edit and recompile)
// ============================================================================
namespace {
// Max number of spheres created for the stress test.
constexpr int kMaxSpheres = 100;
}  // namespace

#include "chrono/input_output/ChUtilsInputOutput.h"
#include "chrono_vehicle/ChVehicleDataPath.h"
#include "chrono_vehicle/driver/ChInteractiveDriver.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono/solver/ChIterativeSolverVI.h"


#include "chrono/physics/ChBodyEasy.h"

#include "chrono/assets/ChVisualSystem.h"

//#ifdef CHRONO_IRRLICHT
//    #include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
//    #include "chrono_vehicle/visualization/ChVehicleVisualSystemIrrlicht.h"
//    using namespace chrono::irrlicht;
//#endif
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

    // Stress test configuration
    const int stress_spheres = std::max(1, kMaxSpheres);

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

    // Stress test spheres (collision-heavy). By default only the first sphere is dynamic/visualized.
    const double sphere_radius = 1.0;
    const double sphere_density = 1000.0;
    std::vector<std::shared_ptr<ChBody>> stress_sphere_bodies;
    stress_sphere_bodies.reserve(stress_spheres);

    auto sphereObject = chrono_types::make_shared<ChBodyEasySphere>(sphere_radius, sphere_density, true, true, patch_mat);
    sys->Add(sphereObject);
    sphereObject->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/spheretexture.png"));
    stress_sphere_bodies.push_back(sphereObject);

    for (int i = 1; i < stress_spheres; ++i) {
        auto s = chrono_types::make_shared<ChBodyEasySphere>(sphere_radius, sphere_density, true, true, patch_mat);
        sys->Add(s);
        s->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/spheretexture.png"));
        stress_sphere_bodies.push_back(s);
    }

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
    auto total_start = std::chrono::high_resolution_clock::now();

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
    auto perlin_gen_start = std::chrono::high_resolution_clock::now();
    
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
    
    auto perlin_gen_end = std::chrono::high_resolution_clock::now();
    double perlin_gen_ms = std::chrono::duration<double, std::milli>(perlin_gen_end - perlin_gen_start).count();
    std::cout << "Perlin noise generation: " << std::fixed << std::setprecision(2) << perlin_gen_ms << " ms\n";
    std::cout << "Perlin noise range: [" << noiseMin << ", " << noiseMax << "] -> scaled to [0, " << perlinAmplitude << "]\n";
    
    // Add the Perlin noise patch using LOCAL heights (heightsAreLocal = true)
    // This is how you'd import terrain from Unreal/Unity where heights are already BASE-relative
    auto patch_add_start = std::chrono::high_resolution_clock::now();
    
    auto perlinPatch = terrain.AddPatch(patch_mat, 
                                         ChCoordsys<>(ChVector3d(0, 0, -6), QUNIT),  // Position offset (BASE sits at z=-6)
                                         perlin_heights,
                                         perlinNx, perlinNy,
                                         perlinWidth, perlinLength,
                                         0.001f,   // sweep sphere radius
                                         true,     // build visual mesh
                                         true);    // heightsAreLocal=true: heights are already BASE-relative [0, amplitude]
    perlinPatch->SetColor(ChColor(0.6f, 0.55f, 0.4f));  // Sandy/dirt color
    
    auto patch_add_end = std::chrono::high_resolution_clock::now();
    double patch_add_ms = std::chrono::duration<double, std::milli>(patch_add_end - patch_add_start).count();
    std::cout << "AddPatch (collision + visual mesh): " << std::fixed << std::setprecision(2) << patch_add_ms << " ms\n";
    
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
    auto init_start = std::chrono::high_resolution_clock::now();
    terrain.Initialize();
    auto init_end = std::chrono::high_resolution_clock::now();
    double init_ms = std::chrono::duration<double, std::milli>(init_end - init_start).count();
    std::cout << "terrain.Initialize(): " << std::fixed << std::setprecision(2) << init_ms << " ms\n";

    // Place stress-test spheres AFTER terrain init (so we can query ground height)
    // - A subset is placed near the surface to generate contacts.
    // - The rest are placed high above the terrain to stress broadphase + early-out paths.
    const int total_spheres = static_cast<int>(stress_sphere_bodies.size());
    const bool disable_sphere_sphere = true;
    const int dynamic_spheres = std::min(total_spheres, 200);
    const int near_spheres = std::min(total_spheres, std::max(256, total_spheres / 10));
    std::mt19937 stress_rng(123);
    std::uniform_real_distribution<double> ux(-perlinWidth * 0.45, perlinWidth * 0.45);
    std::uniform_real_distribution<double> uy(-perlinLength * 0.45, perlinLength * 0.45);

    for (int i = 0; i < total_spheres; ++i) {
        auto& body = stress_sphere_bodies[i];
        if (disable_sphere_sphere) {
            body->GetCollisionModel()->SetFamily(1);
            body->GetCollisionModel()->SetFamily(1);
            body->GetCollisionModel()->CollidesWith(0);  // collide with terrain (family 0)
            body->GetCollisionModel()->DisallowCollisionsWith(1);  // no sphere-sphere collisions

        }

        body->SetFixed(i >= dynamic_spheres);

        // Keep the original demo sphere for consistent behavior/comparison.
        if (i == 0) {
            body->SetPos(ChVector3d(6, -3, 4));
        } else {
            const double x = ux(stress_rng);
            const double y = uy(stress_rng);
            const double base_ground = terrain.GetHeight(ChVector3d(x, y, 0));
            const double z = (i < near_spheres)
                                 ? (base_ground + sphere_radius + 0.05 + 0.01 * (i % 10))
                                 : (50.0 + 0.02 * (i - near_spheres));
            body->SetPos(ChVector3d(x, y, z));
        }
        body->SetPosDt(ChVector3d(0, 0, 0));
        body->SetAngVelParent(ChVector3d(0, 0, 0));
    }

    // Track the first sphere for reporting
    auto sphereStart = sphereObject->GetPos();
    std::cout << "Stress spheres: " << total_spheres
              << " (dynamic=" << dynamic_spheres
              << ", near_surface=" << near_spheres
              << ", sphere_sphere_collisions=OFF"
              << ")\n";
    
    // timer ouput
    auto total_end = std::chrono::high_resolution_clock::now();
    auto build_time_ms = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Print terrain build statistics
    std::cout << "\n========== HEIGHTFIELD TERRAIN STATISTICS ==========\n";
    std::cout << "Build time: " << std::fixed << std::setprecision(2) << build_time_ms << " ms\n";
    std::cout << "Grid resolution: " << perlinNx << " x " << perlinNy << " = " 
              << (perlinNx * perlinNy) << " vertices\n";
    std::cout << "Physical size: " << perlinWidth << " x " << perlinLength << " m\n";
    std::cout << "Cell size: " << std::setprecision(4) << (perlinWidth / (perlinNx - 1)) << " m\n";
    
    // Memory estimate for collision shape (approximate)
    size_t heightDataMB = (perlinNx * perlinNy * sizeof(double)) / (1024 * 1024);
    size_t vertexCacheMB = 0;  // Only for small terrains
    if (perlinNx * perlinNy <= 512 * 512) {
        vertexCacheMB = (perlinNx * perlinNy * 3 * sizeof(double)) / (1024 * 1024);
    }
    // Quad extents cache is always enabled (default m_useQuadExtentsCache = true).
    // Each quad stores min/max height (2 scalars).
    size_t quadExtentsMB = ((perlinNx - 1) * (perlinNy - 1) * 2 * sizeof(double)) / (1024 * 1024);
    std::cout << "Estimated collision memory: " << (heightDataMB + vertexCacheMB + quadExtentsMB) << " MB\n";
    std::cout << "  - Height data: " << heightDataMB << " MB\n";
    std::cout << "  - Vertex cache: " << (vertexCacheMB > 0 ? std::to_string(vertexCacheMB) + " MB" : "disabled (large terrain)") << "\n";
    std::cout << "  - Quad extents: " << quadExtentsMB << " MB\n";
    std::cout << "====================================================\n\n";

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
//
//#ifdef CHRONO_IRRLICHT
//    ChVisualSystem::Type vis_type = ChVisualSystem::Type::IRRLICHT;
//    auto v = chrono_types::make_shared<ChVehicleVisualSystemIrrlicht>();
//    v->SetWindowTitle("Streaming Terrain Demo");
//    v->SetChaseCamera(ChVector3d(0, 0, 0.75), 6.0, 0.75);
//    v->Initialize();
//    v->AddLightDirectional();
//    v->AddSkyBox();
//    v->AddLogo();
//    v->AttachVehicle(&hmmwv.GetVehicle());
//    v->AttachDriver(&driver);
//#elif CHRONO_VSG
#ifdef CHRONO_VSG   
    ChVisualSystem::Type vis_type = ChVisualSystem::Type::VSG;
    auto v = chrono_types::make_shared<ChVehicleVisualSystemVSG>();
    v->SetWindowTitle("Streaming Terrain Demo");
    v->SetChaseCamera(ChVector3d(0, 0, 15.75), 16.0, 0.75);
    v->SetTargetRenderFPS(60); // keep the rendering at 60, unconnected to the step
    v->AttachVehicle(&hmmwv.GetVehicle());
    v->AttachDriver(&driver);
    v->Initialize();
#endif

    vis = v;
    hmmwv.GetVehicle().EnableRealtime(false);
    ///////////////////////////////////////////////////

    // Performance tracking variables
    int frame_count = 0;
    int report_interval = 100;  // Report every 100 frames
    double total_step_time = 0.0;
    double max_step_time = 0.0;
    double min_step_time = 1e9;
    auto last_report_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "\n========== STARTING SIMULATION ==========\n";
    std::cout << "Performance reports every " << report_interval << " frames\n";
    std::cout << "==========================================\n\n";

    // Sim loop
    int collision_count = 0;
while (vis->Run()) {
    auto step_start = std::chrono::high_resolution_clock::now();
    
    double time = sys->GetChTime();
    
    // Get contact count from container
    collision_count = sys->GetContactContainer()->GetNumContacts();
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
        
        //// Performance tracking
        //auto step_end = std::chrono::high_resolution_clock::now();
        //double step_time_ms = std::chrono::duration<double, std::milli>(step_end - step_start).count();
        //total_step_time += step_time_ms;
        //max_step_time = std::max(max_step_time, step_time_ms);
        //min_step_time = std::min(min_step_time, step_time_ms);
        //frame_count++;
        //
        //// Periodic performance report
        //if (frame_count % report_interval == 0) {
        //    auto now = std::chrono::high_resolution_clock::now();
        //    double elapsed_sec = std::chrono::duration<double>(now - last_report_time).count();
        //    double avg_step_ms = total_step_time / report_interval;
        //    double fps = report_interval / elapsed_sec;
        //    double realtime_factor = (step_size * report_interval) / elapsed_sec;
        //    
        //    std::cout << "Frame " << std::setw(6) << frame_count 
        //              << " | Sim time: " << std::fixed << std::setprecision(2) << std::setw(6) << time << "s"
        //              << " | Step: " << std::setprecision(2) << std::setw(5) << avg_step_ms << "ms"
        //              << " (min:" << std::setw(4) << min_step_time << ", max:" << std::setw(5) << max_step_time << ")"
        //              << " | FPS: " << std::setprecision(1) << std::setw(5) << fps
        //              << " | RTF: " << std::setprecision(2) << realtime_factor << "x"
        //          << " | Contacts: " << collision_count
        //          << " | Sphere Z: " << std::setprecision(3) << sphereObject->GetPos().z()
        //          << " | Ground: " << terrain.GetHeight(sphereObject->GetPos())
        //          << std::endl;
        //    
        //    // Reset for next interval
        //    total_step_time = 0.0;
        //    max_step_time = 0.0;
        //    min_step_time = 1e9;
        //    last_report_time = now;
        //}
    }


    std::cout << "Sphere Starting point: "
              << sphereStart << "\nFinal resting place of sphere after rolling: " << sphereObject->GetPos() << std::endl;

    return 0;
}

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
//    #include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtHeightfieldTerrainShape.h"


using namespace chrono::irrlicht;

using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::hmmwv;

// Simple Perlin noise implementation
class PerlinNoise {
  public:
    PerlinNoise(unsigned int seed = 12345) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < 256; ++i) {
            p[i] = p[i + 256] = i;
        }
        std::shuffle(p.begin(), p.begin() + 256, gen);
    }

    double noise(double x, double y) const {
        int X = (int)std::floor(x) & 255;
        int Y = (int)std::floor(y) & 255;
        x -= std::floor(x);
        y -= std::floor(y);
        double u = fade(x);
        double v = fade(y);

        int aa = p[p[X] + Y];
        int ab = p[p[X] + Y + 1];
        int ba = p[p[X + 1] + Y];
        int bb = p[p[X + 1] + Y + 1];

        double grad_aa = grad(p[aa], x, y);
        double grad_ba = grad(p[ba], x - 1, y);
        double grad_ab = grad(p[ab], x, y - 1);
        double grad_bb = grad(p[bb], x - 1, y - 1);

        return lerp(v, lerp(u, grad_aa, grad_ba), lerp(u, grad_ab, grad_bb));
    }

  private:
    std::array<int, 512> p;

    double fade(double t) const { return t * t * t * (t * (t * 6 - 15) + 10); }
    double lerp(double t, double a, double b) const { return a + t * (b - a); }
    double grad(int hash, double x, double y) const {
        int h = hash & 15;
        double u = h < 8 ? x : y;
        double v = h < 4 ? y : h == 12 || h == 14 ? x : 0;
        return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }
};

// Height map class to manage a cached nx×ny height array,
// stored row-major with j=0 at the bottom (y minimum) <-- this is how the heightfield class and Bullet also expect it
class HeightMap {
  public:
    HeightMap(double width, double height, int nx, int ny, double amplitude, unsigned int seed = 987654321)
        : width_(width), height_(height), nx_(nx), ny_(ny), amplitude_(amplitude) {
        // Resize storage: put j=0 is bottom row, j=ny-1 as top row
        heights_.resize(nx * ny);

        PerlinNoise noise(seed);
        double scale = 0.15;  // variation scale

        // fill heights
        for (int j_noise = 0; j_noise < ny; ++j_noise) {
            for (int i = 0; i < nx; ++i) {
                double x = i * width / (nx - 1);
                double y = j_noise * height / (ny - 1);
                // multi-octave perlin noise
                double n = 0.5 * noise.noise(x * scale, y * scale) + 0.25 * noise.noise(x * scale * 2, y * scale * 2) +
                           0.125 * noise.noise(x * scale * 4, y * scale * 4);
                int j_store = (ny - 1) - j_noise;
                heights_[j_store * nx + i] = amplitude * n;
            }
        }
    }

    // getter for full bottom-first array
    const std::vector<double>& getAllHeights() const { return heights_; }

    // Get height at sample (i,j), with j=0 at bottom.
    double getHeight(int i, int j) const {
        if (i < 0 || i >= nx_ || j < 0 || j >= ny_)
            return 0.0;
        return heights_[j * nx_ + i];
    }

  private:
    double width_, height_;
    int nx_, ny_;
    double amplitude_;
    std::vector<double> heights_;  ///< row-major, j=0 at the bottom
};


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

    // Create height map
    HeightMap heightMap(terrainWidth, terrainHeight, heightMapNx, heightMapNy, heightAmp);

        
    // Create vehicle
    HMMWV_Full hmmwv;
    hmmwv.SetContactMethod(ChContactMethod::NSC);
    hmmwv.SetChassisFixed(false);
    hmmwv.SetInitPosition(ChCoordsys<>(ChVector3d(3, 0, 2.75), QuatFromAngleX(65)));
    hmmwv.SetEngineType(EngineModelType::SIMPLE);
    hmmwv.SetTransmissionType(TransmissionModelType::AUTOMATIC_SIMPLE_MAP);
    hmmwv.SetDriveType(DrivelineTypeWV::AWD);
    hmmwv.SetBrakeType(BrakeType::SHAFTS);
    hmmwv.SetTireType(TireModelType::TMEASY);
    hmmwv.SetTireStepSize(tire_step_size);
    hmmwv.SetChassisCollisionType(CollisionType::MESH);
    
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
    my_obstacle->SetPos(ChVector3d(3, -3, 1));
    my_obstacle->SetRot(QuatFromAngleX(30));
    my_obstacle->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/cubetexture_wood.png"));

    // add a sphere for testing
    auto my_obstacle2 = chrono_types::make_shared<ChBodyEasySphere>(1, 100, true, true, patch_mat);
    sys->Add(my_obstacle2);
    my_obstacle2->SetPos(ChVector3d(6, 0, 5));
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

    // Use the entire height map for one large patch
    //auto heights = heightMap.getAllHeights();
    //auto dynPatch = terrain.AddPatch(patch_mat,                                 // NSC or SMC material
    //                                            ChCoordsys<>(ChVector3d(0, 0, 0), QUNIT),  // patch body at world origin, no rotation
    //                                            heightMap.getAllHeights(),                 // row major height array
    //                                            heightMapNx, heightMapNy,     // full grid resolution
    //                                            terrainWidth, terrainHeight,  // full terrain extents (Am x Bm)
    //                                            0.001f,
    //                                            true            // build chtrianglemesh visual shape
    //);

    // -----------------------------------------------------------------------------
    // replace the old Perlin-heightfield call with the new j = 0 image overload
    // -----------------------------------------------------------------------------
    std::string heightmap_file = GetChronoDataFile("vehicle/terrain/height_maps/terrain3.bmp"); // other terrain works but terrain3 flat spots are causing issues
    double hMin = -heightAmp * 0.5;  // map black is -amplitude
    double hMax = heightAmp * 0.5;   // map white up to +amplitude

    auto dynPatch = terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector3d(0, 0, 0), QUNIT),  // patch at world origin
                                     heightmap_file,                                       // grayscale height-map
                                     terrainWidth, terrainHeight,                          // physical X/Y extents (m)
                                     hMin, hMax,                                           // height range (m)
                                     0.001f,                                                // swept-sphere radius
                                     true);                                                // build visual mesh



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

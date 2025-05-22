// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2025 projectchrono.org
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
// As an example, a perlin noise heightmap is generated as an array.
// The Heightfield patch type uses a y=0 bottom approach, and can be used 
// with Unity/Unreal terrain types.
// Note: Visual mesh generation is the primary performance bottleneck
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

    double step_size = 2e-3;
    double tire_step_size = 2e-3;

    // Height map params
    double terrainWidth = 100.0;  // eg. scale across 100m x 100m
    double terrainHeight = 100.0;
    int heightMapNx = 512, heightMapNy = 513;  // 513x513 resolution (higher resolutions still run fast - try 2049x2049 -irrlicht visual mesh is what's slow)
    double heightAmp = 2;  // Scale the nooise-based height map to this height (m)

    // Create height map
    HeightMap heightMap(terrainWidth, terrainHeight, heightMapNx, heightMapNy, heightAmp);

    // Create vehicle
    HMMWV_Full hmmwv;
    hmmwv.SetContactMethod(ChContactMethod::NSC);
    hmmwv.SetChassisFixed(false);
    hmmwv.SetInitPosition(ChCoordsys<>(ChVector3d(5, 0, 2.75), ChQuaterniond(QuatFromAngleX(90))));
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


    auto my_obstacle = chrono_types::make_shared<ChBodyEasyBox>(1, 0.5, 1, 200, true, true, patch_mat);
    sys->Add(my_obstacle);
    my_obstacle->SetPos(ChVector3d(20, 2, 10));
    my_obstacle->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/cubetexture_wood.png"));



    // Create the terrain
    RigidTerrain terrain(sys);
    
    // time the build - note the visual is what takes the longest. Try turning off the visual for a test of the collision shape speed
    auto start = std::chrono::high_resolution_clock::now();

    // Use the entire height map for one large patch
    auto heights = heightMap.getAllHeights();
    auto dynPatch = terrain.AddPatch(patch_mat,                                 // NSC or SMC material
                                                ChCoordsys<>(ChVector3d(0, 0, 0), QUNIT),  // patch body at world origin, no rotation
                                                heightMap.getAllHeights(),                 // row major height array
                                                heightMapNx, heightMapNy,     // full grid resolution
                                                terrainWidth, terrainHeight,  // full terrain extents (Am x Bm)
                                                0.1f,
                                                true            // build chtrianglemesh visual shape
    );

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

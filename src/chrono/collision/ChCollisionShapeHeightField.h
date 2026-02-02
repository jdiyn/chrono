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

#ifndef CH_COLLISION_SHAPE_HEIGHTFIELD_H
#define CH_COLLISION_SHAPE_HEIGHTFIELD_H

#include "chrono/collision/ChCollisionShape.h"
#include <vector>

namespace chrono {

/// @addtogroup chrono_collision
/// @{

/// A heightfield collision shape for terrain represented as a 2D grid of height samples.
///
/// Under Bullet, this becomes a cbtHeightfieldChronoTerrainShape with optimized O(1) collision
/// detection for spheres and efficient support sampling for convex shapes.
///
/// ## Coordinate Conventions
/// - **Local coordinates** are centered: planar origin at (0,0), center of the heightfield
/// - **Height values** use BASE-at-origin: local height=0 corresponds to minHeight
/// - **Grid indexing** is row-major: heights[j * nx + i] where i is along width, j is along length
/// - **Up-axis** determines which local axis is height: 0=X, 1=Y, 2=Z (default Z)
///
/// ## Caching Behavior
/// For small heightfields (≤512×512), vertex caching is automatic. For larger terrains,
/// a tiled LOD cache can be enabled with SetTiledCacheConfig() to efficiently handle
/// million-vertex heightfields with O(1) collision performance.
///
/// ## Performance Notes
/// - Sphere collision: O(1) analytical bilinear interpolation + ClosestPointOnTriangle
/// - Convex collision: O(1) height rejection + support sampling (26-42 directions)
/// - Large terrain: Tiled caching with dirty region tracking
class ChApi ChCollisionShapeHeightField : public ChCollisionShape {
  public:

    /// World space query: returns true if the point projects inside the patch
    bool SampleWorld(const ChCoordsys<>& patchFrame,  // body/patch frame
                     const ChVector3d& worldPos,      // query point
                     double& heightOut,               // abs world height on HF
                     ChVector3d& normalOut) const;

    ChCollisionShapeHeightField();
    /// Construct from a regular grid of heights.
    ChCollisionShapeHeightField(std::shared_ptr<ChContactMaterial> material,    ///< [in] contact material
                                int nx,                                         ///< [in] number of height samples along X
                                int ny,                                         ///< [in] number of height samples along Y
                                double dimX,                                    ///< [in] physical width along X (meters)
                                double dimY,                                    ///< [in] physical length along Y (meters)
                                const std::vector<double>& heights_absolute,    ///< [in] height array (row-major: j*nx + i) ---- always double input, if single bullet, gets coverted in constructor
                                float heightScale,                              ///< [in] multiplier applied to raw height values
                                float minHeight,                                ///< [in] minimum height value (for AABB centring)
                                float maxHeight,                                ///< [in] maximum height value (for AABB centring)
                                int upAxis = 2,                                 ///< [in] 0 is X, 1 is Y, 2 is Z (height axis)
                                float sphere_radius = 0.001f,                   ///< [in] swept sphere radius
                                bool flipQuadEdges = true);                     ///< [in] default quad flips for robust uniformity to the mesh

    ~ChCollisionShapeHeightField() override {}

    /// Getters for internal data
    int               GetWidthSamples()   const { return m_nx; }
    int               GetLengthSamples()  const { return m_ny; }
    double            GetFieldWidth()     const { return m_width; }
    double            GetFieldLength()    const { return m_length; }
    // Return NON-CONST pointer for Bullet shape creation (in ch collision model bullet)
   // double* GetHeights() { return m_heights.data(); }              // Non-const version
    const double* GetHeights() const { return m_heights.data(); }  // ABSOLUTE HEIGHTS! NOT SHIFTED HEIGHTS!
#ifndef BT_USE_DOUBLE_PRECISION
    // if bullet is single precision
   // float* GetHeightsFloat() { return m_heights_f.data(); }
    const float* GetHeightsFloat() const { return m_heights_f.data(); }
 #endif
 
    float             GetHeightScale()    const { return m_heightScale; }
    float             GetMinHeight()      const { return m_minHeight; }
    float             GetMaxHeight()      const { return m_maxHeight; }
    int               GetUpAxis()         const { return m_upAxis; }
    bool              GetFlipQuadEdges()  const { return m_flipQuadEdges; }

    /// Override bounding box to encompass the full heightfield
    ChAABB GetBoundingBox() const override;

    /// Return the thickness as the radius of a sphere-swept mesh.
    double GetRadius() const { return sradius; }

    /// Calculate the height of the patch at the given position and return the height and normal vector.
    /// Uses bilinear interpolation matching the Bullet collision shape's internal algorithm.
    /// @param patch_frame [in] body/patch frame (position and rotation)
    /// @param query_pos [in] query position in world space
    /// @param world_up [in] world up vector (usually Z-up: ChVector3d(0,0,1))
    /// @param out_height [out] height at the query position projected onto world_up
    /// @param out_normal [out] surface normal at the query position (in world space)
    /// @return true if query projects inside patch bounds, false if outside
    bool RayHit(const ChCoordsys<>& patch_frame,    ///< [in] body/patch frame
                const ChVector3d& query_pos,        ///< [in] query position in world space
                const ChVector3d& world_up,         ///< [in] world up vector (usually Z)
                double& out_height,                 ///< [out] height at the query position (in world space)
                ChVector3d& out_normal) const;      ///< [out] normal at the query position (in world space)

    //
    // Tiled Cache Configuration (for large terrains > 512x512)
    // These settings are applied when the collision model is built.
    //

    /// Enable tiled LOD caching for large heightfields.
    /// When enabled, only tiles near active collision queries are cached at full resolution.
    /// @param tileSize Grid cells per tile edge (default 64, must be power of 2)
    /// @param cacheRadius Number of tiles to cache at full resolution around query point
    void SetTiledCacheConfig(int tileSize = 64, int cacheRadius = 2) {
        m_tiledCacheTileSize = tileSize;
        m_tiledCacheRadius = cacheRadius;
    }

    /// Get tile size for tiled cache (0 if not configured)
    int GetTiledCacheTileSize() const { return m_tiledCacheTileSize; }

    /// Get cache radius for tiled cache (0 if not configured)
    int GetTiledCacheRadius() const { return m_tiledCacheRadius; }
        
    void ArchiveOut(ChArchiveOut& archive);
    void ArchiveIn(ChArchiveIn& archive);


  private:
    ChCoordsys<> m_patchFrame;  ///< patch frame in world space
    int                     m_nx, m_ny;
    double                  m_width, m_length;
    double m_heightCentre;  // (minH+maxH)/2
    std::vector<double>     m_heights;
#ifndef BT_USE_DOUBLE_PRECISION
    std::vector<float>      m_heights_f; // float array if not using chrono's bullet double setup
#endif
    float                   m_heightScale, m_minHeight, m_maxHeight;
    int                     m_upAxis;
    bool                    m_flipQuadEdges;
    float                   sradius;

    // Precompute cell sizes (metres) and reciprocals to trasnlate to from the grid
    double m_cellSizeU;  // width  /(nx1)
     double m_cellSizeV;      // length /(ny1)
     double m_invCellSizeU;   // 1/m_cellSizeU
     double m_invCellSizeV;   // 1/m_cellSizeV

    // Tiled cache configuration (applied when collision model is built)
    int m_tiledCacheTileSize = 0;   // 0 = auto (use flat cache for small HFs)
    int m_tiledCacheRadius = 2;     // Number of tiles at full res
};


/// @} chrono_collision

}  // namespace chrono

#endif

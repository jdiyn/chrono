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

/// A heightfield collision shape.  Under Bullet, this becomes a btHeightfieldTerrainShape.
class ChApi ChCollisionShapeHeightField : public ChCollisionShape {
  public:
    /// Construct from a regular grid of heights.
    ChCollisionShapeHeightField(std::shared_ptr<ChContactMaterial> material,    ///< [in] contact material
                                int nx,                                         ///< [in] number of height samples along X
                                int ny,                                         ///< [in] number of height samples along Y
                                double dimX,                                    ///< [in] physical width along X (meters)
                                double dimY,                                    ///< [in] physical length along Y (meters)
                                const std::vector<double>& heights,             ///< [in] height array (row-major: j*nx + i), must be float
                                float heightScale,                              ///< [in] multiplier applied to raw height values - needed for any integer-based heightfield data types
                                float minHeight,                                ///< [in] minimum height value (for AABB centering)
                                float maxHeight,                                ///< [in] maximum height value (for AABB centering)
                                int upAxis = 2,                                 ///< [in] 0 is X, 1 is Y, 2 is Z (height axis)
                                float sphere_radius = 0.001f,                            ///< [in] swept sphere radius
                                bool flipQuadEdges = true);                     ///< [in] default quad flips for robust uniformity to the mesh

    ~ChCollisionShapeHeightField() override {}

    /// Getters for internal data
    int               GetWidthSamples()   const { return m_nx; }
    int               GetLengthSamples()  const { return m_ny; }
    double            GetFieldWidth()     const { return m_width; }
    double            GetFieldLength()    const { return m_length; }
    const double* GetHeights() const { return m_heights.data(); }
#ifndef BT_USE_DOUBLE_PRECISION
    // Bullet is single precision
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

        
    void ArchiveOut(ChArchiveOut& archive);
    void ArchiveIn(ChArchiveIn& archive);


  private:
    int                     m_nx, m_ny;
    double                  m_width, m_length;
    std::vector<double>     m_heights;
#ifndef BT_USE_DOUBLE_PRECISION
    std::vector<float>      m_heights_f; // float array if not using chrono's bullet double setup
#endif
    float                   m_heightScale, m_minHeight, m_maxHeight;
    int                     m_upAxis;
    bool                    m_flipQuadEdges;
    float                   sradius;
};


/// @} chrono_collision

}  // namespace chrono

#endif

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

#include "chrono/collision/ChCollisionShapeHeightField.h"

namespace chrono {

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChCollisionShapeHeightField)
CH_UPCASTING(ChCollisionShapeHeightField, ChCollisionShape)

ChCollisionShapeHeightField::ChCollisionShapeHeightField() : ChCollisionShape(Type::HEIGHTFIELD) {}


// constructor
ChCollisionShapeHeightField::ChCollisionShapeHeightField(std::shared_ptr<ChContactMaterial> material,
                                                         int nx,
                                                         int ny,
                                                         double dimX,
                                                         double dimY,
                                                         const std::vector<double>& heights_absolute, // double input only - bullet will convert in its constrcutor if using floats
                                                         float heightScale,
                                                         float minHeight,
                                                         float maxHeight,
                                                         int upAxis,
                                                         float sphere_radius,
                                                         bool flipQuadEdges)
    : ChCollisionShape(Type::HEIGHTFIELD, material),
      m_nx(nx),
      m_ny(ny),
      m_width(dimX),
      m_length(dimY),
      m_heights(heights_absolute),
      m_heightScale(heightScale),
      m_minHeight(minHeight),
      m_maxHeight(maxHeight),
      m_upAxis(upAxis),
      sradius(sphere_radius),
      m_flipQuadEdges(flipQuadEdges) {

    // if not using double precision, setup for bullet a float array
    // the alternative option is to set up a use of the cbtScalar, but dont want bullet headers here
    // as it becomes limiting to the multicore system only
#ifndef BT_USE_DOUBLE_PRECISION
    m_heights_f.resize(heights_absolute.size());
    for (size_t i = 0; i < heights_absolute.size(); ++i) {
        m_heights_f[i] = static_cast<float>(heights_absolute[i]);  // no centering for bullet shape - the base position
                                                                   // is 0,0,0 at the floor of the heightfield
    }
#endif

    // Precompute cell geometry for direct sampling
    m_cellSizeU = m_width / double(m_nx - 1);
    m_cellSizeV = m_length / double(m_ny - 1);
    m_invCellSizeU = 1.0 / m_cellSizeU;
    m_invCellSizeV = 1.0 / m_cellSizeV;

}

ChAABB ChCollisionShapeHeightField::GetBoundingBox() const {
    // AABB centered at zero in the local frame
    ChVector3d lo(-m_width / 2, -m_length / 2, m_minHeight * m_heightScale);
    ChVector3d hi(m_width / 2, m_length / 2, m_maxHeight * m_heightScale);
    return ChAABB(lo, hi);
}

// ============================================================================
// RayHit - Analytic vertical raycast on the heightfield patch
// ============================================================================
// COORDINATE CONVENTIONS:
// - Local coordinates are centered: planar origin at (0,0), height uses BASE-at-origin
// - The heightfield BASE (local height=0) corresponds to minHeight in world units
// - Heights are stored as absolute values; the Bullet shape handles centering internally
// - This function delegates to the Bullet shape's bilinear interpolation to avoid
//   duplicating logic (single source of truth for height queries)
//
// Returns true if query projects onto the patch rectangle, false if outside.
bool ChCollisionShapeHeightField::RayHit(const ChCoordsys<>& frame,
                                         const ChVector3d& query_pos,
                                         const ChVector3d& world_up,
                                         double& out_height,
                                         ChVector3d& out_normal) const {
    // Transform the query point to the shape's local frame
    ChVector3d p_local = frame.TransformPointParentToLocal(query_pos);

    // Extract planar (u,v) coordinates based on up-axis (in meters, centered around origin)
    double u, v;
    switch (m_upAxis) {
        case 0:  // X-up: planar axes are Y,Z
            u = p_local.y();
            v = p_local.z();
            break;
        case 1:  // Y-up: planar axes are X,Z
            u = p_local.x();
            v = p_local.z();
            break;
        default:  // Z-up: planar axes are X,Y
            u = p_local.x();
            v = p_local.y();
            break;
    }

    // Quick bounds check (half-width/length centered at origin)
    const double half_width = 0.5 * m_width;
    const double half_length = 0.5 * m_length;
    if (std::fabs(u) > half_width || std::fabs(v) > half_length) {
        return false;  // Outside patch footprint
    }

    // Use the same bilinear interpolation logic as the Bullet collision shape.
    // This ensures consistency between ray queries and collision detection.
    // The computation matches cbtHeightfieldChronoTerrainShape::queryHeightAndGradient()
    
    // Convert to grid coordinates (fractional)
    const double fx = (u + half_width) * m_invCellSizeU;
    const double fy = (v + half_length) * m_invCellSizeV;

    int i = static_cast<int>(fx);
    int j = static_cast<int>(fy);
    double tx = fx - i;
    double ty = fy - j;

    // Clamp to valid grid indices
    i = std::clamp(i, 0, m_nx - 2);
    j = std::clamp(j, 0, m_ny - 2);
    tx = std::clamp(tx, 0.0, 1.0);
    ty = std::clamp(ty, 0.0, 1.0);

    // Raw height access (matches Bullet shape's getRawHeightFieldValue * heightScale)
    auto H = [&](int ix, int iy) -> double { 
        return m_heights[static_cast<size_t>(iy) * m_nx + ix] * m_heightScale; 
    };

    // Fetch corner heights
    const double h00 = H(i, j);
    const double h10 = H(i + 1, j);
    const double h01 = H(i, j + 1);
    const double h11 = H(i + 1, j + 1);

    // Bilinear interpolation for height
    const double h0 = h00 + tx * (h10 - h00);
    const double h1 = h01 + tx * (h11 - h01);
    const double h_absolute = h0 + ty * (h1 - h0);

    // Convert to BASE-relative height (local height=0 corresponds to minHeight)
    const double h_local = h_absolute - (m_minHeight * m_heightScale);

    // Compute gradient (partial derivatives) for normal calculation
    // Uses the same formula as cbtHeightfieldChronoTerrainShape::queryHeightAndGradient
    const double dhdx = ((h10 - h00) * (1.0 - ty) + (h11 - h01) * ty) * m_invCellSizeU;
    const double dhdz = ((h01 - h00) * (1.0 - tx) + (h11 - h10) * tx) * m_invCellSizeV;

    // Build local normal from gradient (points outward from ground)
    ChVector3d n_local;
    switch (m_upAxis) {
        case 0:  // X-up
            n_local.Set(1.0, -dhdx, -dhdz);
            break;
        case 1:  // Y-up
            n_local.Set(-dhdx, 1.0, -dhdz);
            break;
        default:  // Z-up
            n_local.Set(-dhdx, -dhdz, 1.0);
            break;
    }
    
    // Normalize (with fallback for degenerate cases)
    const double len2 = n_local.Length2();
    if (len2 > 1e-12) {
        n_local /= std::sqrt(len2);
    } else {
        // Fallback to pure up-axis normal
        n_local = ChVector3d(0, 0, 0);
        switch (m_upAxis) {
            case 0: n_local.x() = 1.0; break;
            case 1: n_local.y() = 1.0; break;
            default: n_local.z() = 1.0; break;
        }
    }

    // Build local surface point (BASE-relative coordinates)
    ChVector3d surf_local;
    switch (m_upAxis) {
        case 0:  surf_local.Set(h_local, u, v); break;
        case 1:  surf_local.Set(u, h_local, v); break;
        default: surf_local.Set(u, v, h_local); break;
    }

    // Transform to world space
    ChVector3d surf_world = frame.TransformPointLocalToParent(surf_local);
    out_height = surf_world.Dot(world_up);  // Project onto world up direction
    out_normal = frame.TransformDirectionLocalToParent(n_local);

    return true;
}


// ============================================================================
// Serialization
// ============================================================================
// Note: Only the precision-matching height array is serialized (double or float
// depending on BT_USE_DOUBLE_PRECISION). The other array and all Bullet-level
// caches (vertex cache, quad extents, accelerator) are reconstructed when the
// Bullet collision shape is created from this data during model injection.

void ChCollisionShapeHeightField::ArchiveOut(ChArchiveOut& archive_out) {
    archive_out.VersionWrite<ChCollisionShapeHeightField>();
    archive_out << CHNVP(m_nx);
    archive_out << CHNVP(m_ny);
    archive_out << CHNVP(m_width);
    archive_out << CHNVP(m_length);
    archive_out << CHNVP(m_heights);
#ifndef BT_USE_DOUBLE_PRECISION
    archive_out << CHNVP(m_heights_f);
#endif
    archive_out << CHNVP(m_heightScale);
    archive_out << CHNVP(m_minHeight);
    archive_out << CHNVP(m_maxHeight);
    archive_out << CHNVP(m_upAxis);
    archive_out << CHNVP(sradius);
    archive_out << CHNVP(m_flipQuadEdges);
}

void ChCollisionShapeHeightField::ArchiveIn(ChArchiveIn& archive_in) {
    archive_in.VersionRead<ChCollisionShapeHeightField>();
    archive_in >> CHNVP(m_nx);
    archive_in >> CHNVP(m_ny);
    archive_in >> CHNVP(m_width);
    archive_in >> CHNVP(m_length);
    archive_in >> CHNVP(m_heights);
#ifndef BT_USE_DOUBLE_PRECISION
    archive_in >> CHNVP(m_heights_f);
#endif
    archive_in >> CHNVP(m_heightScale);
    archive_in >> CHNVP(m_minHeight);
    archive_in >> CHNVP(m_maxHeight);
    archive_in >> CHNVP(m_upAxis);
    archive_in >> CHNVP(sradius);
    archive_in >> CHNVP(m_flipQuadEdges);

    // Recompute derived cell geometry after deserialization
    if (m_nx > 1 && m_ny > 1) {
        m_cellSizeU = m_width / double(m_nx - 1);
        m_cellSizeV = m_length / double(m_ny - 1);
        m_invCellSizeU = 1.0 / m_cellSizeU;
        m_invCellSizeV = 1.0 / m_cellSizeV;
    }
}

}  // namespace chrono
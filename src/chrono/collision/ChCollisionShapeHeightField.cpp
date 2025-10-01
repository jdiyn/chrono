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
    : ChCollisionShape(Type::UNKNOWN_SHAPE, material),
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
    m_type = Type::HEIGHTFIELD;

        // store the vertical centre that bullet does so SampleHeight can restore absolute heights correctly
   // m_heightCentre = 0.5 * (m_minHeight + m_maxHeight);

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

////////////// -- TODO --- possibly no longer needed?
//  Analytic "vertical rayhit" on a heightfield patch
//  returns true - query projects onto the patch rectangle or false - outside
bool ChCollisionShapeHeightField::RayHit(const ChCoordsys<>& frame,
                                         const ChVector3d& query_pos,
                                         const ChVector3d& world_up,
                                         double& out_height,
                                         ChVector3d& out_normal) const {
    // Transform the query point to the shape local frame
    ChVector3d p_local = frame.TransformPointParentToLocal(query_pos);

    // Convert to planar (u,v) coordinates according to up-axis (in meters, centered)
    double u, v;
    switch (m_upAxis) {
        case 0:  // X-up
            u = p_local.y();
            v = p_local.z();
            break;
        case 1:  // Y-up
            u = p_local.x();
            v = p_local.z();
            break;
        default:  // Z-up
            u = p_local.x();
            v = p_local.y();
            break;
    }

    // Check if query projects inside the patch rectangle (half-width/length bounds)
    double half_width = 0.5 * m_width;
    double half_length = 0.5 * m_length;
    if (std::fabs(u) > half_width || std::fabs(v) > half_length) {
        return false;  // Outside patch
    }

    // Convert to grid coordinates (fractional)
    double fx = (u + half_width) * m_invCellSizeU;
    double fy = (v + half_length) * m_invCellSizeV;

    int i = static_cast<int>(fx);
    int j = static_cast<int>(fy);
    double tx = fx - i;
    double ty = fy - j;

    // Clamp to valid grid indices (prevent out-of-bounds)
    i = std::clamp(i, 0, m_nx - 2);
    j = std::clamp(j, 0, m_ny - 2);
    if (tx > 1.0)
        tx = 1.0;
    if (ty > 1.0)
        ty = 1.0;

    // Lambda for raw height access (scaled by m_heightScale, which is 1.0)
    auto H = [&](int ix, int iy) { return m_heights[iy * m_nx + ix] * m_heightScale; };

    // Fetch corner heights (raw, not centered)
    const double h00 = H(i, j);
    const double h10 = H(i + 1, j);
    const double h01 = H(i, j + 1);
    const double h11 = H(i + 1, j + 1);

    // Bilinear interpolation for height (raw height)
    const double h0 = h00 + tx * (h10 - h00);
    const double h1 = h01 + tx * (h11 - h01);
    const double h = h0 + ty * (h1 - h0);

    // Gradients for normal (central differences, scaled by cell size)
    const double dhdu = (h10 - h00 + h11 - h01) * 0.5 * m_invCellSizeU;
    const double dhdv = (h01 - h00 + h11 - h10) * 0.5 * m_invCellSizeV;

    // Build local normal (points out of the ground)
    ChVector3d n_local(0.0, 0.0, 0.0);
    const double epsilon = 1e-8;  // Threshold for flat areas
    if (std::fabs(dhdu) < epsilon && std::fabs(dhdv) < epsilon) {
        // Perfectly flat: force normal along up-axis
        switch (m_upAxis) {
            case 0:
                n_local.Set(1.0, 0.0, 0.0);
                break;
            case 1:
                n_local.Set(0.0, 1.0, 0.0);
                break;
            default:
                n_local.Set(0.0, 0.0, 1.0);
                break;
        }
    } else {
        // Non-flat: compute gradient-based normal
        switch (m_upAxis) {
            case 0:
                n_local.Set(1.0, -dhdu, -dhdv);
                break;
            case 1:
                n_local.Set(-dhdu, 1.0, -dhdv);
                break;
            default:
                n_local.Set(-dhdu, -dhdv, 1.0);
                break;
        }
        n_local.Normalize();
    }

    // Construct local surface point (using raw height h, no centering needed here)
    ChVector3d surf_local(0.0, 0.0, 0.0);
    switch (m_upAxis) {
        case 0:
            surf_local.Set(h, u, v);
            break;
        case 1:
            surf_local.Set(u, h, v);
            break;
        default:
            surf_local.Set(u, v, h);
            break;
    }

    // Transform to world space
    ChVector3d surf_world = frame.TransformPointLocalToParent(surf_local);
    out_height = surf_world.Dot(world_up);  // Project onto world up (handles arbitrary up-axis)
    out_normal = frame.TransformDirectionLocalToParent(n_local);

    return true;
}


void ChCollisionShapeHeightField::ArchiveOut(ChArchiveOut& archive_out) {
    archive_out.VersionWrite<ChCollisionShapeHeightField>();
    archive_out << CHNVP(m_nx);
    archive_out << CHNVP(m_ny);
    archive_out << CHNVP(m_width);
    archive_out << CHNVP(m_length);
#ifdef BT_USE_DOUBLE_PRECISION
    archive_out << CHNVP(m_heights);
#else
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
#ifdef BT_USE_DOUBLE_PRECISION
    archive_in >> CHNVP(m_heights);
#else
    archive_in << CHNVP(m_heights_f);
#endif
    archive_in >> CHNVP(m_heightScale);
    archive_in >> CHNVP(m_minHeight);
    archive_in >> CHNVP(m_maxHeight);
    archive_in >> CHNVP(m_upAxis);
    archive_in >> CHNVP(sradius);
    archive_in >> CHNVP(m_flipQuadEdges);
}

}  // namespace chrono
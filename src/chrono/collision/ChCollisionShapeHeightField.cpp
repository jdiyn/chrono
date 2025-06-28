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
                                                         const std::vector<double>& heights, // double input only - bullet will convert in its constrcutor if using floats
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
      m_heights(heights),
      m_heightScale(heightScale),
      m_minHeight(minHeight),
      m_maxHeight(maxHeight),
      m_upAxis(upAxis),
      sradius(sphere_radius),
      m_flipQuadEdges(flipQuadEdges) {
    m_type = Type::HEIGHTFIELD;

    // if not using double precision, setup for bullet a float array
    // the alternative option is to set up a use of the cbtScalar, but dont want bullet headers here
    // as it becomes limiting to the multicore system only
#ifndef BT_USE_DOUBLE_PRECISION
    m_heights_f.resize(heights.size());
    for (size_t i = 0; i < heights.size(); ++i) {
        m_heights_f[i] = static_cast<float>(heights[i]);
    }
#endif

    // Precompute cell geometry for direct sampling
    m_cellSizeU = m_width / double(m_nx - 1);
    m_cellSizeV = m_length / double(m_ny - 1);
    m_invCellSizeU = 1.0 / m_cellSizeU;
    m_invCellSizeV = 1.0 / m_cellSizeV;

    // store the vertical centre that bullet does so SampleHeight can restore absolute heights correctly
    m_heightCentre = 0.5 * (m_minHeight + m_maxHeight);
}

ChAABB ChCollisionShapeHeightField::GetBoundingBox() const {
    // AABB centered at zero in the local frame
    ChVector3d lo(-m_width / 2, -m_length / 2, m_minHeight * m_heightScale);
    ChVector3d hi(m_width / 2, m_length / 2, m_maxHeight * m_heightScale);
    return ChAABB(lo, hi);
}


//  Analytic "vertical rayhit" on a heightfield patch
//  returns true - query projects onto the patch rectangle or false - outside
bool ChCollisionShapeHeightField::RayHit(const ChCoordsys<>& frame,
                                         const ChVector3d& query_pos,
                                         const ChVector3d& world_up,
                                         double& out_height,
                                         ChVector3d& out_normal) const {
    // Transform the query point to the shape local frame
    ChVector3d p_local = frame.TransformPointParentToLocal(query_pos);

    // NB:  In a purely Chrono context the heightfield is not subject to an
    // additional bulletstyle local scaling, so there's no need for inverse
    // scaling like in the bullet class. Height field values are stored directly

    // Convert to planar (u,v) coordinates according to upaxis
    double u, v;  // metres in the heightfield plane
    switch (m_upAxis) {
        case 0:
            u = p_local.y();
            v = p_local.z();
            break;
        case 1:
            u = p_local.x();
            v = p_local.z();
            break;
        default:
            u = p_local.x();
            v = p_local.y();
            break;
    }

    // Rectangle bounds check (halfwidth / halflength)
    // false if we're outside the patch
    if (std::fabs(u) > 0.5 * m_width || std::fabs(v) > 0.5 * m_length)
        return false;

    // bilinear height evaluation
    const double fx = (u + 0.5 * m_width) * m_invCellSizeU;  // grid coord
    const double fy = (v + 0.5 * m_length) * m_invCellSizeV;

    int i = (int)fx;
    int j = (int)fy;
    double tx = fx - i;
    double ty = fy - j;

    // clamp on right/top border
    if (i >= m_nx - 1) {
        i = m_nx - 2;
        tx = 1.0;
    }
    if (j >= m_ny - 1) {
        j = m_ny - 2;
        ty = 1.0;
    }

    auto H = [&](int ix, int iy) { return m_heights[iy * m_nx + ix] * m_heightScale; };

    const double h00 = H(i, j);
    const double h10 = H(i + 1, j);
    const double h01 = H(i, j + 1);
    const double h11 = H(i + 1, j + 1);

    const double h0 = h00 + tx * (h10 - h00);
    const double h1 = h01 + tx * (h11 - h01);
    const double h = h0 + ty * (h1 - h0);  // bilinear height

    // Gradient dH/du, dH/dv  (cell constant central difference). remember the inverse cell size handled in constructor
    // cannot be changed afterwards
    const double dhdu = (h10 - h00 + h11 - h01) * 0.5 * m_invCellSizeU;
    const double dhdv = (h01 - h00 + h11 - h10) * 0.5 * m_invCellSizeV;

    // Build local normal (unscaled, because Chrono HF stores height directly)
    ChVector3d n_local;
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

    // Return worldspace height & normal
    const double centre = 0.5 * (m_minHeight + m_maxHeight);  // Chrono/Bullet convention
    // convert local surface point to world
    ChVector3d surf_local(0, 0, 0);
    switch (m_upAxis) {
        case 0:
            surf_local.Set(h - centre, u, v);
            break;  // Xup
        case 1:
            surf_local.Set(u, h - centre, v);
            break;  // Yup
        default:
            surf_local.Set(u, v, h - centre);
            break;  // Zup
    }
    ChVector3d surf_world = frame.TransformPointLocalToParent(surf_local);

    out_height = world_up ^ surf_world;  // dot product duplication of the chrono_vehicle chworldframe::height
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
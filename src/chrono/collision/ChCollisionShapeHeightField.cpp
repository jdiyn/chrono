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

// constructor
ChCollisionShapeHeightField::ChCollisionShapeHeightField(std::shared_ptr<ChContactMaterial> material,
                                                         int nx,
                                                         int ny,
                                                         double dimX,
                                                         double dimY,
                                                         const std::vector<double>& heights,
                                                         float heightScale,
                                                         float minHeight,
                                                         float maxHeight,
                                                         int upAxis,
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
}

ChAABB ChCollisionShapeHeightField::GetBoundingBox() const {
    // AABB centered at zero in the local frame
    ChVector3d lo(-m_width / 2, -m_length / 2, m_minHeight * m_heightScale);
    ChVector3d hi(m_width / 2, m_length / 2, m_maxHeight * m_heightScale);
    return ChAABB(lo, hi);
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
    archive_in >> CHNVP(m_flipQuadEdges);
}

}  // namespace chrono
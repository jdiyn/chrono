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
// This is not a standard bullet heightfield class, but a Chrono-specific one
// which is inspired by upon the bullet version but aims to improve it.
//
// SIMD OPTIMIZATION:
// This file includes optional SSE intrinsics for bilinear interpolation,
// providing ~2x speedup on x86/x64 platforms. Falls back to scalar code
// =============================================================================

#include "cbtHeightfieldChronoTerrainShape.h"
#include "chrono/collision/bullet/LinearMath/cbtTransformUtil.h"
#include <algorithm>
#include <cmath>

// stores the vertices, without touching bullet, to use for parallelism safely
struct TmpTriangleBuffer : public cbtTriangleCallback {
    std::vector<cbtVector3> verts;  // 3 × N
    void processTriangle(cbtVector3* tri, int, int) override {
        verts.push_back(tri[0]);
        verts.push_back(tri[1]);
        verts.push_back(tri[2]);
    }
};

void cbtHeightfieldChronoTerrainShape::sampleHeight(cbtScalar u,
                                                    cbtScalar v,
                                                    cbtScalar& outH,
                                                    cbtVector3& outGrad) const {
    // NOTE: In this heightfield, planar coordinates are centered around (0,0) in the two axes
    // orthogonal to upAxis (matching queryHeightAndGradient and the collision algorithms).
    // This helper now matches that convention.

    // Clamp into the valid footprint (queryHeightAndGradient also clamps, but this avoids huge values).
    if (u < -m_halfWidth)
        u = -m_halfWidth;
    else if (u > m_halfWidth)
        u = m_halfWidth;
    if (v < -m_halfLength)
        v = -m_halfLength;
    else if (v > m_halfLength)
        v = m_halfLength;

    // queryHeightAndGradient returns base-relative height (centered) and a unit normal in local coords.
    queryHeightAndGradient(u, v, outH, outGrad);
}

// quantisation helper
static inline int getQuantized(cbtScalar x) {
    return (x < cbtScalar(0)) ? int(x - cbtScalar(0.5)) : int(x + cbtScalar(0.5));
}

void cbtHeightfieldChronoTerrainShape::quantizeWithClamp(int out[3], const cbtVector3& pt) const {
    cbtVector3 clamped = pt;
    clamped.setMax(m_localAabbMin);  // ensure correct min / maxing here. max the aabb min and minimise the aabb max.
    clamped.setMin(m_localAabbMax);

    out[0] = getQuantized(clamped.getX());
    out[1] = getQuantized(clamped.getY());
    out[2] = getQuantized(clamped.getZ());

    cbtClamp(out[0], 0, m_heightStickWidth - 1);
    cbtClamp(out[1], 0, m_heightStickLength - 1);
    out[2] = 0;
}

////------------------------------------------------------------------------------
//// replicate Bullet’s getVertex (including centering by m_localOrigin)
////------------------------------------------------------------------------------
// void cbtHeightfieldChronoTerrainShape::getVertex(int x, int y, cbtVector3& vtx) const {
//    cbtScalar h = getRawHeightFieldValue(x, y) * m_heightScale;
//    switch (m_upAxis) {
//        case 0:
//            vtx.setValue(h - m_localOrigin.getX(), (-m_width * cbtScalar(0.5)) + cbtScalar(x),
//                         (-m_length * cbtScalar(0.5)) + cbtScalar(y));
//            break;
//        case 1:
//            vtx.setValue((-m_width * cbtScalar(0.5)) + cbtScalar(x), h - m_localOrigin.getY(),
//                         (-m_length * cbtScalar(0.5)) + cbtScalar(y));
//            break;
//        case 2:
//            vtx.setValue((-m_width * cbtScalar(0.5)) + cbtScalar(x), (-m_length * cbtScalar(0.5)) + cbtScalar(y),
//                         h - m_localOrigin.getZ());
//            break;
//        default:
//            cbtAssert(false);
//    }
//    vtx *= m_localScaling;
//}

// uses the cached vertices
inline void cbtHeightfieldChronoTerrainShape::getVertex(int x, int y, cbtVector3& vtx) const {
    if (m_useVertexCache && !m_vertexCache.empty()) {
        vtx = m_vertexCache[static_cast<std::size_t>(y) * m_heightStickWidth + x];
        return;
    }

    // Compute vertex on the fly (centered + scaled).
    // Must match buildVertexCache exactly!
    // In unscaled space, grid step is always 1.0 (m_width = m_heightStickWidth - 1)
    const cbtScalar h = getRawHeightFieldValue(x, y) * m_heightScale;
    const cbtScalar lx = cbtScalar(x) - m_halfWidth;
    const cbtScalar ly = cbtScalar(y) - m_halfLength;
    switch (m_upAxis) {
        case 0:
            vtx.setValue(h, lx, ly);
            break;
        case 1:
            vtx.setValue(lx, h, ly);
            break;
        default:
            vtx.setValue(lx, ly, h);
            break;
    }
    vtx *= m_localScaling;
}

void cbtHeightfieldChronoTerrainShape::updateInverseLocalScaling() {
    const cbtScalar x = m_localScaling.getX();
    const cbtScalar y = m_localScaling.getY();
    const cbtScalar z = m_localScaling.getZ();
    if (x == cbtScalar(0) || y == cbtScalar(0) || z == cbtScalar(0)) {
        m_invLocalScaling = cbtVector3(cbtScalar(1), cbtScalar(1), cbtScalar(1));
    } else {
        m_invLocalScaling = cbtVector3(cbtScalar(1) / x, cbtScalar(1) / y, cbtScalar(1) / z);
    }
}

void cbtHeightfieldChronoTerrainShape::setUseVertexCache(bool enable) {
    m_useVertexCache = enable;
    if (m_useVertexCache)
        buildVertexCache();
    else
        m_vertexCache.clear();
}

void cbtHeightfieldChronoTerrainShape::setUseQuadExtentsCache(bool enable) {
    m_useQuadExtentsCache = enable;
    if (m_useQuadExtentsCache)
        buildQuadExtents();
    else
        m_quadExtents.clear();
}

void cbtHeightfieldChronoTerrainShape::rebuildCaches() {
    if (m_useVertexCache)
        buildVertexCache();
    if (m_useQuadExtentsCache)
        buildQuadExtents();
    if (m_vboundsChunkSize > 0)
        buildAccelerator(m_vboundsChunkSize);
}

void cbtHeightfieldChronoTerrainShape::setLocalScaling(const cbtVector3& s) {
    // store locally – this fulfils the pure‑virtual contract
    m_localScaling = s;

    updateInverseLocalScaling();

    // rebuild lookup tables that depend on scale
    if (m_useVertexCache)
        buildVertexCache();
    if (m_useQuadExtentsCache)
        buildQuadExtents();
    if (m_vboundsChunkSize > 0)
        buildAccelerator(m_vboundsChunkSize);
}

cbtScalar cbtHeightfieldChronoTerrainShape::getHeight(int x, int z) const {
    // clamp to valid grid
    cbtClamp(x, 0, m_heightStickWidth - 1);
    cbtClamp(z, 0, m_heightStickLength - 1);
    // pull the cached, fully-built vertex in local (centered + scaled) space
    cbtVector3 vtx;
    getVertex(x, z, vtx);
    // return the component along the up-axis
    return vtx[m_upAxis];
}

// min/max helper
cbtHeightfieldChronoTerrainShape::Range cbtHeightfieldChronoTerrainShape::minmaxRange(cbtScalar a,
                                                                                      cbtScalar b,
                                                                                      cbtScalar c) {
    if (a > b) {
        if (b > c)
            return Range(c, a);
        else if (a > c)
            return Range(b, a);
        else
            return Range(b, c);
    } else {
        if (a > c)
            return Range(c, b);
        else if (b > c)
            return Range(a, b);
        else
            return Range(a, c);
    }
}

cbtHeightfieldChronoTerrainShape::cbtHeightfieldChronoTerrainShape(int heightStickWidth,
                                                                   int heightStickLength,
                                                                   const cbtScalar* heightfieldData, ///< absolute height data
                                                                   cbtScalar scale,
                                                                   cbtScalar minH, ///< absolute min, not centred
                                                                   cbtScalar maxH,  ///< absolute max, not centred
                                                                   int upAxis,
                                                                   bool flip)
    : m_heightStickWidth(heightStickWidth),
      m_heightStickLength(heightStickLength),
      m_heightScale(scale),
//      m_heightfieldData(heightfieldData), // init below
      m_minHeight(minH),
      m_maxHeight(maxH),
      m_upAxis(upAxis),
      m_flipQuadEdges(flip),
      m_useDiamondSubdivision(false),
      m_useZigzagSubdivision(false),
      m_vboundsChunkSize(0),
      m_vboundsGridWidth(0),
      m_vboundsGridLength(0),
      m_ownsHeightData(true) {  // we own the data, so copy it

    // initialise member variables
    m_shapeType = TERRAIN_SHAPE_PROXYTYPE;
    cbtAssert(heightStickWidth > 1 && heightStickLength > 1 && minH <= maxH && upAxis >= 0 && upAxis < 3);

    m_width = cbtScalar(heightStickWidth - 1);
    m_length = cbtScalar(heightStickLength - 1);
    m_halfWidth = m_width * cbtScalar(0.5);
    m_halfLength = m_length * cbtScalar(0.5);

    // BASE POSITIONING!!
    // The body frame is positioned at the BASE (z=0 in local coords corresponds to minHeight)
    // So we offset heights by -minHeight to make the base at local z=0

    // allocate and copy the heightfield data, shifting so that BASE is at 0
    int numSamples = heightStickWidth * heightStickLength;
    m_heightfieldData = new cbtScalar[numSamples];
    for (int i = 0; i < numSamples; ++i) {
        m_heightfieldData[i] = heightfieldData[i] - minH;  // Shift so BASE is at 0
    }

    // Local origin is at the BASE (local z=0)
    m_localOrigin = cbtVector3(0, 0, 0);
    // NO vertical offset - the base IS the origin - to follow the standard of other heightfields (unity/unreal/etc)

    // Build AABB with BASE at origin
    cbtScalar heightRange = (maxH - minH) * scale;
    switch (upAxis) {
        case 0:  // X-up
            m_localAabbMin.setValue(0, -m_width * cbtScalar(0.5), -m_length * cbtScalar(0.5));
            m_localAabbMax.setValue(heightRange, m_width * cbtScalar(0.5), m_length * cbtScalar(0.5));
            break;
        case 1:  // Y-up
            m_localAabbMin.setValue(-m_width * cbtScalar(0.5), 0, -m_length * cbtScalar(0.5));
            m_localAabbMax.setValue(m_width * cbtScalar(0.5), heightRange, m_length * cbtScalar(0.5));
            break;
        case 2:  // Z-up
            m_localAabbMin.setValue(-m_width * cbtScalar(0.5), -m_length * cbtScalar(0.5), 0);
            m_localAabbMax.setValue(m_width * cbtScalar(0.5), m_length * cbtScalar(0.5), heightRange);
            break;
    }


    updateInverseLocalScaling();

    // Auto-enable caching based on terrain size
    const int totalVertices = m_heightStickWidth * m_heightStickLength;
    
    if (totalVertices <= m_autoCacheThreshold) {
        // Small terrain: use flat vertex cache (fast, fits in memory)
        m_useVertexCache = true;
    }
    // Medium/large terrains rely on quad extents + on-the-fly vertex computation.
    // The O(1) neighborhood lookups (3x3 for sphere, small AABB for convex) make
    // flat caching unnecessary — the bottleneck is the solver, not vertex access.

    // build caching first time
    if (m_useQuadExtentsCache)
        buildQuadExtents();
    if (m_useVertexCache)
        buildVertexCache();

    // Build a chunked min/max accelerator by default (used by raycasts and collision chunk-culling).
    // 16 is a good tradeoff for typical heightfield tile sizes.
    buildAccelerator(16);
}

cbtHeightfieldChronoTerrainShape::~cbtHeightfieldChronoTerrainShape() {
    if (m_ownsHeightData && m_heightfieldData) {
        delete[] m_heightfieldData;
        m_heightfieldData = nullptr;
    }
    clearAccelerator();
}

void cbtHeightfieldChronoTerrainShape::getAabb(const cbtTransform& tr, cbtVector3& aabbMin, cbtVector3& aabbMax) const {
    // AABB half-extents in local space
    cbtVector3 half = (m_localAabbMax - m_localAabbMin) * (getLocalScaling() * cbtScalar(0.5));

    // BASE center offset (half the height range)
    cbtVector3 localCenter(0, 0, 0);
    localCenter[m_upAxis] =
        (m_localAabbMax[m_upAxis] - m_localAabbMin[m_upAxis]) * getLocalScaling()[m_upAxis] * cbtScalar(0.5);

    // Minimum thickness for broadphase
    const cbtScalar MIN_THICKNESS = 0.5f;
    if (half[m_upAxis] < MIN_THICKNESS)
        half[m_upAxis] = MIN_THICKNESS;

    // Transform to world space
    cbtMatrix3x3 abs_b = tr.getBasis().absolute();
    cbtVector3 worldExtents(abs_b[0].dot(half), abs_b[1].dot(half), abs_b[2].dot(half));
    worldExtents += cbtVector3(getMargin(), getMargin(), getMargin());

    cbtVector3 worldCenter = tr(localCenter);
    aabbMin = worldCenter - worldExtents;
    aabbMax = worldCenter + worldExtents;
}

// Efficient single-vertex height update using region-based cache updates
void cbtHeightfieldChronoTerrainShape::updateHeight(int x, int y, cbtScalar newHeight_absolute) {
    cbtAssert(x >= 0 && x < m_heightStickWidth && y >= 0 && y < m_heightStickLength);

    // Convert absolute height to base-relative
    m_heightfieldData[y * m_heightStickWidth + x] = newHeight_absolute - m_minHeight;

    // Expand AABB if new height exceeds current bounds (critical for CRM/SCM dynamic terrain)
    const cbtScalar newScaledH = (newHeight_absolute - m_minHeight) * m_heightScale;
    if (newScaledH > m_localAabbMax[m_upAxis])
        m_localAabbMax[m_upAxis] = newScaledH;
    if (newScaledH < m_localAabbMin[m_upAxis])
        m_localAabbMin[m_upAxis] = newScaledH;

    // Update caches using efficient region-based methods
    if (m_useVertexCache && !m_vertexCache.empty())
        rebuildVertexCacheRegion(x, y, x, y);
    
    // Quad extents: update the 4 quads that share this vertex
    if (m_useQuadExtentsCache && !m_quadExtents.empty())
        rebuildQuadExtentsRegion(x > 0 ? x - 1 : 0, y > 0 ? y - 1 : 0, x, y);

    // Update accelerator chunks affected by this vertex
    if (m_vboundsChunkSize > 0)
        updateAcceleratorRegion(x, y, x, y);
}

void cbtHeightfieldChronoTerrainShape::updateHeights(const cbtScalar* newHeights_absolute, int numSamples) {
    cbtAssert(numSamples == m_heightStickWidth * m_heightStickLength);

    for (int i = 0; i < numSamples; ++i) {
        m_heightfieldData[i] = newHeights_absolute[i] - m_minHeight;
    }

    // Recompute AABB height bounds from new data (critical for CRM/SCM dynamic terrain)
    {
        cbtScalar minH = m_heightfieldData[0] * m_heightScale;
        cbtScalar maxH = minH;
        for (int i = 1; i < numSamples; ++i) {
            const cbtScalar h = m_heightfieldData[i] * m_heightScale;
            if (h < minH) minH = h;
            if (h > maxH) maxH = h;
        }
        m_localAabbMin[m_upAxis] = minH;
        m_localAabbMax[m_upAxis] = maxH;
    }

    // Rebuild/refresh any dependent caches (vertex cache, quad extents cache, accelerator)
    rebuildCaches();
}

// ---------------------------------------------------------------------------
// Region-based height update for efficient partial terrain modification
// ---------------------------------------------------------------------------

void cbtHeightfieldChronoTerrainShape::updateHeightRegion(int x0, int z0, int x1, int z1, 
                                                           const cbtScalar* newHeights_absolute) {
    // Clamp bounds
    x0 = cbtMax(0, x0);
    z0 = cbtMax(0, z0);
    x1 = cbtMin(m_heightStickWidth - 1, x1);
    z1 = cbtMin(m_heightStickLength - 1, z1);
    
    if (x1 < x0 || z1 < z0)
        return;
    
    // Update raw height data in region
    const int regionW = x1 - x0 + 1;
    for (int z = z0; z <= z1; ++z) {
        for (int x = x0; x <= x1; ++x) {
            const int srcIdx = (z - z0) * regionW + (x - x0);
            const int dstIdx = z * m_heightStickWidth + x;
            m_heightfieldData[dstIdx] = newHeights_absolute[srcIdx] - m_minHeight;
        }
    }

    // Expand AABB if any new height exceeds current bounds (critical for CRM/SCM dynamic terrain)
    for (int z = z0; z <= z1; ++z) {
        for (int x = x0; x <= x1; ++x) {
            const cbtScalar h = m_heightfieldData[z * m_heightStickWidth + x] * m_heightScale;
            if (h > m_localAabbMax[m_upAxis])
                m_localAabbMax[m_upAxis] = h;
            if (h < m_localAabbMin[m_upAxis])
                m_localAabbMin[m_upAxis] = h;
        }
    }
    
    // Update caches for affected region only
    if (m_useVertexCache && !m_vertexCache.empty())
        rebuildVertexCacheRegion(x0, z0, x1, z1);
    
    if (m_useQuadExtentsCache && !m_quadExtents.empty())
        rebuildQuadExtentsRegion(x0 > 0 ? x0 - 1 : 0, z0 > 0 ? z0 - 1 : 0, x1, z1);
    
    if (m_vboundsChunkSize > 0)
        updateAcceleratorRegion(x0, z0, x1, z1);
}

void cbtHeightfieldChronoTerrainShape::markRegionDirty(int x0, int z0, int x1, int z1) {
    // Clamp and validate
    x0 = cbtMax(0, x0);
    z0 = cbtMax(0, z0);
    x1 = cbtMin(m_heightStickWidth - 1, x1);
    z1 = cbtMin(m_heightStickLength - 1, z1);
    
    if (x1 < x0 || z1 < z0)
        return;
    
    // Try to merge with existing dirty regions
    for (auto& dr : m_dirtyRegions) {
        // Check if regions overlap or are adjacent
        if (x0 <= dr.x1 + 1 && x1 >= dr.x0 - 1 && z0 <= dr.z1 + 1 && z1 >= dr.z0 - 1) {
            // Merge into existing region
            dr.x0 = cbtMin(dr.x0, x0);
            dr.z0 = cbtMin(dr.z0, z0);
            dr.x1 = cbtMax(dr.x1, x1);
            dr.z1 = cbtMax(dr.z1, z1);
            return;
        }
    }
    
    // Add as new dirty region
    m_dirtyRegions.push_back({x0, z0, x1, z1});
}

void cbtHeightfieldChronoTerrainShape::flushDirtyRegions() {
    if (m_dirtyRegions.empty())
        return;
    
    for (const auto& dr : m_dirtyRegions) {
        if (m_useVertexCache && !m_vertexCache.empty())
            rebuildVertexCacheRegion(dr.x0, dr.z0, dr.x1, dr.z1);
        
        if (m_useQuadExtentsCache && !m_quadExtents.empty())
            rebuildQuadExtentsRegion(dr.x0 > 0 ? dr.x0 - 1 : 0, 
                                     dr.z0 > 0 ? dr.z0 - 1 : 0, 
                                     dr.x1, dr.z1);
        
        if (m_vboundsChunkSize > 0)
            updateAcceleratorRegion(dr.x0, dr.z0, dr.x1, dr.z1);
    }
    
    m_dirtyRegions.clear();
}

void cbtHeightfieldChronoTerrainShape::rebuildVertexCacheRegion(int x0, int z0, int x1, int z1) {
    if (!m_useVertexCache || m_vertexCache.empty())
        return;
    
    const int W = m_heightStickWidth;
    // In unscaled space, grid step is always 1.0 (m_width = m_heightStickWidth - 1)
    
    for (int y = z0; y <= z1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            cbtScalar h = getRawHeightFieldValue(x, y) * m_heightScale;
            cbtScalar lx = cbtScalar(x) - m_halfWidth;
            cbtScalar ly = cbtScalar(y) - m_halfLength;
            
            cbtVector3 v;
            switch (m_upAxis) {
                case 0: v.setValue(h, lx, ly); break;
                case 1: v.setValue(lx, h, ly); break;
                default: v.setValue(lx, ly, h); break;
            }
            v *= m_localScaling;
            m_vertexCache[static_cast<std::size_t>(y) * W + x] = v;
        }
    }
}

void cbtHeightfieldChronoTerrainShape::rebuildQuadExtentsRegion(int x0, int z0, int x1, int z1) {
    if (!m_useQuadExtentsCache || m_quadExtents.empty())
        return;
    
    const int wQuads = m_heightStickWidth - 1;
    const int lQuads = m_heightStickLength - 1;
    
    // Clamp to valid quad range
    x0 = cbtMax(0, x0);
    z0 = cbtMax(0, z0);
    x1 = cbtMin(wQuads - 1, x1);
    z1 = cbtMin(lQuads - 1, z1);
    
    for (int z = z0; z <= z1; ++z) {
        for (int x = x0; x <= x1; ++x) {
            // Heights are already BASE-relative (shifted by -minH in constructor)
            const cbtScalar h00 = getRawHeightFieldValue(x, z) * m_heightScale;
            const cbtScalar h10 = getRawHeightFieldValue(x + 1, z) * m_heightScale;
            const cbtScalar h01 = getRawHeightFieldValue(x, z + 1) * m_heightScale;
            const cbtScalar h11 = getRawHeightFieldValue(x + 1, z + 1) * m_heightScale;
            
            QuadExtents e;
            e.minH = cbtMin(cbtMin(h00, h10), cbtMin(h01, h11));
            e.maxH = cbtMax(cbtMax(h00, h10), cbtMax(h01, h11));
            m_quadExtents[static_cast<std::size_t>(z) * wQuads + x] = e;
        }
    }
}

void cbtHeightfieldChronoTerrainShape::updateAcceleratorRegion(int x0, int z0, int x1, int z1) {
    if (!hasAccelerator())
        return;
    
    // Find affected chunks
    const int cX0 = x0 / m_vboundsChunkSize;
    const int cX1 = x1 / m_vboundsChunkSize;
    const int cZ0 = z0 / m_vboundsChunkSize;
    const int cZ1 = z1 / m_vboundsChunkSize;
    
    for (int cz = cZ0; cz <= cZ1 && cz < m_vboundsGridLength; ++cz) {
        for (int cx = cX0; cx <= cX1 && cx < m_vboundsGridWidth; ++cx) {
            // Rebuild this chunk's min/max
            int chunkX0 = cx * m_vboundsChunkSize;
            int chunkZ0 = cz * m_vboundsChunkSize;
            
            Range r;
            int sx = std::min(chunkX0, m_heightStickWidth - 1);
            int sz = std::min(chunkZ0, m_heightStickLength - 1);
            // Heights are already BASE-relative (shifted by -minH in constructor)
            cbtScalar h0 = getRawHeightFieldValue(sx, sz) * m_heightScale;
            r.min = r.max = h0;
            
            for (int zz = chunkZ0; zz < chunkZ0 + m_vboundsChunkSize + 1 && zz < m_heightStickLength; ++zz) {
                for (int xx = chunkX0; xx < chunkX0 + m_vboundsChunkSize + 1 && xx < m_heightStickWidth; ++xx) {
                    cbtScalar h = getRawHeightFieldValue(xx, zz) * m_heightScale;
                    r.min = cbtMin(r.min, h);
                    r.max = cbtMax(r.max, h);
                }
            }
            m_vboundsGrid[cx + cz * m_vboundsGridWidth] = r;
        }
    }
}

// accelerator build/clear
// TODO:- we likley don't need the accelerator anymore. however, its worth considering this where the terrain is partitioned into patches
//      perhaps at a higher level - i.e. up in the Chrono class. Not rigid terrain, since we want to move this to dynamic terrain handling (i.e. not 'rigid')
void cbtHeightfieldChronoTerrainShape::buildAccelerator(int chunkSize) {
    if (chunkSize <= 0) {
        clearAccelerator();
        return;
    }
    m_vboundsChunkSize = chunkSize;
    m_vboundsGridWidth = (m_heightStickWidth + chunkSize - 1) / chunkSize;
    m_vboundsGridLength = (m_heightStickLength + chunkSize - 1) / chunkSize;
    m_vboundsGrid.resize(m_vboundsGridWidth * m_vboundsGridLength);
    for (int cz = 0; cz < m_vboundsGridLength; ++cz) {
        int z0 = cz * chunkSize;
        for (int cx = 0; cx < m_vboundsGridWidth; ++cx) {
            int x0 = cx * chunkSize;
            Range r;
            // Heights are already BASE-relative (shifted by -minH in constructor)
            int sx = std::min(x0, m_heightStickWidth - 1), sz = std::min(z0, m_heightStickLength - 1);
            cbtScalar h0 = getRawHeightFieldValue(sx, sz) * m_heightScale;
            r.min = r.max = h0;
            for (int zz = z0; zz < z0 + chunkSize + 1 && zz < m_heightStickLength; ++zz) {
                for (int xx = x0; xx < x0 + chunkSize + 1 && xx < m_heightStickWidth; ++xx) {
                    cbtScalar h = getRawHeightFieldValue(xx, zz) * m_heightScale;
                    r.min = cbtMin(r.min, h);
                    r.max = cbtMax(r.max, h);
                }
            }
            m_vboundsGrid[cx + cz * m_vboundsGridWidth] = r;
        }
    }
}

void cbtHeightfieldChronoTerrainShape::clearAccelerator() {
    m_vboundsGrid.clear();
    m_vboundsGridWidth = 0;
    m_vboundsGridLength = 0;
    m_vboundsChunkSize = 0;
}

void cbtHeightfieldChronoTerrainShape::buildQuadExtents() {
    if (!m_useQuadExtentsCache)
        return;
    const int w = m_heightStickWidth - 1;
    const int l = m_heightStickLength - 1;
    m_quadExtents.resize(static_cast<std::size_t>(w) * l);

    for (int z = 0, idx = 0; z < l; ++z) {
        for (int x = 0; x < w; ++x, ++idx) {
            // Heights are already BASE-relative (shifted by -minH in constructor)
            const cbtScalar h00 = getRawHeightFieldValue(x, z) * m_heightScale;
            const cbtScalar h10 = getRawHeightFieldValue(x + 1, z) * m_heightScale;
            const cbtScalar h01 = getRawHeightFieldValue(x, z + 1) * m_heightScale;
            const cbtScalar h11 = getRawHeightFieldValue(x + 1, z + 1) * m_heightScale;

            QuadExtents e;
            e.minH = cbtMin(cbtMin(h00, h10), cbtMin(h01, h11));
            e.maxH = cbtMax(cbtMax(h00, h10), cbtMax(h01, h11));
            m_quadExtents[idx] = e;
        }
    }
}

// ---------------------------------------------------------------------------
// cbtHeightfieldChronoTerrainShape.cpp
// ---------------------------------------------------------------------------

void cbtHeightfieldChronoTerrainShape::buildVertexCache() {
    if (!m_useVertexCache)
        return;
    const int W = m_heightStickWidth;
    const int L = m_heightStickLength;
    m_vertexCache.resize(static_cast<std::size_t>(W) * L);

    // In unscaled space, grid step is always 1.0 (m_width = m_heightStickWidth - 1)

    for (int y = 0; y < L; ++y)
        for (int x = 0; x < W; ++x) {
            // Get height (BASE is at 0, apply heightScale once)
            cbtScalar h = getRawHeightFieldValue(x, y) * m_heightScale;

            // grid → local centered coords (step = 1.0 in unscaled space)
            cbtScalar lx = cbtScalar(x) - m_halfWidth;
            cbtScalar ly = cbtScalar(y) - m_halfLength;

            cbtVector3 v;
            switch (m_upAxis) {
                case 0:
                    v.setValue(h, lx, ly);
                    break;
                case 1:
                    v.setValue(lx, h, ly);
                    break;
                default:
                    v.setValue(lx, ly, h);
                    break;  // Z‑up
            }

        // Apply localScaling exactly once.
        v *= m_localScaling;
        m_vertexCache[static_cast<std::size_t>(y) * W + x] = v;  
    }
}


/// return centered+scaled height and normalized world‐gradient.
void cbtHeightfieldChronoTerrainShape::getHeightAndNormalAtGrid(
    const cbtScalar gridU,             // e.g. vertexTerrain[upAxis==1]? x : y
                              const cbtScalar gridV,             // e.g. z  or y
                              cbtScalar& outHeight,              // out: local‐centered, scaled
                              cbtVector3& outNormalLocal) const  // out: gradient in local coords
{
    // Reuse your existing queryHeightAndGradient, which expects (u,v)
    // in centered‐meters units and returns centered height + raw gradient.
    queryHeightAndGradient(gridU, gridV, outHeight, outNormalLocal);
    // Correct normal under non-uniform scaling: normals transform with inverse-transpose.
    // For axis-aligned scaling, this is equivalent to component-wise multiply by inverse scale.
    const cbtVector3 invS = getInverseLocalScaling();
    outNormalLocal = cbtVector3(outNormalLocal.x() * invS.x(), outNormalLocal.y() * invS.y(), outNormalLocal.z() * invS.z());
    if (outNormalLocal.length2() > SIMD_EPSILON)
        outNormalLocal.normalize();
}


void cbtHeightfieldChronoTerrainShape::processAllTriangles(cbtTriangleCallback* cb,
                                                           const cbtVector3& aabbMinWorld,
                                                           const cbtVector3& aabbMaxWorld) const {

    // Convert the world‑space AABB into the *un‑scaled* local grid frame
    const cbtVector3 invS = getInverseLocalScaling();   /// TODO:- should this be cached? with handling of cache - updates?
    const cbtVector3 aabbMin = aabbMinWorld * invS;
    const cbtVector3 aabbMax = aabbMaxWorld * invS;

    // Now proceed exactly as before ─ but all distances are in grid metres
    const cbtScalar dx = m_width / cbtScalar(m_heightStickWidth - 1);
    const cbtScalar dz = m_length / cbtScalar(m_heightStickLength - 1);
    const cbtScalar halfW = m_width * cbtScalar(0.5);
    const cbtScalar halfL = m_length * cbtScalar(0.5);

    int gx0, gx1, gz0, gz1;
    switch (m_upAxis) {
        case 0:
            gx0 = int(floor((aabbMin.y() + halfW) / dx));
            gx1 = int(ceil((aabbMax.y() + halfW) / dx));
            gz0 = int(floor((aabbMin.z() + halfL) / dz));
            gz1 = int(ceil((aabbMax.z() + halfL) / dz));
            break;
        case 1:
            gx0 = int(floor((aabbMin.x() + halfW) / dx));
            gx1 = int(ceil((aabbMax.x() + halfW) / dx));
            gz0 = int(floor((aabbMin.z() + halfL) / dz));
            gz1 = int(ceil((aabbMax.z() + halfL) / dz));
            break;
        default:
            gx0 = int(floor((aabbMin.x() + halfW) / dx));
            gx1 = int(ceil((aabbMax.x() + halfW) / dx));
            gz0 = int(floor((aabbMin.y() + halfL) / dz));
            gz1 = int(ceil((aabbMax.y() + halfL) / dz));
            break;
    }
    gx0 = cbtMax(0, cbtMin(gx0, m_heightStickWidth - 2));
    gx1 = cbtMax(0, cbtMin(gx1, m_heightStickWidth - 1));
    gz0 = cbtMax(0, cbtMin(gz0, m_heightStickLength - 2));
    gz1 = cbtMax(0, cbtMin(gz1, m_heightStickLength - 1));
    if (gx1 <= gx0 || gz1 <= gz0)
        return;

    // ----- step 2 : emit triangles (NO vertical filter) ----------------
    const int wQuads = m_heightStickWidth - 1;
    int triIndex = 0;

    auto emit = [&](cbtVector3 a, cbtVector3 b, cbtVector3 c) {
        // ensure up‑facing winding without re‑scaling
        cbtVector3 n = (b - a).cross(c - a);
        if (n[m_upAxis] < 0)
            std::swap(b, c);
        cbtVector3 tri[3] = {a, b, c};  // cached verts already include localScaling
        cb->processTriangle(tri, 0, triIndex++);
        
    };


    for (int z = gz0; z < gz1; ++z)
        for (int x = gx0; x < gx1; ++x) {
            cbtVector3 v00, v10, v01, v11;
            getVertex(x, z, v00);
            getVertex(x + 1, z, v10);
            getVertex(x, z + 1, v01);
            getVertex(x + 1, z + 1, v11);

            bool alt = m_flipQuadEdges || (m_useDiamondSubdivision && (((x + z) & 1) != 0)) ||
                       (m_useZigzagSubdivision && ((z & 1) != 0));

            if (alt) {
                emit(v00, v10, v11);
                emit(v00, v11, v01);
            } else {
                emit(v00, v10, v01);
                emit(v10, v11, v01);
            }
        }
}

// Extract planar (u,v) from local unscaled Pl based on upAxis
void cbtHeightfieldChronoTerrainShape::getUV(const cbtVector3& Pl, cbtScalar& u, cbtScalar& v) const {
    switch (getUpAxis()) {
        case 0:
            u = Pl.y();
            v = Pl.z();
            break;  // X-up
        case 1:
            u = Pl.x();
            v = Pl.z();
            break;  // Y-up
        default:
            u = Pl.x();
            v = Pl.y();
            break;  // Z-up
    }
}


// Simple bilinear height interpolation given integer cell indices and fractional offsets
void cbtHeightfieldChronoTerrainShape::getBilinearHeight(int cellX,
                                                         int cellZ,
                                                         cbtScalar fracX,
                                                         cbtScalar fracZ,
                                                         cbtScalar& outHeight) const {
    // Fetch the four corner heights (scaled by heightScale)
    auto heightAt = [&](int x, int z) { return getRawHeightFieldValue(x, z) * m_heightScale; };
    cbtScalar h00 = heightAt(cellX, cellZ);
    cbtScalar h10 = heightAt(cellX + 1, cellZ);
    cbtScalar h01 = heightAt(cellX, cellZ + 1);
    cbtScalar h11 = heightAt(cellX + 1, cellZ + 1);

    // Use the class's SIMD-accelerated bilinear interpolation
    outHeight = BilinearHeight(h00, h10, h01, h11, fracX, fracZ);
}


// Full height and gradient query at any (u,v) in the heightfield grid
void cbtHeightfieldChronoTerrainShape::queryHeightAndGradient(cbtScalar coordU,
                                                              cbtScalar coordV,
                                                              cbtScalar& outHeight,
                                                              cbtVector3& outNormalLocal) const {
    // In unscaled local space, grid step is always 1.0 (m_width = m_heightStickWidth - 1)
    // So mapping from centered coords to grid indices is a simple offset

    // Map (coordU,coordV) in centered grid units to [0..width-1]/[0..length-1]
    cbtScalar gridX = coordU + m_halfWidth;
    cbtScalar gridZ = coordV + m_halfLength;

    // Clamp to valid cell range
    cbtClamp(gridX, cbtScalar(0), cbtScalar(m_heightStickWidth - 1) - cbtScalar(1e-6));
    cbtClamp(gridZ, cbtScalar(0), cbtScalar(m_heightStickLength - 1) - cbtScalar(1e-6));

    int cellX = int(std::floor(gridX));
    int cellZ = int(std::floor(gridZ));
    cbtScalar fracX = gridX - cellX;
    cbtScalar fracZ = gridZ - cellZ;

    // Get the interpolated height (already BASE-relative from stored data)
    cbtScalar heightLocal;
    getBilinearHeight(cellX, cellZ, fracX, fracZ, heightLocal);

    // Heights are already BASE-relative (shifted by -minH in constructor)
    outHeight = heightLocal;

    // Compute local partial derivatives ∂H/∂u, ∂H/∂v in meters/meter
    // We can reuse the raw data directly:
    auto rawH = [&](int x, int z) { return getRawHeightFieldValue(x, z) * m_heightScale; };
    // Grid step in unscaled space is 1.0, so gradient = height difference directly
    cbtScalar dhdx = (rawH(cellX + 1, cellZ) - rawH(cellX, cellZ)) * (1 - fracZ) +
                     (rawH(cellX + 1, cellZ + 1) - rawH(cellX, cellZ + 1)) * fracZ;
    cbtScalar dhdz = (rawH(cellX, cellZ + 1) - rawH(cellX, cellZ)) * (1 - fracX) +
                     (rawH(cellX + 1, cellZ + 1) - rawH(cellX + 1, cellZ)) * fracX;

    // Build the local‐space normal (upAxis defines which component is “height”)
    switch (m_upAxis) {
        case 0:  // X‐up: height ∈ X, so normal = [ 1, -∂H/∂u, -∂H/∂v ]
            outNormalLocal.setValue(cbtScalar(1), -dhdx, -dhdz);
            break;
        case 1:  // Y‐up
            outNormalLocal.setValue(-dhdx, cbtScalar(1), -dhdz);
            break;
        default:  // Z‐up
            outNormalLocal.setValue(-dhdx, -dhdz, cbtScalar(1));
            break;
    }

    // Robust normalize (avoid NaNs on extreme values)
    auto finite3 = [](const cbtVector3& v) {
        return std::isfinite((double)v.x()) && std::isfinite((double)v.y()) && std::isfinite((double)v.z());
    };
    if (!finite3(outNormalLocal) || outNormalLocal.length2() < SIMD_EPSILON) {
        // Fallback to pure up-axis normal
        outNormalLocal.setValue(0, 0, 0);
        outNormalLocal[m_upAxis] = cbtScalar(1);
    } else {
        outNormalLocal.normalize();
    }
}

//  Grid‐raycast machinery in an unnamed namespace
// Essentially copies the Bullet latest approach here... thus far its not used
namespace {
struct GridRaycastState {
    int x, z, prev_x, prev_z;
    cbtScalar param, prevParam;
    cbtScalar maxDistanceFlat, maxDistance3d;
};

// Amanatides–Woo Bresenham‐style grid traversal
template <typename Action_T>
void gridRaycast(Action_T& action, const cbtVector3& from, const cbtVector3& to, int indices[3]) {
    GridRaycastState rs;
    rs.maxDistance3d = from.distance(to);
    if (rs.maxDistance3d < cbtScalar(1e-4))
        return;

    cbtScalar dx = to[indices[0]] - from[indices[0]];
    cbtScalar dz = to[indices[2]] - from[indices[2]];
    rs.maxDistanceFlat = cbtSqrt(dx * dx + dz * dz);

    cbtScalar dirX = (rs.maxDistanceFlat > cbtScalar(1e-4)) ? dx / rs.maxDistanceFlat : cbtScalar(0);
    cbtScalar dirZ = (rs.maxDistanceFlat > cbtScalar(1e-4)) ? dz / rs.maxDistanceFlat : cbtScalar(0);

    int stepX = (dirX > 0 ? 1 : (dirX < 0 ? -1 : 0));
    int stepZ = (dirZ > 0 ? 1 : (dirZ < 0 ? -1 : 0));

    const cbtScalar inf = cbtScalar(1e20);
    cbtScalar deltaX = (stepX != 0 ? cbtScalar(1) / cbtFabs(dirX) : inf);
    cbtScalar deltaZ = (stepZ != 0 ? cbtScalar(1) / cbtFabs(dirZ) : inf);

    cbtScalar crossX = (stepX > 0 ? (ceil(from[indices[0]]) - from[indices[0]]) * deltaX
                                  : (from[indices[0]] - floor(from[indices[0]])) * deltaX);
    cbtScalar crossZ = (stepZ > 0 ? (ceil(from[indices[2]]) - from[indices[2]]) * deltaZ
                                  : (from[indices[2]] - floor(from[indices[2]])) * deltaZ);

    rs.x = int(floor(from[indices[0]]));
    rs.z = int(floor(from[indices[2]]));
    rs.prev_x = rs.x;
    rs.prev_z = rs.z;
    rs.param = cbtScalar(0);

    while (true) {
        rs.prevParam = rs.param;
        if (crossX < crossZ) {
            rs.x += stepX;
            rs.param = crossX;
            crossX += deltaX;
        } else {
            rs.z += stepZ;
            rs.param = crossZ;
            crossZ += deltaZ;
        }
        if (rs.param > rs.maxDistanceFlat) {
            rs.param = rs.maxDistanceFlat;
            action(rs);
            break;
        } else {
            action(rs);
        }
    }
}
}  // unnamed namespace

// performRaycast: supports ANY upAxis, with winding-corrected triangles
// and hierarchical Bresenham + chunked culling
// This is a slightly modified version of the Bullet code
void cbtHeightfieldChronoTerrainShape::performRaycast(cbtTriangleCallback* callback,
                                                      const cbtVector3& raySource,
                                                      const cbtVector3& rayTarget) const {
    // Explicit AABB check and early exit
    cbtVector3 shapeAabbMin, shapeAabbMax;
    getAabb(cbtTransform::getIdentity(), shapeAabbMin, shapeAabbMax);
    cbtVector3 rayMin(cbtMin(raySource.getX(), rayTarget.getX()), cbtMin(raySource.getY(), rayTarget.getY()),
                      cbtMin(raySource.getZ(), rayTarget.getZ()));
    cbtVector3 rayMax(cbtMax(raySource.getX(), rayTarget.getX()), cbtMax(raySource.getY(), rayTarget.getY()),
                      cbtMax(raySource.getZ(), rayTarget.getZ()));
    if (rayMax.getX() < shapeAabbMin.getX() || rayMin.getX() > shapeAabbMax.getX() ||
        rayMax.getY() < shapeAabbMin.getY() || rayMin.getY() > shapeAabbMax.getY() ||
        rayMax.getZ() < shapeAabbMin.getZ() || rayMin.getZ() > shapeAabbMax.getZ())
        return;

    // Transform into cell-local space (m_localOrigin is always (0,0,0) for BASE-at-origin mode)
    cbtVector3 beginPos = raySource / getLocalScaling();
    cbtVector3 endPos = rayTarget / getLocalScaling();

    // Functor to emit the two triangles in each cell (with winding & subdivision)
    struct ProcessTrianglesAction {
        const cbtHeightfieldChronoTerrainShape* s;
        bool flipQuadEdges, useDiamondSubdivision, useZigzagSubdivision;
        int width, length, upAxis;
        cbtTriangleCallback* cb;
        void exec(int x, int z) const {
            if (x < 0 || z < 0 || x >= width || z >= length)
                return;
            cbtVector3 v00, v10, v01, v11;
            s->getVertex(x, z, v00);
            s->getVertex(x + 1, z, v10);
            s->getVertex(x, z + 1, v01);
            s->getVertex(x + 1, z + 1, v11);
            bool useAlt = flipQuadEdges || (useDiamondSubdivision && (((z + x) & 1) != 0)) ||
                          (useZigzagSubdivision && ((z & 1) != 0));
            cbtVector3 upVec(0, 0, 0);
            upVec[upAxis] = cbtScalar(1);

            if (useAlt) {
                cbtVector3 tri0[3] = {v00, v10, v11};
                if (((v10 - v00).cross(v11 - v00)).dot(upVec) < 0)
                    std::swap(tri0[1], tri0[2]);
                cb->processTriangle(tri0, x, z);
                cbtVector3 tri1[3] = {v00, v11, v01};
                if (((v11 - v00).cross(v01 - v00)).dot(upVec) < 0)
                    std::swap(tri1[1], tri1[2]);
                cb->processTriangle(tri1, x, z);
            } else {
                cbtVector3 tri0[3] = {v00, v10, v01};
                if (((v10 - v00).cross(v01 - v00)).dot(upVec) < 0)
                    std::swap(tri0[1], tri0[2]);
                cb->processTriangle(tri0, x, z);
                cbtVector3 tri1[3] = {v10, v11, v01};
                if (((v11 - v10).cross(v01 - v10)).dot(upVec) < 0)
                    std::swap(tri1[1], tri1[2]);
                cb->processTriangle(tri1, x, z);
            }
        }
        void operator()(const GridRaycastState& rs) const { exec(rs.prev_x, rs.prev_z); }
    } proc = {this,
              m_flipQuadEdges,
              m_useDiamondSubdivision,
              m_useZigzagSubdivision,
              m_heightStickWidth - 1,
              m_heightStickLength - 1,
              m_upAxis,
              callback};

    // upAxis to grid-axis mapping
    int indices[3] = {0, 1, 2};
    indices[0] = (m_upAxis == 0 ? 1 : 0);
    indices[2] = (m_upAxis == 2 ? 1 : 2);

    // Single-cell “vertical” ray optimisation
    int iBX = int(floor(beginPos[indices[0]]));
    int iBZ = int(floor(beginPos[indices[2]]));
    int iEX = int(floor(endPos[indices[0]]));
    int iEZ = int(floor(endPos[indices[2]]));
    if (iBX == iEX && iBZ == iEZ) {
        proc.exec(iBX, iBZ);
        return;
    }

    // Hierarchical Bresenham walk with the chunk culling
    if (m_vboundsGrid.size() == 0) {
        gridRaycast(proc, beginPos, endPos, indices);
    } else {
        cbtVector3 diff = endPos - beginPos;
        cbtScalar flat2 = diff[indices[0]] * diff[indices[0]] + diff[indices[2]] * diff[indices[2]];
        if (flat2 < m_vboundsChunkSize * m_vboundsChunkSize) {
            gridRaycast(proc, beginPos, endPos, indices);
        } else {
            struct ProcessVBoundsAction {
                const cbtAlignedObjectArray<Range>& vb;
                int width, length, chunkSize;
                cbtVector3 rayB, rayE, rayD;
                int idx[3];
                ProcessTrianglesAction proc;
                void operator()(const GridRaycastState& rs) const {
                    int gx = rs.prev_x, gz = rs.prev_z;
                    if (gx < 0 || gz < 0 || gx >= width || gz >= length)
                        return;
                    const Range& r = vb[gx + gz * width];
                    // compute entry/exit in 3D
                    cbtScalar scale = chunkSize * rs.maxDistance3d / rs.maxDistanceFlat;
                    cbtVector3 enterP = rayB + rayD * (rs.prevParam * scale);
                    cbtVector3 exitP = rayB + rayD * (rs.param * scale);
                    int h = idx[1];  // the “height” axis
                    if ((enterP[h] > r.max && exitP[h] > r.max) || (enterP[h] < r.min && exitP[h] < r.min))
                        return;
                    // pass &idx[0] to match gridRaycast(int indices[3])
                    gridRaycast(proc, enterP, exitP, const_cast<int*>(idx));
                }
            } fb = {m_vboundsGrid,
                    m_vboundsGridWidth,
                    m_vboundsGridLength,
                    m_vboundsChunkSize,
                    beginPos,
                    endPos,
                    (endPos - beginPos).normalized(),
                    {indices[0], indices[1], indices[2]},
                    proc};
            gridRaycast(fb, beginPos / m_vboundsChunkSize, endPos / m_vboundsChunkSize, fb.idx);
        }
    }
}

// sampleWorld with minimised transform operations
bool cbtHeightfieldChronoTerrainShape::sampleWorld(const cbtTransform& terrainFrame,
                                                   const cbtVector3& queryPointWorld,
                                                   cbtVector3& outSurfacePointWorld,
                                                   cbtVector3& outSurfaceNormalWorld) const {
    // Transform to local unscaled space
    const cbtVector3 invScale = getInverseLocalScaling();
    const cbtVector3 queryLocal = terrainFrame.invXform(queryPointWorld) * invScale;

    // Extract planar coordinates
    cbtScalar u, v;
    switch (m_upAxis) {
        case 0:
            u = queryLocal.y();
            v = queryLocal.z();
            break;
        case 1:
            u = queryLocal.x();
            v = queryLocal.z();
            break;
        default:
            u = queryLocal.x();
            v = queryLocal.y();
            break;
    }

    // Bounds check
    if (u < -m_halfWidth || u > m_halfWidth || v < -m_halfLength || v > m_halfLength)
        return false;

    // Query height and gradient (height is relative to BASE, i.e., local z=0)
    cbtScalar terrainHeight_fromBase;
    cbtVector3 gradientLocal;
    queryHeightAndGradient(u, v, terrainHeight_fromBase, gradientLocal);

    // Build surface point in local space (BASE at origin)
    cbtVector3 surfacePointLocal;
    switch (m_upAxis) {
        case 0:
            surfacePointLocal.setValue(terrainHeight_fromBase, u, v);
            break;
        case 1:
            surfacePointLocal.setValue(u, terrainHeight_fromBase, v);
            break;
        default:
            surfacePointLocal.setValue(u, v, terrainHeight_fromBase);
            break;
    }

    // Apply scaling
    surfacePointLocal *= getLocalScaling();

    // Transform to world
    outSurfacePointWorld = terrainFrame * surfacePointLocal;

    // Correct normal under non-uniform scaling.
    // queryHeightAndGradient returns a unit normal in *unscaled* local coordinates.
    // With axis-aligned scaling S, the normal transforms as n' = normalize(S^{-T} n).
    const cbtVector3 invS = getInverseLocalScaling();
    cbtVector3 normalLocalScaled(gradientLocal.x() * invS.x(), gradientLocal.y() * invS.y(), gradientLocal.z() * invS.z());
    if (normalLocalScaled.length2() < SIMD_EPSILON)
        return false;
    normalLocalScaled.normalize();
    outSurfaceNormalWorld = terrainFrame.getBasis() * normalLocalScaled;
    outSurfaceNormalWorld.normalize();

    return true;
}
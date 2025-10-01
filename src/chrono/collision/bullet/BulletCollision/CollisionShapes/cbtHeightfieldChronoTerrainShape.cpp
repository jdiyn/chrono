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
// which is inspired by upon the bullet version but aims to improve it
// =============================================================================

#include "cbtHeightfieldChronoTerrainShape.h"
#include "chrono/collision/bullet/LinearMath/cbtTransformUtil.h"
#include <algorithm>
#include <cmath>
#include <iostream>

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
    // Convert query coordinates (u,v) into integer grid cell indices
    int iu, iv;
    cbtScalar fu, fv;

    switch (m_upAxis) {
        case 0:  // X-up: (u,v) = (y,z)
            iu = std::clamp(int(std::floor(u)), 0, m_heightStickWidth - 2);
            iv = std::clamp(int(std::floor(v)), 0, m_heightStickLength - 2);
            fu = u - iu;
            fv = v - iv;
            break;
        case 1:  // Y-up: (u,v) = (x,z)
            iu = std::clamp(int(std::floor(u)), 0, m_heightStickWidth - 2);
            iv = std::clamp(int(std::floor(v)), 0, m_heightStickLength - 2);
            fu = u - iu;
            fv = v - iv;
            break;
        default:  // Z-up: (u,v) = (x,y)
            iu = std::clamp(int(std::floor(u)), 0, m_heightStickWidth - 2);
            iv = std::clamp(int(std::floor(v)), 0, m_heightStickLength - 2);
            fu = u - iu;
            fv = v - iv;
            break;
    }

    // Heights at grid corners (already scaled)
    auto H = [this](int x, int y) { return getRawHeightFieldValue(x, y) * m_heightScale; };

    const cbtScalar h00 = H(iu, iv);
    const cbtScalar h10 = H(iu + 1, iv);
    const cbtScalar h01 = H(iu, iv + 1);
    const cbtScalar h11 = H(iu + 1, iv + 1);

    // Bilinear interpolation
    const cbtScalar h0 = h00 + fu * (h10 - h00);
    const cbtScalar h1 = h01 + fu * (h11 - h01);
    outH = h0 + fv * (h1 - h0);

    // Compute gradients in u,v
    cbtScalar du = ((h10 - h00) * (1 - fv) + (h11 - h01) * fv);
    cbtScalar dv = ((h01 - h00) * (1 - fu) + (h11 - h10) * fu);

    // Convert gradients to proper world axes
    switch (m_upAxis) {
        case 0:
            outGrad.setValue(1.0, -du, -dv);
            break;  // X-up
        case 1:
            outGrad.setValue(-du, 1.0, -dv);
            break;  // Y-up
        default:
            outGrad.setValue(-du, -dv, 1.0);
            break;  // Z-up
    }
    outGrad.normalize();
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
    vtx = m_vertexCache[static_cast<std::size_t>(y) * m_heightStickWidth + x];
}

void cbtHeightfieldChronoTerrainShape::setLocalScaling(const cbtVector3& s) {
    // store locally – this fulfils the pure‑virtual contract
    m_localScaling = s;

    // rebuild lookup tables that depend on scale
    buildVertexCache();
    buildQuadExtents();
    if (m_vboundsChunkSize > 0)
        buildAccelerator(m_vboundsChunkSize);
}

// raw height from user array (uncentered; heightScale & centering applied elsewhere)
inline cbtScalar cbtHeightfieldChronoTerrainShape::getRawHeightFieldValue(int x, int y) const {
    return m_heightfieldData[y * m_heightStickWidth + x];
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
    cbtAssert(heightStickWidth > 1 && heightStickLength > 1 && minH <= maxH && up >= 0 && up < 3);

    m_width = cbtScalar(heightStickWidth - 1);
    m_length = cbtScalar(heightStickLength - 1);

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


    // build caching first time
    buildQuadExtents();
    buildVertexCache();
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

// Untested - dynamic height update
void cbtHeightfieldChronoTerrainShape::updateHeight(int x, int y, cbtScalar newHeight_absolute) {
    cbtAssert(x >= 0 && x < m_heightStickWidth && y >= 0 && y < m_heightStickLength);

    // Convert absolute height to base-relative
    m_heightfieldData[y * m_heightStickWidth + x] = newHeight_absolute - m_minHeight;

    // Update cached data
    buildVertexCache();  // Rebuild affected regions
    if (m_vboundsChunkSize > 0) {
        // Optionally: rebuild only affected chunks
        buildAccelerator(m_vboundsChunkSize);
    }
}

void cbtHeightfieldChronoTerrainShape::updateHeights(const cbtScalar* newHeights_absolute, int numSamples) {
    cbtAssert(numSamples == m_heightStickWidth * m_heightStickLength);

    for (int i = 0; i < numSamples; ++i) {
        m_heightfieldData[i] = newHeights_absolute[i] - m_minHeight;
    }

    buildVertexCache();
    if (m_vboundsChunkSize > 0) {
        buildAccelerator(m_vboundsChunkSize);
    }
}

// accelerator build/clear
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
            // init
            int sx = std::min(x0, m_heightStickWidth - 1), sz = std::min(z0, m_heightStickLength - 1);
            cbtScalar h0 = getRawHeightFieldValue(sx, sz) * m_heightScale - m_localOrigin[m_upAxis];
            r.min = r.max = h0;
            for (int zz = z0; zz < z0 + chunkSize + 1 && zz < m_heightStickLength; ++zz) {
                for (int xx = x0; xx < x0 + chunkSize + 1 && xx < m_heightStickWidth; ++xx) {
                    cbtScalar h = getRawHeightFieldValue(xx, zz) * m_heightScale - m_localOrigin[m_upAxis];
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
    const int w = m_heightStickWidth - 1;
    const int l = m_heightStickLength - 1;
    m_quadExtents.resize(static_cast<std::size_t>(w) * l);

    for (int z = 0, idx = 0; z < l; ++z) {
        for (int x = 0; x < w; ++x, ++idx) {
            const cbtScalar h00 = getRawHeightFieldValue(x, z) * m_heightScale - m_localOrigin[m_upAxis];
            const cbtScalar h10 = getRawHeightFieldValue(x + 1, z) * m_heightScale - m_localOrigin[m_upAxis];
            const cbtScalar h01 = getRawHeightFieldValue(x, z + 1) * m_heightScale - m_localOrigin[m_upAxis];
            const cbtScalar h11 = getRawHeightFieldValue(x + 1, z + 1) * m_heightScale - m_localOrigin[m_upAxis];

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
    const int W = m_heightStickWidth;
    const int L = m_heightStickLength;
    m_vertexCache.resize(static_cast<std::size_t>(W) * L);

    const cbtVector3 S = m_localScaling;  // cache once

    const cbtScalar dx = m_width / cbtScalar(W - 1);
    const cbtScalar dy = m_length / cbtScalar(L - 1);
    const cbtScalar halfW = m_width * cbtScalar(0.5);
    const cbtScalar halfL = m_length * cbtScalar(0.5);

    for (int y = 0; y < L; ++y)
        for (int x = 0; x < W; ++x) {
            // Get height (already offset so BASE is at 0)
            cbtScalar h = getRawHeightFieldValue(x, y) * m_heightScale; // do we even need to scale again here??

            // grid → local before scaling
            cbtScalar lx = x * dx - halfW;
            cbtScalar ly = y * dy - halfL;

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

        // **apply scaling exactly once**
        v *= S; // <<---- we're scaling again here, is that correct?
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
    // Then scale the gradient into world by shape‐scaling:
    outNormalLocal = outNormalLocal * getLocalScaling();
    outNormalLocal.normalize();
}


void cbtHeightfieldChronoTerrainShape::processAllTriangles(cbtTriangleCallback* cb,
                                                           const cbtVector3& aabbMinWorld,
                                                           const cbtVector3& aabbMaxWorld) const {

    // Convert the world‑space AABB into the *un‑scaled* local grid frame
    const cbtVector3 invS = getInverseLocalScaling();
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

    // Bilinear interpolation
    cbtScalar h0 = h00 + fracX * (h10 - h00);
    cbtScalar h1 = h01 + fracX * (h11 - h01);
    outHeight = h0 + fracZ * (h1 - h0);
}


// Full height and gradient query at any (u,v) in the heightfield grid
void cbtHeightfieldChronoTerrainShape::queryHeightAndGradient(cbtScalar coordU,
                                                              cbtScalar coordV,
                                                              cbtScalar& outHeight,
                                                              cbtVector3& outNormalLocal) const {
    // Grid spacing and half‑extents
    cbtScalar dx = m_width / (m_heightStickWidth - 1);
    cbtScalar dz = m_length / (m_heightStickLength - 1);
    cbtScalar halfW = m_width * cbtScalar(0.5);
    cbtScalar halfL = m_length * cbtScalar(0.5);

    // Map (coordU,coordV) in centered meters to [0..width-1]/[0..length-1]
    cbtScalar gridX = (coordU + halfW) / dx;
    cbtScalar gridZ = (coordV + halfL) / dz;

    // Clamp to valid cell range
    cbtClamp(gridX, cbtScalar(0), cbtScalar(m_heightStickWidth - 1) - cbtScalar(1e-6));
    cbtClamp(gridZ, cbtScalar(0), cbtScalar(m_heightStickLength - 1) - cbtScalar(1e-6));

    int cellX = int(std::floor(gridX));
    int cellZ = int(std::floor(gridZ));
    cbtScalar fracX = gridX - cellX;
    cbtScalar fracZ = gridZ - cellZ;

    // Get the interpolated height (still centered, un‑scaled)
    cbtScalar heightCentered;
    getBilinearHeight(cellX, cellZ, fracX, fracZ, heightCentered);

    // Center around the stored origin (same as getVertex)
    outHeight = heightCentered - m_localOrigin[m_upAxis];

    // Compute local partial derivatives ∂H/∂u, ∂H/∂v in meters/meter
    // We can reuse the raw data directly:
    auto rawH = [&](int x, int z) { return getRawHeightFieldValue(x, z) * m_heightScale; };
    cbtScalar dhdx = ((rawH(cellX + 1, cellZ) - rawH(cellX, cellZ)) * (1 - fracZ) +
                      (rawH(cellX + 1, cellZ + 1) - rawH(cellX, cellZ + 1)) * fracZ) /
                     dx;
    cbtScalar dhdz = ((rawH(cellX, cellZ + 1) - rawH(cellX, cellZ)) * (1 - fracX) +
                      (rawH(cellX + 1, cellZ + 1) - rawH(cellX + 1, cellZ)) * fracX) /
                     dz;

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
    outNormalLocal.normalize();
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

    // Transform into cell-local space
    cbtVector3 beginPos = raySource / getLocalScaling() + m_localOrigin;
    cbtVector3 endPos = rayTarget / getLocalScaling() + m_localOrigin;

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
    const cbtScalar halfW = m_width * cbtScalar(0.5);
    const cbtScalar halfL = m_length * cbtScalar(0.5);
    if (u < -halfW || u > halfW || v < -halfL || v > halfL)
        return false;

    // Query height and gradient (height is relative to BASE, i.e., local z=0)
    cbtScalar terrainHeight_fromBase;
    cbtVector3 gradientLocal;
    queryHeightAndGradient(u, v, terrainHeight_fromBase, gradientLocal);

    if (gradientLocal.length2() < SIMD_EPSILON)
        return false;

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

    cbtVector3 scaledGradient = gradientLocal * getLocalScaling();
    outSurfaceNormalWorld = (terrainFrame.getBasis() * scaledGradient).normalized();

    return true;
}
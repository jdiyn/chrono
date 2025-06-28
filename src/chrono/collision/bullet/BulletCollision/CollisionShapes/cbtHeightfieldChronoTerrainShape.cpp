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
// which draws upon the bullet version
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



//  Heightfield fast bilinear sampler - using cells
void cbtHeightfieldChronoTerrainShape::sampleHeight(cbtScalar xu,
                                                    cbtScalar zv,
                                                    cbtScalar& outH,
                                                    cbtVector3& outGrad) const {
    // locate the cell in the heightfield grid
    int ix = int(std::floor(xu));
    int iz = int(std::floor(zv));
    ix = std::clamp(ix, 0, m_heightStickWidth - 2);
    iz = std::clamp(iz, 0, m_heightStickLength - 2);

    // fractional location inside that cell
    cbtScalar fu = xu - ix;
    cbtScalar fv = zv - iz;

    // fetch the four corner heights (raw, un-scaled)
    cbtScalar h00 = getRawHeightFieldValue(ix, iz);
    cbtScalar h10 = getRawHeightFieldValue(ix + 1, iz);
    cbtScalar h01 = getRawHeightFieldValue(ix, iz + 1);
    cbtScalar h11 = getRawHeightFieldValue(ix + 1, iz + 1);

    // bilinear interpolate height
    cbtScalar h0 = h00 + fu * (h10 - h00);
    cbtScalar h1 = h01 + fu * (h11 - h01);
    cbtScalar h = h0 + fv * (h1 - h0);

    // compute the partials (du = ∂H/∂u, dv = ∂H/∂v)
    cbtScalar du = (h10 - h00) * (1 - fv) + (h11 - h01) * fv;
    cbtScalar dv = (h01 - h00) * (1 - fu) + (h11 - h10) * fu;

    // apply the heightScale
    h *= m_heightScale;
    du *= m_heightScale;
    dv *= m_heightScale;

    // pack into a gradient vector in the correct axes
    cbtVector3 grad(0, 0, 0);
    switch (m_upAxis) {
        case 0:  // X-up
            grad.setValue(0, du, dv);
            break;
        case 1:  // Y-up
            grad.setValue(du, 0, dv);
            break;
        default:  // Z-up
            grad.setValue(du, dv, 0);
            break;
    }
    // send back the results
    outH = h;
    outGrad = grad;
}



// quantisation helper
static inline int getQuantized(cbtScalar x) {
    return (x < cbtScalar(0)) ? int(x - cbtScalar(0.5)) : int(x + cbtScalar(0.5));
}

void cbtHeightfieldChronoTerrainShape::quantizeWithClamp(int out[3], const cbtVector3& pt) const {
    cbtVector3 clamped = pt;
    clamped.setMax(m_localAabbMin); // ensure correct min / maxing here. max the aabb min and minimise the aabb max.
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
                                                                   const cbtScalar* data,
                                                                   cbtScalar scale,
                                                                   cbtScalar minH,
                                                                   cbtScalar maxH,
                                                                   int up,
                                                                   bool flip)
    : m_heightStickWidth(heightStickWidth),
      m_heightStickLength(heightStickLength),
      m_heightScale(scale),
      m_heightfieldData(data),
      m_minHeight(minH),
      m_maxHeight(maxH),
      m_upAxis(up),
      m_flipQuadEdges(flip),
      m_useDiamondSubdivision(false),
      m_useZigzagSubdivision(false),
      m_vboundsChunkSize(0),
      m_vboundsGridWidth(0),
      m_vboundsGridLength(0) {
    // initialise member variables
    m_shapeType = TERRAIN_SHAPE_PROXYTYPE;
    cbtAssert(heightStickWidth > 1 && heightStickLength > 1 && minH <= maxH && up >= 0 && up < 3);
    m_width = cbtScalar(heightStickWidth - 1);
    m_length = cbtScalar(heightStickLength - 1);
    // build unscaled local AABB & center
    switch (up) {
        case 0:
            m_localAabbMin.setValue(minH, 0, 0);
            m_localAabbMax.setValue(maxH, m_width, m_length);
            break;
        case 1:
            m_localAabbMin.setValue(0, minH, 0);
            m_localAabbMax.setValue(m_width, maxH, m_length);
            break;
        case 2:
            m_localAabbMin.setValue(0, 0, minH);
            m_localAabbMax.setValue(m_width, m_length, maxH);
            break;
    }
    m_localOrigin = (m_localAabbMin + m_localAabbMax) * cbtScalar(0.5);
    // build caching first time
    buildQuadExtents();
    buildVertexCache();
}


cbtHeightfieldChronoTerrainShape::~cbtHeightfieldChronoTerrainShape() {
    clearAccelerator();
}

// world‐space AABB
void cbtHeightfieldChronoTerrainShape::getAabb(const cbtTransform& t, cbtVector3& aabbMin, cbtVector3& aabbMax) const {
    // half‐extents in unscaled local
    cbtVector3 half = (m_localAabbMax - m_localAabbMin) * (getLocalScaling() * cbtScalar(0.5));
    // local center offset
    cbtVector3 localCenter(0, 0, 0);
    localCenter[m_upAxis] = (m_minHeight + m_maxHeight) * cbtScalar(0.5);
    localCenter *= getLocalScaling();
    // rotate half‐extents
    cbtMatrix3x3 abs_b = t.getBasis().absolute();
    cbtVector3 worldExtents(abs_b[0].dot(half), abs_b[1].dot(half), abs_b[2].dot(half));
    // add margin
    worldExtents += cbtVector3(getMargin(), getMargin(), getMargin());
    // center
    cbtVector3 worldCenter = t(localCenter);
    aabbMin = worldCenter - worldExtents;
    aabbMax = worldCenter + worldExtents;
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
    const int w = m_heightStickWidth - 1;
    const int l = m_heightStickLength - 1;
    m_quadExtents.resize(static_cast<std::size_t>(w) * l);

    for (int z = 0, idx = 0; z < l; ++z) {
        for (int x = 0; x < w; ++x, ++idx) {
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

void cbtHeightfieldChronoTerrainShape::buildVertexCache() {
    const int W = m_heightStickWidth;
    const int L = m_heightStickLength;
    m_vertexCache.resize(static_cast<std::size_t>(W) * L);

    for (int y = 0; y < L; ++y)
        for (int x = 0; x < W; ++x) {
            cbtScalar h = getRawHeightFieldValue(x, y) * m_heightScale;
            cbtVector3 v;

            // original getVertex math, but executed only once to build cache
            switch (m_upAxis) {
                case 0:
                    v.setValue(h - m_localOrigin.getX(), (-m_width * cbtScalar(0.5)) + x,
                               (-m_length * cbtScalar(0.5)) + y);
                    break;
                case 1:
                    v.setValue((-m_width * cbtScalar(0.5)) + x, h - m_localOrigin.getY(),
                               (-m_length * cbtScalar(0.5)) + y);
                    break;
                case 2:
                    v.setValue((-m_width * cbtScalar(0.5)) + x, (-m_length * cbtScalar(0.5)) + y,
                               h - m_localOrigin.getZ());
                    break;
            }
            v *= getLocalScaling();
            m_vertexCache[static_cast<std::size_t>(y) * W + x] = v;
        }
}

// Not needed with new convexheightfield collision algo
// TODO: fix inheritance to remove this method
void cbtHeightfieldChronoTerrainShape::processAllTriangles(cbtTriangleCallback* callback,
                                                           const cbtVector3& aabbMin,  // local-space AABB min?
                                                           const cbtVector3& aabbMax) const {

    // do nothing for now
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

// This is duplicated in the Chrono Shape.
bool cbtHeightfieldChronoTerrainShape::sampleWorld(const cbtTransform& terrainFrame,
                                                   const cbtVector3& QueryPoint,
                                                   cbtVector3& SurfacePoint,
                                                   cbtVector3& SurfaceNormal) const {
    // world to local coords
    const cbtVector3 Pl = terrainFrame.invXform(QueryPoint);  // to shape local
    const cbtVector3 invScale = getInverseLocalScaling();
    const cbtVector3 Pcell = Pl * invScale;  // HF‑local coords

    // planar (u,v) according to up‑axis
    cbtScalar pu, pv;
    switch (m_upAxis) {
        case 0:
            pu = Pcell.getY();
            pv = Pcell.getZ();
            break;  // X‑up
        case 1:
            pu = Pcell.getX();
            pv = Pcell.getZ();
            break;  // Y‑up
        default:
            pu = Pcell.getX();
            pv = Pcell.getY();
            break;  // Z‑up
    }

    // reject outside rectangle (half extents)
    if (cbtFabs(pu) > m_width * 0.5 || cbtFabs(pv) > m_length * 0.5)
        return false;

    // locate the grid cell and bilinear weights
    const cbtScalar cellU = (pu + m_width * 0.5) * (m_heightStickWidth - 1) / m_width;
    const cbtScalar cellV = (pv + m_length * 0.5) * (m_heightStickLength - 1) / m_length;

    int ix = int(cellU);
    int iz = int(cellV);
    const cbtScalar fu = cellU - ix;
    const cbtScalar fv = cellV - iz;

    // clamp on right / top border
    if (ix >= m_heightStickWidth - 1) {
        ix = m_heightStickWidth - 2;
    }
    if (iz >= m_heightStickLength - 1) {
        iz = m_heightStickLength - 2;
    }

    // corner heights (already × heightScale)
    auto H = [&](int x, int z) { return getRawHeightFieldValue(x, z) * m_heightScale; };
    const cbtScalar h00 = H(ix, iz);
    const cbtScalar h10 = H(ix + 1, iz);
    const cbtScalar h01 = H(ix, iz + 1);
    const cbtScalar h11 = H(ix + 1, iz + 1);

    const cbtScalar h0 = h00 + fu * (h10 - h00);
    const cbtScalar h1 = h01 + fu * (h11 - h01);
    const cbtScalar h = h0 + fv * (h1 - h0);

    // surface normal (central diff, same as Chrono)
    const cbtScalar du = (h10 - h00 + h11 - h01) * 0.5 * (m_heightStickWidth - 1) / m_width;
    const cbtScalar dv = (h01 - h00 + h11 - h10) * 0.5 * (m_heightStickLength - 1) / m_length;

    cbtVector3 nL;
    switch (m_upAxis) {
        case 0:
            nL.setValue(1, -du, -dv);
            break;  // X‑up
        case 1:
            nL.setValue(-du, 1, -dv);
            break;  // Y‑up
        default:
            nL.setValue(-du, -dv, 1);
            break;  // Z‑up
    }

    // scale-correct the gradient BEFORE normalising (same as convex path)
    const cbtVector3& S = getLocalScaling();  // (sx,sy,sz)
    nL.setValue(nL.x() * S.x(), nL.y() * S.y(), nL.z() * S.z());

    nL.normalize();
    SurfaceNormal = terrainFrame.getBasis() * nL;  // to world

    // world‑space surface point
    cbtVector3 PsL(0, 0, 0);
    switch (m_upAxis) {
        case 0:
            PsL.setValue(h - m_localOrigin.getX(), Pcell.getY(), Pcell.getZ());
            break;
        case 1:
            PsL.setValue(Pcell.getX(), h - m_localOrigin.getY(), Pcell.getZ());
            break;
        default:
            PsL.setValue(Pcell.getX(), Pcell.getY(), h - m_localOrigin.getZ());
            break;
    }
    PsL *= getLocalScaling();  // apply localScaling
    SurfacePoint = terrainFrame * PsL;              // to world

    return true;  // if requested, true if the point is on the terrain
}

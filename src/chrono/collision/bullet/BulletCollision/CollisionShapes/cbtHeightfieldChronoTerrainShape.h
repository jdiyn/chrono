// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
// =============================================================================
// Authors: Josh Diyn
// =============================================================================

#ifndef CBT_HEIGHTFIELD_CHRONO_TERRAIN_SHAPE_H
#define CBT_HEIGHTFIELD_CHRONO_TERRAIN_SHAPE_H

#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtConcaveShape.h"
#include "chrono/collision/bullet/LinearMath/cbtAlignedObjectArray.h"
#include "chrono/collision/bullet/LinearMath/cbtVector3.h"
#include "chrono/collision/bullet/LinearMath/cbtTransform.h"
#include "chrono/collision/bullet/LinearMath/cbtScalar.h"
#include <vector>

ATTRIBUTE_ALIGNED16(class)
cbtHeightfieldChronoTerrainShape : public cbtConcaveShape {
  public:
    BT_DECLARE_ALIGNED_ALLOCATOR();

    /// Simple min/max range used by the chunked accelerator
    struct Range {
        cbtScalar min, max;
        Range() {}
        Range(cbtScalar _min, cbtScalar _max) : min(_min), max(_max) {}
        bool overlaps(const Range& o) const { return !(min > o.max || max < o.min); }
    };

    void sampleHeight(cbtScalar xu, cbtScalar zv, cbtScalar & outH, cbtVector3 & outGrad) const;
    // (optionally – tiny helpers for width/length used by the algo)
    int getWidth() const {
        return m_heightStickWidth;
    }
    int getLength() const {
        return m_heightStickLength;
    }
    const cbtVector3& getLocalOrigin() const {
        return m_localOrigin;
    }

    /// constructor generic with cbtscalar - shoul dhandle either doubles or float
    cbtHeightfieldChronoTerrainShape(int heightStickWidth, int heightStickLength,
                                     const cbtScalar* heightfieldData,  // length = width*length
                                     cbtScalar heightScale, cbtScalar minHeight, cbtScalar maxHeight,
                                     int upAxis,  // 0=X,1=Y,2=Z
                                     bool flipQuadEdges);


    virtual ~cbtHeightfieldChronoTerrainShape();

    /// Diamond/zigzag subdivision control (Bullet 3.2+)
    inline void setUseDiamondSubdivision(bool d) {
        m_useDiamondSubdivision = d;
    }
    inline void setUseZigzagSubdivision(bool z) {
        m_useZigzagSubdivision = z;
    }

    // hook into the bullet local scaling to ensure cache is rebuilt
    void setLocalScaling(const cbtVector3& scaling) override;
    const cbtVector3& getLocalScaling() const override {
        return m_localScaling;
    }

    inline cbtVector3 getInverseLocalScaling() const {
        cbtScalar x = m_localScaling.getX();
        cbtScalar y = m_localScaling.getY();
        cbtScalar z = m_localScaling.getZ();

        // if any axis is zero, return identity to avoid divid by zero
        if (x == cbtScalar(0) || y == cbtScalar(0) || z == cbtScalar(0)) {
            return cbtVector3(cbtScalar(1), cbtScalar(1), cbtScalar(1));
        }

        return cbtVector3(cbtScalar(1) / x, cbtScalar(1) / y, cbtScalar(1) / z);
    }


    /// Build or clear the chunked min/max accelerator
    void buildAccelerator(int chunkSize);
    void clearAccelerator();

    // in cbtHeightfieldChronoTerrainShape.h (public)
    const Range& GetVBoundsChunk(int cx, int cz) const {
        return m_vboundsGrid[cx + cz * m_vboundsGridWidth];
    }



    //------------------------------------------------------------------------
    // cbtCollisionShape interface
    //------------------------------------------------------------------------

    /// World-space AABB (accounts for localScaling & rotation)
    virtual void getAabb(const cbtTransform& t, cbtVector3& aabbMin, cbtVector3& aabbMax) const override;

    /// Narrow-phase: triangles overlapping [aabbMin,aabbMax]
    virtual void processAllTriangles(cbtTriangleCallback * callback, const cbtVector3& aabbMin,
                                     const cbtVector3& aabbMax) const override;

    /// Hierarchical Bresenham raycast (Bullet 3.2+)
    void performRaycast(cbtTriangleCallback * callback, const cbtVector3& raySource, const cbtVector3& rayTarget) const;

    /// Concave shapes have zero local inertia
    virtual void calculateLocalInertia(cbtScalar /*mass*/, cbtVector3 & inertia) const override {
        inertia.setValue(0, 0, 0);
    }

    virtual const char* getName() const override {
        return "HEIGHTFIELD_CHRONO";
    }

    /// Which axis is “up” (0=X,1=Y,2=Z)
    int getUpAxis() const {
        return m_upAxis;
    }

    /// Scaled height at integer grid coordinates (x,z)
    cbtScalar getHeight(int x, int z) const;


    /// Sample height‑field at an arbitrary world point.
    /// Return false when (Pw) projects outside the X/Z (or Y/Z etc.) extent.
    bool sampleWorld(const cbtTransform& terrainFrame,  // world‑space frame of the shape
                     const cbtVector3& QueryPoint,              // query point in world
                     cbtVector3& SurfacePoint,                    // surface point (world)
                     cbtVector3& SurfaceNormal) const;         // surface normal (world)


  protected:
    // raw height data (must outlive this shape)
    const cbtScalar* m_heightfieldData;

    int m_heightStickWidth; // number of height samples
    int m_heightStickLength; // number of height samples
    cbtScalar m_heightScale;
    cbtScalar m_minHeight;
    cbtScalar m_maxHeight;
    cbtScalar m_width;   // = heightStickWidth  – 1
    cbtScalar m_length;  // = heightStickLength – 1
    int m_upAxis;
    bool m_flipQuadEdges;
    bool m_useDiamondSubdivision;
    bool m_useZigzagSubdivision;

    // AABB in local (unscaled) coords
    cbtVector3 m_localAabbMin;
    cbtVector3 m_localAabbMax;
    cbtVector3 m_localOrigin;  // = 0.5*(min+max)

    // local scaling as just 1,1,1 for initialisation
    cbtVector3 m_localScaling{1, 1, 1};

    // chunked min/max grid for raycast accelerator
    cbtAlignedObjectArray<Range> m_vboundsGrid;
    int m_vboundsGridWidth;
    int m_vboundsGridLength;
    int m_vboundsChunkSize;

    // caching quads
    struct QuadExtents {
        cbtScalar minH, maxH;
    };
    std::vector<QuadExtents> m_quadExtents;
    void buildQuadExtents();  // call whenever the heightfield changes

    // caching vertices
    std::vector<cbtVector3> m_vertexCache;  // (width  × length) grid
    void buildVertexCache();                // rebuild when data or scaling changes! user need to manage
    // TODO:: set a public function rebuildcache to build quad extents and verteces
    //// replicate Bullet’s getVertex (including centering by m_localOrigin)
    void getVertex(int x, int y, cbtVector3& vertex) const;

    // helper routines (defined in .cpp)
    void quantizeWithClamp(int out[3], const cbtVector3& pt) const;
    cbtScalar getRawHeightFieldValue(int x, int y) const;
    static Range minmaxRange(cbtScalar a, cbtScalar b, cbtScalar c);
};

#endif  // CBT_HEIGHTFIELD_CHRONO_TERRAIN_SHAPE_H

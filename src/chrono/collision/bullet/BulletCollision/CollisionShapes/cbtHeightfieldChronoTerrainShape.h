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
//
// Optimised heightfield terrain shape for Chrono physics
//
// ## Intended for high performance analytical collision response with large terrains
// - Sphere collision: O(1) analytical with bilinear interpolation + ClosestPointOnTriangle
// - Convex collision: O(1) height rejection + support sampling (26-42 directions)
// - Ray queries: O(1) bilinear interpolation with gradient-based normals
//
// ## Coordinate Conventions (BASE-at-origin standard)
// - **Height storage**: Heights in m_heightfieldData are BASE-relative (shifted by -minH
//   in constructor). This means local height=0 corresponds to minHeight.
// - **Planar origin**: Center of the heightfield patch (planar dimensions centered at 0)
// - **Height range**: 0 to (maxHeight - minHeight) in local coordinates
// - **Grid indexing**: Row-major [j * width + i] where i is along width axis
// - **Up-axis**: 0=X-up, 1=Y-up, 2=Z-up (default Z)
// - **Local scaling**: Applied after height lookup via m_localScaling vector
//
// Terrain Base/origin/ref sits at local z=0
//
// ## Caching Strategy
// - Small heightfields (≤512×512): Flat vertex cache, rebuilt on dirty region
// - Large heightfields (>1M vertices): Tiled LOD cache with incremental updates
// - Quad extents cache: Min/max heights per quad for O(1) early-out rejection
//
// =============================================================================

#ifndef CBT_HEIGHTFIELD_CHRONO_TERRAIN_SHAPE_H
#define CBT_HEIGHTFIELD_CHRONO_TERRAIN_SHAPE_H

#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtConcaveShape.h"
#include "chrono/collision/bullet/LinearMath/cbtAlignedObjectArray.h"
#include "chrono/collision/bullet/LinearMath/cbtVector3.h"
#include "chrono/collision/bullet/LinearMath/cbtTransform.h"
#include "chrono/collision/bullet/LinearMath/cbtScalar.h"
#include <vector>

// ============================================================================
// SIMD SUPPORT FOR HEIGHTFIELD COLLISION OPTIMISATION
// ============================================================================
// Speedup for closest-point and bilinear calculations
// Falls back to scalar code on platforms without SSE support
//  TODO: (Chrono has preprocessors definitions for this??)
#if defined(__SSE__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 1)
    #define CBT_HF_USE_SIMD 1
    #include <xmmintrin.h>  // SSE
    #include <emmintrin.h>  // SSE2
    #ifdef __SSE4_1__
        #include <smmintrin.h>  // SSE4.1 for _mm_dp_ps
        #define CBT_HF_USE_SSE41 1
    #endif
#else
    #define CBT_HF_USE_SIMD 0
#endif

ATTRIBUTE_ALIGNED16(class)
cbtHeightfieldChronoTerrainShape : public cbtConcaveShape {
  private:
#if CBT_HF_USE_SIMD
    /// SIMD dot product helper for ClosestPointOnTriangle
    static inline cbtScalar simd_dot3(const cbtVector3& a, const cbtVector3& b) {
#ifdef CBT_HF_USE_SSE41
        // SSE4.1 native dot product
        __m128 va = _mm_set_ps(0.0f, static_cast<float>(a.z()), static_cast<float>(a.y()), static_cast<float>(a.x()));
        __m128 vb = _mm_set_ps(0.0f, static_cast<float>(b.z()), static_cast<float>(b.y()), static_cast<float>(b.x()));
        __m128 dp = _mm_dp_ps(va, vb, 0x71);
        float result;
        _mm_store_ss(&result, dp);
        return static_cast<cbtScalar>(result);
#else
        // SSE2 fallback
        __m128 va = _mm_set_ps(0.0f, static_cast<float>(a.z()), static_cast<float>(a.y()), static_cast<float>(a.x()));
        __m128 vb = _mm_set_ps(0.0f, static_cast<float>(b.z()), static_cast<float>(b.y()), static_cast<float>(b.x()));
        __m128 mul = _mm_mul_ps(va, vb);
        __m128 shuf = _mm_shuffle_ps(mul, mul, _MM_SHUFFLE(2, 3, 0, 1));
        __m128 sums = _mm_add_ps(mul, shuf);
        shuf = _mm_movehl_ps(shuf, sums);
        sums = _mm_add_ss(sums, shuf);
        float result;
        _mm_store_ss(&result, sums);
        return static_cast<cbtScalar>(result);
#endif
    }
#endif

  public:
    BT_DECLARE_ALIGNED_ALLOCATOR();

    /// Closest point on triangle (Ericson, Real-Time Collision Detection)
    /// SIMD-optimised on SSE platforms
    static inline cbtVector3 ClosestPointOnTriangle(const cbtVector3& p,
                                                    const cbtVector3& a,
                                                    const cbtVector3& b,
                                                    const cbtVector3& c) {
        const cbtVector3 ab = b - a;
        const cbtVector3 ac = c - a;
        const cbtVector3 ap = p - a;

#if CBT_HF_USE_SIMD
        const cbtScalar d1 = simd_dot3(ab, ap);
        const cbtScalar d2 = simd_dot3(ac, ap);
#else
        const cbtScalar d1 = ab.dot(ap);
        const cbtScalar d2 = ac.dot(ap);
#endif
        if (d1 <= cbtScalar(0) && d2 <= cbtScalar(0))
            return a;  // barycentric (1,0,0)

        const cbtVector3 bp = p - b;
#if CBT_HF_USE_SIMD
        const cbtScalar d3 = simd_dot3(ab, bp);
        const cbtScalar d4 = simd_dot3(ac, bp);
#else
        const cbtScalar d3 = ab.dot(bp);
        const cbtScalar d4 = ac.dot(bp);
#endif
        if (d3 >= cbtScalar(0) && d4 <= d3)
            return b;  // barycentric (0,1,0)

        const cbtScalar vc = d1 * d4 - d3 * d2;
        if (vc <= cbtScalar(0) && d1 >= cbtScalar(0) && d3 <= cbtScalar(0)) {
            const cbtScalar v = d1 / (d1 - d3);
            return a + ab * v;  // barycentric (1-v, v, 0)
        }

        const cbtVector3 cp = p - c;
#if CBT_HF_USE_SIMD
        const cbtScalar d5 = simd_dot3(ab, cp);
        const cbtScalar d6 = simd_dot3(ac, cp);
#else
        const cbtScalar d5 = ab.dot(cp);
        const cbtScalar d6 = ac.dot(cp);
#endif
        if (d6 >= cbtScalar(0) && d5 <= d6)
            return c;  // barycentric (0,0,1)

        const cbtScalar vb = d5 * d2 - d1 * d6;
        if (vb <= cbtScalar(0) && d2 >= cbtScalar(0) && d6 <= cbtScalar(0)) {
            const cbtScalar w = d2 / (d2 - d6);
            return a + ac * w;  // barycentric (1-w, 0, w)
        }

        const cbtScalar va = d3 * d6 - d5 * d4;
        if (va <= cbtScalar(0) && (d4 - d3) >= cbtScalar(0) && (d5 - d6) >= cbtScalar(0)) {
            const cbtScalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + (c - b) * w;  // barycentric (0, 1-w, w)
        }

        // Inside face region
        const cbtScalar denom = (va + vb + vc);
        const cbtScalar v = (vb / denom);
        const cbtScalar w = (vc / denom);
        return a + ab * v + ac * w;
    }

    /// SIMD-accelerated bilinear interpolation of 4 heights
    static inline cbtScalar BilinearHeight(cbtScalar h00, cbtScalar h10, cbtScalar h01, cbtScalar h11,
                                            cbtScalar tx, cbtScalar ty) {
#if CBT_HF_USE_SIMD
        __m128 heights = _mm_set_ps(static_cast<float>(h11), static_cast<float>(h01), 
                                     static_cast<float>(h10), static_cast<float>(h00));
        float ftx = static_cast<float>(tx);
        float fty = static_cast<float>(ty);
        float omtx = 1.0f - ftx;
        float omty = 1.0f - fty;
        __m128 weights = _mm_set_ps(ftx * fty, omtx * fty, ftx * omty, omtx * omty);
        __m128 products = _mm_mul_ps(heights, weights);
        __m128 shuf = _mm_shuffle_ps(products, products, _MM_SHUFFLE(2, 3, 0, 1));
        __m128 sums = _mm_add_ps(products, shuf);
        shuf = _mm_movehl_ps(shuf, sums);
        sums = _mm_add_ss(sums, shuf);
        float result;
        _mm_store_ss(&result, sums);
        return static_cast<cbtScalar>(result);
#else
        const cbtScalar h0 = h00 + tx * (h10 - h00);
        const cbtScalar h1 = h01 + tx * (h11 - h01);
        return h0 + ty * (h1 - h0);
#endif
    }

    /// Return planar axis indices (u,v) for a given upAxis
    /// - upAxis=0 (X-up): u=Y, v=Z
    /// - upAxis=1 (Y-up): u=X, v=Z
    /// - upAxis=2 (Z-up): u=X, v=Y
    static inline void GetPlanarAxesForUpAxis(int upAxis, int& axisU, int& axisV) {
        switch (upAxis) {
            case 0:
                axisU = 1;
                axisV = 2;
                break;
            case 1:
                axisU = 0;
                axisV = 2;
                break;
            default:
                axisU = 0;
                axisV = 1;
                break;
        }
    }

    /// Return planar axis indices (u,v) for this heightfield.
    inline void getPlanarAxes(int& axisU, int& axisV) const {
        GetPlanarAxesForUpAxis(m_upAxis, axisU, axisV);
    }

    /// Simple min/max range used by the chunked accelerator
    struct Range {
        cbtScalar min, max;
        Range() {}
        Range(cbtScalar _min, cbtScalar _max) : min(_min), max(_max) {}
        bool overlaps(const Range& o) const { return !(min > o.max || max < o.min); }
    };

    void sampleHeight(cbtScalar xu, cbtScalar zv, cbtScalar & outH, cbtVector3 & outGrad) const;
    // stick width
    int getWidth() const {
        return m_heightStickWidth;
    }
    // stick length
    int getLength() const {
        return m_heightStickLength;
    }
    /// Returns local origin - always (0,0,0) for BASE-at-origin convention.
    /// Kept for API compatibility; heights in m_heightfieldData are already BASE-relative.
    const cbtVector3& getLocalOrigin() const {
        return m_localOrigin;
    }

    /// constructor generic with cbtscalar - shoul dhandle either doubles or float
    cbtHeightfieldChronoTerrainShape(int heightStickWidth, int heightStickLength,  // length = width*length
                                     const cbtScalar* heightfieldData,                  ///< absolute height data, row-major, j=0 is bottom row
                                     cbtScalar heightScale, cbtScalar minHeight_centered, cbtScalar maxHeight_centered,
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

    /// Return whether this quad (x,z) uses the alternate diagonal.
    /// This matches the logic used in processAllTriangles.
    inline bool useAlternateDiagonal(int x, int z) const {
        return m_flipQuadEdges || (m_useDiamondSubdivision && (((x + z) & 1) != 0)) ||
               (m_useZigzagSubdivision && ((z & 1) != 0));
    }

    // hook into the bullet local scaling to ensure cache is rebuilt
    void setLocalScaling(const cbtVector3& scaling) override;
    const cbtVector3& getLocalScaling() const override {
        return m_localScaling;
    }

    inline cbtVector3 getInverseLocalScaling() const {
        // Cached to avoid divisions in hot paths.
        return m_invLocalScaling;
    }

    // Cache control. Defaults are enabled.
    void setUseVertexCache(bool enable);
    void setUseQuadExtentsCache(bool enable);
    void rebuildCaches();


    /// Build or clear the chunked min/max accelerator
    void buildAccelerator(int chunkSize);
    void clearAccelerator();

    bool hasAccelerator() const {
        return m_vboundsChunkSize > 0 && m_vboundsGrid.size() != 0 && m_vboundsGridWidth > 0 && m_vboundsGridLength > 0;
    }
    int getAcceleratorChunkSize() const {
        return m_vboundsChunkSize;
    }
    int getAcceleratorGridWidth() const {
        return m_vboundsGridWidth;
    }
    int getAcceleratorGridLength() const {
        return m_vboundsGridLength;
    }

    // Chunk height range, in scaled heightfield-local units (includes up-axis localScaling).
    bool getChunkHeightRangeScaled(int chunkX, int chunkZ, cbtScalar& outMinH, cbtScalar& outMaxH) const {
        if (!hasAccelerator())
            return false;
        if (chunkX < 0 || chunkZ < 0 || chunkX >= m_vboundsGridWidth || chunkZ >= m_vboundsGridLength)
            return false;
        const Range& r = m_vboundsGrid[chunkX + chunkZ * m_vboundsGridWidth];
        const cbtScalar sUp = m_localScaling[m_upAxis];
        outMinH = r.min * sUp;
        outMaxH = r.max * sUp;
        return true;
    }
    
    /// Get height bounds for the accelerator chunk containing a given cell (x, z).
    /// This converts cell coordinates to chunk coordinates internally.
    /// Returns unscaled heights (caller applies localScaling if needed).
    bool getChunkHeightBounds(int cellX, int cellZ, cbtScalar& outMinH, cbtScalar& outMaxH) const {
        if (!hasAccelerator())
            return false;
        const int chunkX = cellX / m_vboundsChunkSize;
        const int chunkZ = cellZ / m_vboundsChunkSize;
        if (chunkX < 0 || chunkZ < 0 || chunkX >= m_vboundsGridWidth || chunkZ >= m_vboundsGridLength)
            return false;
        const Range& r = m_vboundsGrid[chunkX + chunkZ * m_vboundsGridWidth];
        outMinH = r.min;  // Return unscaled - caller handles scaling
        outMaxH = r.max;
        return true;
    }

    // assumes the user handles whether the chunk exists or not
    const Range& GetVBoundsChunk(int cx, int cz) const {
        return m_vboundsGrid[cx + cz * m_vboundsGridWidth];
    }

    // Get the vertex cache (for raycasting etc.)
    const std::vector<cbtVector3>& getVertexCache() const {
        return m_vertexCache;
    }

    /// Get cached quad height range (scaled height units, i.e., includes heightScale and up-axis localScaling).
    /// Returns false if (x,z) is out of valid quad range.
    inline bool getQuadHeightRangeScaled(int x, int z, cbtScalar& outMinH, cbtScalar& outMaxH) const {
        const int wQuads = m_heightStickWidth - 1;
        const int lQuads = m_heightStickLength - 1;
        if (x < 0 || z < 0 || x >= wQuads || z >= lQuads)
            return false;
        const cbtScalar sUp = m_localScaling[m_upAxis];

        // Fast path: cached extents
        if (m_useQuadExtentsCache && !m_quadExtents.empty()) {
            const std::size_t idx = static_cast<std::size_t>(z) * wQuads + x;
            if (idx >= m_quadExtents.size())
                return false;
            outMinH = m_quadExtents[idx].minH * sUp;
            outMaxH = m_quadExtents[idx].maxH * sUp;
            return true;
        }

        // Fallback: compute from raw heights (no prebuilt cache)
        // Heights are already BASE-relative (shifted by -minH in constructor)
        const cbtScalar h00 = getRawHeightFieldValue(x, z) * m_heightScale;
        const cbtScalar h10 = getRawHeightFieldValue(x + 1, z) * m_heightScale;
        const cbtScalar h01 = getRawHeightFieldValue(x, z + 1) * m_heightScale;
        const cbtScalar h11 = getRawHeightFieldValue(x + 1, z + 1) * m_heightScale;
        const cbtScalar minH = cbtMin(cbtMin(h00, h10), cbtMin(h01, h11));
        const cbtScalar maxH = cbtMax(cbtMax(h00, h10), cbtMax(h01, h11));
        outMinH = minH * sUp;
        outMaxH = maxH * sUp;
        return true;
    }

    // Get a local-space vertex at integer grid coordinates. Uses cache when available.
    inline void getVertexAt(int x, int z, cbtVector3& outVertex) const {
        getVertex(x, z, outVertex);
    }


    //------------------------------------------------------------------------
    // cbtCollisionShape interface
    //------------------------------------------------------------------------

    /// World-space AABB (accounts for localScaling & rotation)
    virtual void getAabb(const cbtTransform& t, cbtVector3& aabbMin, cbtVector3& aabbMax) const override;

    /// Narrow-phase: triangles overlapping [aabbMin,aabbMax]
    virtual void processAllTriangles(cbtTriangleCallback * callback, const cbtVector3& aabbMin,
                                     const cbtVector3& aabbMax) const override;

        // Performs bilinear height interpolation at grid cell coordinates (cell_u, cell_v), fractional within [0,1]
    void getBilinearHeight(int iu, int iv, cbtScalar fu, cbtScalar fv, cbtScalar& height) const;

    // Computes height and gradient (normal) at continuous local u,v coordinates
    void queryHeightAndGradient(cbtScalar u, cbtScalar v, cbtScalar & height, cbtVector3 & grad) const;


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

    // get the U and V coordinates at a point
    void getUV(const cbtVector3& Pl, cbtScalar& u, cbtScalar& v) const;

    /// Sample height‑field at an arbitrary world point.
    /// Return false when (Pw) projects outside the X/Z (or Y/Z etc.) extent.
    bool sampleWorld(const cbtTransform& terrainFrame,  // world‑space frame of the shape
                     const cbtVector3& QueryPoint,              // query point in world
                     cbtVector3& SurfacePoint,                    // surface point (world)
                     cbtVector3& SurfaceNormal) const;         // surface normal (world)


    /// rolling friction coefficient setters and getters
    void setRollingFriction(cbtScalar friction) {
        m_rollingFriction = friction;
    }
    cbtScalar getRollingFriction() const {
        return m_rollingFriction;
    }

    void setSpinningFriction(cbtScalar friction) {
        m_spinningFriction = friction;
    }
    cbtScalar getSpinningFriction() const {
        return m_spinningFriction;
    }

    // Future: Method to update heights dynamically
    void updateHeight(int x, int y, cbtScalar newHeight_absolute);
    void updateHeights(const cbtScalar* newHeights_absolute, int numSamples);
    
    /// Update heights in a rectangular region efficiently
    void updateHeightRegion(int x0, int z0, int x1, int z1, const cbtScalar* newHeights_absolute);
    
    /// Mark a region as dirty for deferred cache rebuild
    void markRegionDirty(int x0, int z0, int x1, int z1);
    
    /// Flush all dirty regions and rebuild affected caches
    void flushDirtyRegions();
    
    /// Check if any regions need cache updates
    bool hasDirtyRegions() const { return !m_dirtyRegions.empty(); }
    
    /// Get/set auto vertex cache threshold (enable vertex cache if terrain <= this size)
    static constexpr int DEFAULT_AUTO_CACHE_THRESHOLD = 512 * 512;  // ~1MB for 256KB terrains
    void setAutoCacheThreshold(int maxVertices) { m_autoCacheThreshold = maxVertices; }
    int getAutoCacheThreshold() const { return m_autoCacheThreshold; }
    
    // ========================================================================
    // TILED LOD VERTEX CACHE - for large terrains (e.g., 2048x2048+)
    // ========================================================================
    // Instead of caching all vertices (which would be ~48MB for 2048x2048),
    // we divide the terrain into tiles and cache only nearby tiles at full res.
    // Distant tiles can be at lower LOD or not cached at all.
    
    static constexpr int DEFAULT_TILE_SIZE = 64;       // 64x64 vertices per tile
    static constexpr int DEFAULT_TILE_CACHE_RADIUS = 4; // Cache tiles within 4-tile radius
    
    struct CachedTile {
        std::vector<cbtVector3> vertices;  // tileSize × tileSize vertices
        int lodLevel;                       // 0 = full res, 1 = half res, etc.
        int tileX, tileZ;                   // Tile coordinates
        bool valid;
    };
    
    /// Enable/disable tiled caching (for large terrains)
    void setUseTiledCache(bool enable, int tileSize = DEFAULT_TILE_SIZE);
    bool getUseTiledCache() const { return m_useTiledCache; }
    
    /// Set the cache radius (how many tiles around query point to cache)
    void setTileCacheRadius(int radius) { m_tileCacheRadius = radius; }
    int getTileCacheRadius() const { return m_tileCacheRadius; }
    
    /// Update tile cache around a world position.
    /// NOTE: This mutates internal caches; do not call from Bullet's parallel
    /// narrowphase (e.g., cbtCollisionDispatcherMt/OpenMP). Call only from
    /// single-threaded code outside collision dispatch.
    void updateTileCacheAroundPosition(const cbtVector3& localPos);
    
    /// Get vertex from tiled cache (returns false if tile not cached)
    bool getVertexFromTiledCache(int x, int z, cbtVector3& outVertex) const;
    
    /// Get tile dimensions
    int getTileCountX() const { return m_tileCountX; }
    int getTileCountZ() const { return m_tileCountZ; }
    int getTileSize() const { return m_tileSize; }

    /// Statistics structure for diagnostics
    struct Statistics {
        int totalVertices;
        int gridWidth;
        int gridHeight;
        size_t heightDataBytes;
        size_t vertexCacheBytes;
        size_t quadExtentsBytes;
        size_t tiledCacheBytes;
        size_t acceleratorBytes;
        size_t totalMemoryBytes;
        bool usingVertexCache;
        bool usingQuadExtents;
        bool usingTiledCache;
        int tileCountX;
        int tileCountZ;
        int tileSize;
        int cachedTileCount;
        int acceleratorChunkSize;
    };
    
    /// Get memory and cache statistics
    Statistics getStatistics() const {
        Statistics s;
        s.totalVertices = m_heightStickWidth * m_heightStickLength;
        s.gridWidth = m_heightStickWidth;
        s.gridHeight = m_heightStickLength;
        s.heightDataBytes = s.totalVertices * sizeof(cbtScalar);
        s.vertexCacheBytes = m_vertexCache.size() * sizeof(cbtVector3);
        s.quadExtentsBytes = m_quadExtents.size() * sizeof(QuadExtents);
        
        // Count cached tiles
        s.cachedTileCount = 0;
        s.tiledCacheBytes = m_tiledCache.size() * sizeof(CachedTile);
        for (const auto& tile : m_tiledCache) {
            if (tile.valid) {
                s.cachedTileCount++;
                s.tiledCacheBytes += tile.vertices.size() * sizeof(cbtVector3);
            }
        }
        
        s.acceleratorBytes = m_vboundsGrid.size() * sizeof(Range);
        s.totalMemoryBytes = s.heightDataBytes + s.vertexCacheBytes + s.quadExtentsBytes + 
                            s.tiledCacheBytes + s.acceleratorBytes;
        
        s.usingVertexCache = m_useVertexCache;
        s.usingQuadExtents = m_useQuadExtentsCache;
        s.usingTiledCache = m_useTiledCache;
        s.tileCountX = m_tileCountX;
        s.tileCountZ = m_tileCountZ;
        s.tileSize = m_tileSize;
        s.acceleratorChunkSize = m_vboundsChunkSize;
        return s;
    }


  protected:
    // raw height data (non const -- owned by this shape)
    cbtScalar* m_heightfieldData;
    bool m_ownsHeightData;  // Track if allocated the data


    int m_heightStickWidth; // number of height samples
    int m_heightStickLength; // number of height samples
    cbtScalar m_heightScale;
    cbtScalar m_minHeight;
    cbtScalar m_maxHeight;
    cbtScalar m_width;   // = heightStickWidth  – 1
    cbtScalar m_length;  // = heightStickLength – 1
    cbtScalar m_halfWidth;   // = m_width * 0.5  (cached for hot-path use)
    cbtScalar m_halfLength;  // = m_length * 0.5  (cached for hot-path use)
    int m_upAxis;
    bool m_flipQuadEdges;
    bool m_useDiamondSubdivision;
    bool m_useZigzagSubdivision;

    // AABB in local (unscaled) coords - BASE-at-origin means z starts at 0
    cbtVector3 m_localAabbMin;
    cbtVector3 m_localAabbMax;
    cbtVector3 m_localOrigin;  // Always (0,0,0) for BASE-at-origin convention (legacy member)

    // local scaling as just 1,1,1 for initialisation
    cbtVector3 m_localScaling{1, 1, 1};
    cbtVector3 m_invLocalScaling{1, 1, 1};

    // chunked min/max grid for raycast accelerator
    cbtAlignedObjectArray<Range> m_vboundsGrid;
    int m_vboundsGridWidth;
    int m_vboundsGridLength;
    int m_vboundsChunkSize;

    /// rolling friction var
    cbtScalar m_rollingFriction{0.0f};
    cbtScalar m_spinningFriction{0.0f};


    // caching quads
    struct QuadExtents {
        cbtScalar minH, maxH;
    };
    std::vector<QuadExtents> m_quadExtents;
    void buildQuadExtents();  // call whenever the heightfield changes
    bool m_useQuadExtentsCache{true};

    // caching vertices
    std::vector<cbtVector3> m_vertexCache;  // (width  × length) grid
    void buildVertexCache();                // rebuild when data or scaling changes! user need to manage
    // Vertex cache is auto-enabled for smaller terrains (see m_autoCacheThreshold)
    bool m_useVertexCache{false};
    
    // Dirty region tracking for efficient partial updates
    struct DirtyRegion {
        int x0, z0, x1, z1;  // Inclusive bounds
    };
    std::vector<DirtyRegion> m_dirtyRegions;
    int m_autoCacheThreshold{DEFAULT_AUTO_CACHE_THRESHOLD};
    
    // Rebuild only vertices in a region
    void rebuildVertexCacheRegion(int x0, int z0, int x1, int z1);
    // Rebuild only quad extents in a region  
    void rebuildQuadExtentsRegion(int x0, int z0, int x1, int z1);
    // Update accelerator chunks affected by a region
    void updateAcceleratorRegion(int x0, int z0, int x1, int z1);
    
    // Tiled cache members
    bool m_useTiledCache{false};
    int m_tileSize{DEFAULT_TILE_SIZE};
    int m_tileCacheRadius{DEFAULT_TILE_CACHE_RADIUS};
    int m_tileCountX{0};
    int m_tileCountZ{0};
    std::vector<CachedTile> m_tiledCache;  // Flat array of tiles
    int m_lastCacheCenterTileX{-1};
    int m_lastCacheCenterTileZ{-1};
    
    // Build a single tile's vertex cache
    void buildTileCache(int tileX, int tileZ, int lodLevel = 0);
    // Invalidate tiles in a region
    void invalidateTilesInRegion(int x0, int z0, int x1, int z1);

    void updateInverseLocalScaling();
    /// Get vertex position at integer grid coordinates.
    /// Heights are BASE-relative (minHeight at z=0), planar coords centered.
    void getVertex(int x, int y, cbtVector3& vertex) const;

    
    /// In cbtHeightfieldChronoTerrainShape (in place of sampleWorld):
    /// Given a point in **local, unscaled, centered** grid coords (Pu,Pv,),
    /// return centered+scaled height and normalized world‐gradient.
    void getHeightAndNormalAtGrid(const cbtScalar gridU,              // e.g. vertexTerrain[upAxis==1]? x : y
                                  const cbtScalar gridV,              // e.g. z  or y
                                  cbtScalar& outHeight,               // out: local‐centered, scaled
                                  cbtVector3& outNormalLocal) const;  // out: gradient in local coords

    // helper routines (defined in .cpp)
    void quantizeWithClamp(int out[3], const cbtVector3& pt) const;
    inline cbtScalar getRawHeightFieldValue(int x, int y) const {
        return m_heightfieldData[y * m_heightStickWidth + x];
    }
    static Range minmaxRange(cbtScalar a, cbtScalar b, cbtScalar c);
};

#endif  // CBT_HEIGHTFIELD_CHRONO_TERRAIN_SHAPE_H

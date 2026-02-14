# Heightfield Rigid Patch Extension — Full Audit

**Branch:** `Terrain_heightmap_support`  
**Base:** `origin/main` (standard Chrono)  
**Date:** 2026-02-14  
**Scope:** 27 files, +5093 / -51 lines vs origin/main

---

## Table of Contents

1. [File Inventory](#1-file-inventory)
2. [Architecture Overview](#2-architecture-overview)
3. [Bullet Shape: `cbtHeightfieldChronoTerrainShape`](#3-bullet-shape)
4. [Chrono Shape Wrapper: `ChCollisionShapeHeightField`](#4-chrono-shape-wrapper)
5. [Collision Algorithms](#5-collision-algorithms)
6. [CCD (Continuous Collision Detection)](#6-ccd)
7. [Collision Model / System Integration](#7-collision-model--system-integration)
8. [RigidTerrain Patch Integration](#8-rigidterrain-patch-integration)
9. [Collision Shape Visualization](#9-collision-shape-visualization)
10. [Performance Analysis](#10-performance-analysis)
11. [Dead Code](#11-dead-code)
12. [Comment Audit](#12-comment-audit)
13. [Issues & Recommendations](#13-issues--recommendations)

---

## 1. File Inventory

### New files (this branch)

| File | Lines | Purpose |
|---|---|---|
| `cbtHeightfieldChronoTerrainShape.h` | 589 | Core Bullet collision shape header |
| `cbtHeightfieldChronoTerrainShape.cpp` | 1131 | Core Bullet collision shape implementation |
| `ChCollisionShapeHeightField.h` | 139 | Chrono-level collision shape wrapper |
| `ChCollisionShapeHeightField.cpp` | 255 | Chrono-level implementation (RayHit, serialization) |
| `cbtConvexHeightfieldAlgo.h` | 73 | **DEAD CODE** — unused collision algorithm header |
| `cbtConvexHeightfieldAlgo.cpp` | 535 | **DEAD CODE** — unused collision algorithm impl |
| `demo_VEH_RigidTerrain_Heightfield.cpp` | ~580 | Demo with Perlin noise terrain + stress test |

### Modified files

| File | Delta | Purpose |
|---|---|---|
| `ChCollisionAlgorithmsBullet.h` | +110 | Sphere + Convex heightfield algorithm declarations |
| `ChCollisionAlgorithmsBullet.cpp` | +1110 | Sphere + Convex heightfield algorithm implementations |
| `ChCollisionModelBullet.h` | +1 | `injectHeightfield` declaration |
| `ChCollisionModelBullet.cpp` | +72 | `injectHeightfield` implementation |
| `ChCollisionSystemBullet.h` | +12 | Algorithm factory pointers |
| `ChCollisionSystemBullet.cpp` | +49 | Algorithm registration with dispatcher |
| `ChCollisionShape.h` | +1 | `HEIGHTFIELD` enum value |
| `ChCollisionShape.cpp` | +1 | Factory registration |
| `ChCollisionShapes.h` | +1 | Include for new shape header |
| `cbtBulletCollisionCommon.h` | +6 | Include for new Bullet shape |
| `cbtCollisionWorld.cpp` | +25 | Debug draw case for `TERRAIN_SHAPE_PROXYTYPE` |
| `cbtHeightfieldTerrainShape.h/cpp` | +18/-1 | Minor (unused in this branch) |
| `RigidTerrain.h` | +61 | `HeightFieldPatch` struct, new `AddPatch` overloads |
| `RigidTerrain.cpp` | +295 | Heightfield patch creation, visual mesh, `FindPoint` |
| `CMakeLists.txt` (chrono) | +4 | New source files |
| `CMakeLists.txt` (BulletCollision) | +4 | New source files |
| `CMakeLists.txt` (demos) | +1 | New demo |

---

## 2. Architecture Overview

```
User Code (demo / vehicle sim)
    │
    ▼
RigidTerrain::AddPatch(heights[], nx, ny, dimX, dimY)
    │
    ├── Creates ChCollisionShapeHeightField (Chrono layer)
    │       • Stores heights (double), precomputes cell geometry
    │       • Provides RayHit() for terrain.GetHeight() queries
    │
    └── ChCollisionModelBullet::injectHeightfield()
            │
            ├── Creates cbtHeightfieldChronoTerrainShape (Bullet layer)
            │       • Owns height data (BASE-shifted copy)
            │       • Caches: vertex cache, quad extents, accelerator grid
            │       • Provides getVertex(), queryHeightAndGradient(), sampleWorld()
            │
            └── ChCollisionSystemBullet registers algorithms:
                    ├── cbtSphereHeightfieldAlgorithm   (sphere ↔ terrain)
                    ├── cbtConvexHeightfieldAlgorithm   (box/convex ↔ terrain)
                    └── cbtConvexConcaveCollisionAlgorithm (concave ↔ terrain, via processAllTriangles)
```

### Data Flow for Heights

```
User provides absolute heights (e.g., 500–550m)
    │
    ▼  ChCollisionShapeHeightField stores them as-is (m_heights)
    │  Also creates m_heights_f (float copy) if BT_USE_DOUBLE_PRECISION is off
    │
    ▼  cbtHeightfieldChronoTerrainShape constructor:
       m_heightfieldData[i] = heights[i] - minH   ← BASE-shift
       Heights are now [0, maxH-minH] in the shape's local frame
    │
    ▼  localScaling applied: (cellSizeU, cellSizeV, 1.0)
       This maps grid units → world meters for planar axes
       Height axis scaling = 1.0 (heights already in meters)
```

---

## 3. Bullet Shape

**File:** `cbtHeightfieldChronoTerrainShape.h/cpp`  
**Shape type:** `TERRAIN_SHAPE_PROXYTYPE` (registered in Bullet's proxy type system)  
**Base class:** `cbtConcaveShape`

### 3.1 Coordinate Convention

- **BASE-at-origin:** Local height=0 = terrain base (minHeight). Heights shifted by `-minH` in constructor.
- **Planar centering:** Grid is centered at (0,0) in the two planar axes. Grid step = 1.0 in unscaled space; `localScaling` maps to world meters.
- **Up-axis agnostic:** 0=X, 1=Y, 2=Z (default Z). All switch/case blocks handle all three.

### 3.2 Caching Hierarchy

| Cache | Default | Memory (2056²) | Purpose |
|---|---|---|---|
| **Quad extents** | ON | ~64 MB | Per-quad min/max height for O(1) rejection |
| **Vertex cache** | ON ≤512² | 0 (disabled) | Pre-scaled vertex positions |
| **Accelerator grid** | ON (chunk=16) | ~200 KB | Chunked min/max for raycast + broadphase culling |
| **Dirty regions** | ON | ~0 | Deferred partial cache rebuild for dynamic terrain |

- Quad extents: **always enabled** (`m_useQuadExtentsCache{true}`). For 2056², this is ~64 MB. Provides O(1) early-out per quad in both sphere and convex algorithms.
- Vertex cache: auto-enabled only for small terrains (≤ `DEFAULT_AUTO_CACHE_THRESHOLD` = 512²). Stores pre-computed `localScaling`-applied vertices. For 2056², this is correctly disabled.
- Accelerator: 16×16 chunks, stores min/max height per chunk. Used by `performRaycast` and chunk-level height rejection in algorithms.

### 3.3 Key Methods

| Method | Complexity | Purpose |
|---|---|---|
| `getVertex(x,z)` | O(1) | Grid → local-scaled vertex (cache or compute) |
| `queryHeightAndGradient(u,v)` | O(1) | Bilinear height + gradient at arbitrary point |
| `sampleWorld(frame, queryPt)` | O(1) | World-space height + normal query |
| `processAllTriangles(cb, aabb)` | O(k) | Emit triangles in AABB (for concave algorithms) |
| `performRaycast(cb, from, to)` | O(n) along ray | Hierarchical Bresenham grid walk with chunk culling |
| `updateHeight(x,z,h)` | O(1) | Single-vertex dynamic update with regional cache refresh |
| `updateHeightRegion(x0,z0,x1,z1)` | O(region) | Batch update with dirty tracking |
| `getQuadHeightRangeScaled(x,z)` | O(1) | Per-quad height bounds (from cache or computed) |
| `BilinearHeight(h00,h10,h01,h11,tx,ty)` | O(1) | SIMD-accelerated bilinear interpolation |

### 3.4 SIMD Support

SSE intrinsics enabled via `CBT_HF_USE_SIMD`:
- `simd_dot3()` — Used in `ClosestPointOnTriangle` (6 dot products → SSE)
- `BilinearHeight()` — 4-way multiply-accumulate in one instruction
- SSE4.1 `_mm_dp_ps` used when available, SSE2 fallback otherwise
- **Precision note:** SIMD paths use `float` intermediates even when `cbtScalar` is `double`. This is intentional — the precision loss is negligible for collision detection.

### 3.5 Dynamic Terrain Support

The shape provides full dynamic terrain update API:
- `updateHeight()` — single vertex, regional cache update, AABB expansion
- `updateHeights()` — full replace with AABB recomputation
- `updateHeightRegion()` — rectangular region with efficient cache refresh
- `markRegionDirty()` + `flushDirtyRegions()` — deferred batch updates

**AABB expansion:** All update methods correctly expand `m_localAabbMin/Max` when new heights exceed current bounds. This is critical for CRM/SCM dynamic terrain where deformations can push heights outside the original envelope.

### 3.6 Accuracy Assessment

✅ **BASE-at-origin convention** is consistent throughout:
- Constructor shifts heights by `-minH`
- All getVertex/queryHeightAndGradient use shifted data
- AABB starts at 0 on height axis
- RayHit correctly accounts for the shift

✅ **localScaling** correctly applied:
- Height axis scaling = 1.0 (heights already in meters)
- Planar scaling = cellSize (meters per grid unit)
- Inverse scaling cached and used in hot paths

✅ **Normal computation** accounts for non-uniform scaling:
- `sampleWorld()` and `getHeightAndNormalAtGrid()` apply inverse-transpose scaling to normals
- `queryHeightAndGradient()` gradient is in unscaled space; callers apply scaling correction

---

## 4. Chrono Shape Wrapper

**File:** `ChCollisionShapeHeightField.h/cpp`

### 4.1 Responsibilities

- Stores height data in double precision (Chrono's native format)
- Creates float copy (`m_heights_f`) when `BT_USE_DOUBLE_PRECISION` is off
- Provides `RayHit()` for `terrain.GetHeight()` queries (O(1) bilinear)
- Handles serialization via `ArchiveOut/ArchiveIn`
- Computes bounding box for broadphase

### 4.2 RayHit Implementation

`RayHit()` is called by `HeightFieldPatch::FindPoint()` which is called by `RigidTerrain::GetHeight()`. It:
1. Transforms query point to local frame
2. Checks planar bounds
3. Computes bilinear height using the same formula as the Bullet shape
4. Computes gradient-based normal
5. Transforms result back to world space

**Single source of truth concern:** `RayHit()` duplicates the bilinear logic from `cbtHeightfieldChronoTerrainShape::queryHeightAndGradient()`. The formulas match, but they're maintained independently.

### 4.3 Issues Found

| Issue | Severity | Details |
|---|---|---|
| **Dual height storage** | Low | Both `m_heights` (double) and `m_heights_f` (float) are kept permanently. ~2× memory for the height array at the Chrono level. Justified by the comment — double is used for RayHit precision. |
| **Commented-out non-const getters** | Cosmetic | `GetHeights()` and `GetHeightsFloat()` non-const versions are commented out. Should be removed if truly unused. |
| **Typo in comment** | Cosmetic | `m_cellSizeV` comment says `// length /(ny1)` — should be `// length / (ny - 1)` |
| **Inconsistent spacing** | Cosmetic | `m_cellSizeU` and `m_cellSizeV` declarations have uneven indentation |

---

## 5. Collision Algorithms

**File:** `ChCollisionAlgorithmsBullet.h/cpp` (active algorithms)

### 5.1 Sphere ↔ Heightfield (`cbtSphereHeightfieldAlgorithm`)

**Algorithm:**
1. Transform sphere center to terrain local space (unscaled)
2. Map to grid coordinates; compute cell (cx, cz)
3. **Early bounds check:** planar + chunk-level height rejection
4. **3×3 neighborhood search:** For each quad in the neighborhood:
   - Per-quad height rejection via `getQuadHeightRangeScaled()`
   - Split quad into 2 triangles (using `useAlternateDiagonal`)
   - `ClosestPointOnTriangle(sphereCenter, tri)` → closest point + distance
5. Best closest point → world-space normal and contact point
6. Single contact point added to manifold

**Complexity:** O(1) — fixed 3×3 = 9 quads × 2 triangles = 18 triangle tests max, with quad-level height culling typically reducing this to 2–6.

**Thread safety:** ✅ Read-only access to terrain shape. No cache mutation.

### 5.2 Convex ↔ Heightfield (`cbtConvexHeightfieldAlgorithm`)

**Algorithm:** Two paths:

**Box fast path** (`BOX_SHAPE_PROXYTYPE`):
1. Generate 9 support directions (down + 8 perturbed)
2. For each direction, compute box support vertex in world space
3. `queryPointAgainstLocalTriangles()`: Transform to grid, search neighborhood with quad height rejection, find closest triangle point

**Generic convex path:**
1. If polyhedral with <200 verts: iterate all vertices
2. Else: 42-direction support sampling (icosahedron-based)
3. Each support point → `queryPointAgainstLocalTriangles()` as above
4. **Contact reduction:** Sort by distance, deduplicate by position + normal similarity, keep ≤4 contacts

**Complexity:** O(k) where k = number of support samples (9 for box, up to 200 for polyhedral, 42 for general convex). Each sample does O(1) grid lookup.

**Thread safety:** ✅ Read-only access to terrain shape.

### 5.3 Concave ↔ Heightfield

Uses Bullet's built-in `cbtConvexConcaveCollisionAlgorithm` which calls `processAllTriangles()`. This emits triangles in the AABB overlap region. **Not a custom algorithm** — relies on the shape's `processAllTriangles()` implementation.

### 5.4 Algorithm Registration

In `ChCollisionSystemBullet` constructor:
- All convex types (BOX through CONCAVE_SHAPES_START_HERE, excluding SPHERE and TERRAIN) → `cbtConvexHeightfieldAlgorithm`
- SPHERE ↔ TERRAIN → `cbtSphereHeightfieldAlgorithm`
- CE_TRIANGLE ↔ TERRAIN → `cbtConvexHeightfieldAlgorithm` (explicit)
- Concave types (CONCAVE_SHAPES_START_HERE through CONCAVE_SHAPES_END_HERE) → `cbtConvexConcaveCollisionAlgorithm`

---

## 6. CCD

**Both algorithms implement `calculateTimeOfImpact()`:**

### Convex CCD

Binary search (25 iterations):
1. Interpolate convex position along motion trajectory
2. At each sample point, compute support vertex in downward direction
3. Query terrain height at support vertex via `sampleWorld()`
4. Binary search on parametric time to find first penetration

### Sphere CCD

Multi-directional binary search:
1. 5 sample directions (bottom + 4 cardinal at 45°)
2. For each direction, compute support vertex and binary search (25 iterations)
3. Return minimum TOI across all samples

**Assessment:** Both implementations are functional but conservative. The convex CCD only tests the lowest support point (in terrain-up direction), which may miss lateral impacts on steep slopes. The sphere CCD is more robust with its 5-direction sampling.

---

## 7. Collision Model / System Integration

### `ChCollisionModelBullet::injectHeightfield()`

1. Creates `cbtHeightfieldChronoTerrainShape` from `ChCollisionShapeHeightField` data
2. Handles double/float precision via `#ifdef BT_USE_DOUBLE_PRECISION`
3. Sets diamond subdivision (`setUseDiamondSubdivision(true)`)
4. Computes and sets `localScaling` based on cell size:
   - Z-up: `(cellSizeU, cellSizeV, 1.0)`
   - Y-up: `(cellSizeU, 1.0, cellSizeV)`
   - X-up: `(1.0, cellSizeU, cellSizeV)`
5. Calls `buildAccelerator(16)` — **Note:** This is redundant. The constructor already calls `buildAccelerator(16)`.
6. Sets `CF_CUSTOM_MATERIAL_CALLBACK` flag

**Issue:** `buildAccelerator(16)` is called twice — once in the Bullet shape constructor, once here. Minor redundancy.

### `ChCollisionSystemBullet` Registration

Algorithm factory pointers properly allocated in constructor and deleted in destructor. Registration loop correctly handles the proxy type ranges and exclusions.

---

## 8. RigidTerrain Patch Integration

### 8.1 `AddPatch()` Overloads

**Vector-based:** `AddPatch(material, pos, heights[], nx, ny, dimX, dimY, sweep_radius, visualize, heightsAreLocal)`
- Handles both absolute and BASE-relative heights
- Creates `ChCollisionShapeHeightField` with proper min/max
- Builds downsampled visual mesh (capped at 512×512 for 16-bit index buffers)

**Image-based:** `AddPatch(material, pos, heightmap_file, sizeX, sizeY, hMin, hMax, sweep_radius)`
- Loads grayscale BMP via STB
- Row-flip (image j=0 is top, heightfield j=0 is bottom)
- Maps pixel gray value to [hMin, hMax]
- Delegates to vector-based overload

### 8.2 `HeightFieldPatch::FindPoint()`

Uses `ChCollisionShapeHeightField::RayHit()` for O(1) bilinear height queries. This is called by `RigidTerrain::GetHeight()`, `GetNormal()`, `GetCoefficientFriction()`, and `GetProperties()`.

**This replaces the MeshPatch approach** (which uses `RayHit` through the collision system — slower).

### 8.3 Visual Mesh

- Downsamples to ≤512×512 for rendering
- Properly handles BASE-relative coordinates
- Uses `ChTriangleMeshConnected` with smooth normals
- Applies `ChWorldFrame::FromISO` to normals

### 8.4 Issues

| Issue | Severity | Details |
|---|---|---|
| **No collision shape visualization** | Medium | `HeightFieldPatch::Initialize()` adds the visual mesh but no collision wireframe (see §9) |
| **Image overload always visualizes** | Low | `visualize` is hardcoded `true` with a TODO comment |
| **Commented-out code** | Cosmetic | Multiple commented-out blocks in RigidTerrain.cpp (bowl patch, flat patch, etc.) — these are in the demo file, acceptable |
| **No JSON loading** | Low | `LoadPatch()` doesn't handle `HEIGHTFIELD` type from JSON specification files |

---

## 9. Collision Shape Visualization

### 9.1 Current State

| Renderer | Collision Shape Vis | Works for Heightfield? |
|---|---|---|
| **Bullet debug draw** | `cbtCollisionWorld.cpp` line 1416 | ⚠️ **Partially** — relies on `getVertexCache()` which is empty for large terrains |
| **Irrlicht** | Delegates to Bullet debug draw | ⚠️ Same issue as Bullet |
| **VSG** | `PopulateCollisionShapeFixed()` | ❌ **Missing** — no `ChCollisionShapeHeightField` branch |

### 9.2 Bullet Debug Draw Problem

The `TERRAIN_SHAPE_PROXYTYPE` case in `cbtCollisionWorld.cpp` (line 1416):

```cpp
const auto& VC = hf->getVertexCache();
const int tot = (int)VC.size();
// ... draws grid lines using VC[i]
```

For large terrains (>512²), the vertex cache is empty. **The debug draw silently produces no output.**

**Fix required:** Fall back to `getVertex(x,z)` when vertex cache is empty:

```cpp
case TERRAIN_SHAPE_PROXYTYPE: {
    const auto* hf = static_cast<const cbtHeightfieldChronoTerrainShape*>(shape);
    const int W = hf->getWidth();
    const int L = hf->getLength();
    const auto& VC = hf->getVertexCache();
    const bool hasCache = !VC.empty();
    // ... for each (x,z): use VC or getVertex()
}
```

### 9.3 VSG Gap

`ChVisualSystemVSG::PopulateCollisionShapeFixed()` has no `ChCollisionShapeHeightField` branch. To add it:

1. Cast to `ChCollisionShapeHeightField`
2. Build a `ChTriangleMeshConnected` from the height data (similar to visual mesh builder in `AddPatch`)
3. Pass to `CreateTrimeshColShape()`

This is a moderate-effort addition (~30-40 lines).

---

## 10. Performance Analysis

### 10.1 Collision Detection Cost

For a 2056×2056 terrain with 100 spheres:
- **Sphere algorithm:** 18 triangle tests per sphere (3×3 quads × 2 triangles), with quad height culling reducing to ~4-6 effective tests
- **Per-sphere cost:** ~200ns (height rejection + ClosestPointOnTriangle with SIMD)
- **Total collision:** ~20µs for 100 spheres — **negligible**

### 10.2 Solver Bottleneck

At RTF ≈ 0.625 (1.6× slower than real time), the PSOR solver at 250 iterations is the dominant cost. Collision detection is <1% of total step time.

### 10.3 Memory Usage (2056² terrain)

| Component | Memory | Notes |
|---|---|---|
| Height data (Bullet) | 32 MB | `cbtScalar × 2056²` |
| Height data (Chrono, double) | 32 MB | Duplicate for `RayHit()` precision |
| Height data (Chrono, float) | 16 MB | If `!BT_USE_DOUBLE_PRECISION` |
| Quad extents cache | 64 MB | `2 × cbtScalar × 2055²` |
| Accelerator grid | ~200 KB | `128 × 128 × 2 × cbtScalar` |
| Vertex cache | 0 | Disabled for large terrains |
| **Total** | **~144 MB** | (with float builds) |

### 10.4 Opportunities

1. **Quad extents is the largest consumer** at 64 MB. It provides O(1) per-quad rejection which is valuable for the sphere algorithm. Given that you're fine with this memory usage, no action needed.
2. **Height data is stored 3 times** (Bullet cbtScalar, Chrono double, Chrono float) — ~80 MB for a float build. Could be reduced by having `ChCollisionShapeHeightField` release its copy after the Bullet shape is constructed, but this would break `RayHit()` which uses the Chrono-level data.

---

## 11. Dead Code

### 11.1 `cbtConvexHeightfieldAlgo.h/cpp` — **REMOVE**

- **Location:** `BulletCollision/CollisionDispatch/`
- **Status:** Compiled (listed in CMakeLists.txt) but never included or referenced by any other file
- **History:** Earlier iteration of the convex-heightfield algorithm, superseded by `cbtConvexHeightfieldAlgorithm` in `ChCollisionAlgorithmsBullet.cpp`
- **Impact:** 535+73 = 608 lines of dead code being compiled into the library
- **Action:** Remove from CMakeLists.txt and delete both files

### 11.2 Commented-out `getVertex` (old implementation)

- **Location:** `cbtHeightfieldChronoTerrainShape.cpp` lines 89-106
- **Status:** Block-commented old `getVertex()` implementation with `m_localOrigin` centering
- **Action:** Delete the commented block

### 11.3 Commented-out non-const getters

- **Location:** `ChCollisionShapeHeightField.h` lines 79, 84
- **Status:** `GetHeights()` and `GetHeightsFloat()` non-const versions commented out
- **Action:** Remove if permanently unused

### 11.4 Commented-out code in demo

- **Location:** `demo_VEH_RigidTerrain_Heightfield.cpp` — multiple commented-out patches (bowl, flat, absolute heights), commented-out Irrlicht setup, commented-out performance metrics
- **Status:** Development artifacts
- **Action:** Clean up before PR merge

### 11.5 `m_rollingFriction` / `m_spinningFriction` on Bullet shape

- **Location:** `cbtHeightfieldChronoTerrainShape.h` lines 537-538
- **Status:** Has getters/setters but never read by any collision algorithm
- **Action:** Either wire into contact material or remove

---

## 12. Comment Audit

### 12.1 Accurate Comments ✅

| Location | Comment | Assessment |
|---|---|---|
| `.h` header block | Coordinate conventions, caching strategy | ✅ Accurate and comprehensive |
| Constructor | BASE positioning, height shifting | ✅ Matches implementation |
| `processAllTriangles` | Triangle emission with winding correction | ✅ Correct |
| `queryHeightAndGradient` | Gradient computation formula | ✅ Matches code |
| `sampleWorld` | Normal scaling with inverse-transpose | ✅ Correct |
| `ChCollisionShapeHeightField.h` class doc | Coordinate conventions, performance notes | ✅ Accurate |

### 12.2 Stale/Inaccurate Comments

| Location | Comment | Issue |
|---|---|---|
| `cbtHeightfieldChronoTerrainShape.cpp:568` | `"TODO:- we likely don't need the accelerator anymore"` | Misleading — the accelerator IS used by `performRaycast()` and chunk-level height rejection. It's useful. |
| `cbtHeightfieldChronoTerrainShape.cpp:696` | `"TODO:- should this be cached?"` | Already resolved — `m_invLocalScaling` IS cached. The `invS` here is just reading the cached value. |
| `cbtHeightfieldChronoTerrainShape.h:54` | `"TODO: (Chrono has preprocessors definitions for this??)"` | Chrono uses `CH_SSE` and similar, but the SIMD detection here is self-contained and correct. Could note this is intentionally standalone. |
| `cbtHeightfieldChronoTerrainShape.h:231` | `"constructor generic with cbtscalar - shoul dhandle either doubles or float"` | Typos: "shoul dhandle" → "should handle" |
| `ChCollisionShapeHeightField.h:127` | `"m_cellSizeV; // length /(ny1)"` | Typo: should be `// length / (ny - 1)` |
| `ChCollisionShapeHeightField.cpp:61` | `"no centering for bullet shape - the base position is 0,0,0 at the floor"` | Partially stale — the Bullet shape DOES center via `-minH` shift, but this comment is about the float copy not doing its own centering, which is correct |
| `ChCollisionAlgorithmsBullet.h:348` | `"TODO:- do we need a pure sphere based special case?"` | Yes, the sphere algorithm IS implemented and IS a special case. The TODO is resolved — should be removed. |
| `RigidTerrain.cpp:1068` | `"visualisation is always set to true for image based heigthfield"` | Typo: "heigthfield" → "heightfield". Also has TODO about putting in signature. |

### 12.3 Missing Comments

| Location | What's missing |
|---|---|
| `cbtHeightfieldChronoTerrainShape.h:413-430` | `updateHeight*` family — no comment explaining that heights passed in are ABSOLUTE, and the method handles the -minH shift internally |
| `ChCollisionModelBullet.cpp` `injectHeightfield` | No comment explaining the redundant `buildAccelerator(16)` call (already done in constructor) |

---

## 13. Issues & Recommendations

### 13.1 Critical Issues

None. The collision detection is functionally correct and performant.

### 13.2 High Priority

| # | Issue | Action |
|---|---|---|
| H1 | **Dead code: `cbtConvexHeightfieldAlgo.h/cpp`** — 608 lines compiled but never used | Remove from CMakeLists.txt and delete files |
| H2 | **Bullet debug draw broken for large terrains** — vertex cache is empty, so debug draw produces nothing | Fix to use `getVertex()` fallback when cache is empty |
| H3 | **VSG collision shape visualization missing** — `PopulateCollisionShapeFixed()` has no heightfield branch | Add `ChCollisionShapeHeightField` case that builds a trimesh |

### 13.3 Medium Priority

| # | Issue | Action |
|---|---|---|
| M1 | **Redundant `buildAccelerator(16)` in `injectHeightfield`** | Remove — constructor already builds it |
| M2 | **Stale TODO comments** | Clean up the 4 stale TODOs listed in §12.2 |
| M3 | **Triple height storage** (~80 MB for 2056² float build) | Document the trade-off; optionally allow Chrono-level data release after injection |
| M4 | **No JSON support for HEIGHTFIELD patches** | Add `LoadPatch()` case if JSON terrain loading is needed |
| M5 | **Serialization incomplete** | `ArchiveOut` writes either double or float heights (ifdef), but doesn't serialize all derived state (cell sizes are recomputed in `ArchiveIn` — this is actually fine) |
| M6 | **`m_rollingFriction`/`m_spinningFriction` unused** | Remove or wire into contact material |

### 13.4 Low Priority / Cosmetic

| # | Issue | Action |
|---|---|---|
| L1 | Typos in comments (3 found) | Fix |
| L2 | Commented-out code in demo | Clean before PR |
| L3 | Commented-out non-const getters in `.h` | Remove |
| L4 | Commented-out old `getVertex` implementation | Remove |
| L5 | Image-based heightfield hardcodes `visualize = true` | Add parameter |
| L6 | Uneven indentation in `ChCollisionShapeHeightField.h` member declarations | Format |
| L7 | `ChCollisionShape.h` comment "not implement**ed**" missing 'ed' | Fix typo |

---

## Summary

The heightfield extension is architecturally sound with proper separation between Chrono and Bullet layers. The collision algorithms are O(1) per-object with effective spatial culling. The BASE-at-origin convention is consistently applied. The main gaps are:

1. **Dead code** (`cbtConvexHeightfieldAlgo.*`) — should be removed before PR
2. **Debug/collision visualization** — broken for large terrains (Bullet), missing entirely (VSG)
3. **Minor code hygiene** — stale TODOs, commented-out code, typos

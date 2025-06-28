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

#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>

#include "cbtConvexHeightfieldAlgo.h"
#include "BulletCollision/CollisionDispatch/cbtCollisionDispatcher.h"
#include "LinearMath/cbtVector3.h"

#include "BulletCollision/CollisionShapes/cbtSphereShape.h"
#include "BulletCollision/CollisionShapes/cbtBoxShape.h"

using chrono_hf = cbtHeightfieldChronoTerrainShape;

// ---------------------------------------------------------------------
// ctor / dtor
// ---------------------------------------------------------------------
cbtConvexHeightfieldAlgo::cbtConvexHeightfieldAlgo(const cbtCollisionAlgorithmConstructionInfo& ci,
                                                   const cbtCollisionObjectWrapper* bodyA,
                                                   const cbtCollisionObjectWrapper* bodyB)
    : cbtCollisionAlgorithm(ci) {
    m_manifold = ci.m_dispatcher1->getNewManifold(bodyA->getCollisionObject(), bodyB->getCollisionObject());
}
cbtConvexHeightfieldAlgo::~cbtConvexHeightfieldAlgo() {
    m_dispatcher->releaseManifold(m_manifold);
}

static inline void gridCoords(const cbtVector3& p, int up, cbtScalar& u, cbtScalar& v) {
    switch (up) {  // return the two “planar” axes in HF order
        case 0:
            u = p.y();
            v = p.z();
            break;  // Xup, so grid=Y,Z
        case 1:
            u = p.x();
            v = p.z();
            break;  // Y up is grid=X,Z
        default:
            u = p.x();
            v = p.y();
            break;  // Z up is grid=X,Y
    }
}

// sphere handling is sorted in the processcollsion
static void collectSupportPoints(const cbtConvexShape* s, cbtAlignedObjectArray<cbtVector3>& out) {

    if (s->getShapeType() == BOX_SHAPE_PROXYTYPE) {
        const auto* box = static_cast<const cbtBoxShape*>(s);
        const cbtVector3 e = box->getHalfExtentsWithMargin();  // WITH margin
        out.resize(8);
        int i = 0;
        for (int x = -1; x <= 1; x += 2)
            for (int y = -1; y <= 1; y += 2)
                for (int z = -1; z <= 1; z += 2)
                    out[i++].setValue(x * e.x(), y * e.y(), z * e.z());
        return;  // nothing else to do for a box
    }

    // other shapes
    if (s->isPolyhedral()) {
        const auto* poly = static_cast<const cbtPolyhedralConvexShape*>(s);
        int nv = poly->getNumVertices();
        constexpr int kMaxDenseVerts = 128;
        int n = std::min(nv, kMaxDenseVerts);
        out.resize(n);
        for (int i = 0; i < n; ++i)
            poly->getVertex(i, out[i]);
        return;
    }

    // Sphere handling of ring-based contact points
    // TODO: alter this to skip non-likely points and save compute time
    if (s->getShapeType() == SPHERE_SHAPE_PROXYTYPE) {
        const cbtScalar R = static_cast<const cbtSphereShape*>(s)->getRadius();

        const cbtScalar deg2rad = SIMD_RADS_PER_DEG;
        //  15deg ring
        constexpr cbtScalar kSmall = 15.f;                 
        const cbtScalar sin15 = cbtSin(kSmall * deg2rad);  
        const cbtScalar cos15 = cbtCos(kSmall * deg2rad);

        // 45 deg ring
        const cbtScalar sin45 = SIMDSQRT12;
        const cbtScalar cos45 = SIMDSQRT12;

        out.resize(16);  // 8 + 8
        int idx = 0;
        auto push = [&](cbtScalar x, cbtScalar y, cbtScalar z) { out[idx++].setValue(x, y, z); };

        // 0‑7: 15deg ring (every 45deg)
        const cbtScalar r15 = R * sin15;
        const cbtScalar z15 = -R * cos15;
        for (int i = 0; i < 8; ++i) {
            const cbtScalar a = i * SIMD_PI / 4;  // 0 … 315deg
            push(r15 * cbtCos(a), r15 * cbtSin(a), z15);
        }

        // 8‑15: 45deg ring (every 45)
        const cbtScalar r45 = R * sin45;
        const cbtScalar z45 = -R * cos45;
        for (int i = 0; i < 8; ++i) {
            const cbtScalar a = i * SIMD_PI / 4;
            push(r45 * cbtCos(a), r45 * cbtSin(a), z45);
        }
        return;
    }

    // fallback – 26 direction GJK support map sample for other shapes
    static const cbtVector3 dir[26] = {{1, 0, 0},   {-1, 0, 0},  {0, 1, 0},  {0, -1, 0},  {0, 0, 1},  {0, 0, -1},
                                       {1, 1, 0},   {-1, 1, 0},  {1, -1, 0}, {-1, -1, 0}, {1, 0, 1},  {-1, 0, 1},
                                       {1, 0, -1},  {-1, 0, -1}, {0, 1, 1},  {0, -1, 1},  {0, 1, -1}, {0, -1, -1},
                                       {1, 1, 1},   {-1, 1, 1},  {1, -1, 1}, {-1, -1, 1}, {1, 1, -1}, {-1, 1, -1},
                                       {1, -1, -1}, {-1, -1, -1}};
    out.resize(26);
    for (int i = 0; i < 26; ++i)
        out[i] = s->localGetSupportingVertexWithoutMargin(dir[i]);
}


/// TODO: - a lot of the calcs in here need to be shifted to their own shapes processcollision to keep this cleaner

void cbtConvexHeightfieldAlgo::processCollision(const cbtCollisionObjectWrapper* bodyA,
                                                const cbtCollisionObjectWrapper* bodyB,
                                                const cbtDispatcherInfo&,
                                                cbtManifoldResult* result) {
    // Identify which body is the Chrono height-field
    const bool terrainIsA = bodyA->getCollisionShape()->getShapeType() == TERRAIN_SHAPE_PROXYTYPE;
    const auto* terrain = terrainIsA ? static_cast<const chrono_hf*>(bodyA->getCollisionShape())
                                     : static_cast<const chrono_hf*>(bodyB->getCollisionShape());
    const auto* shape = terrainIsA ? bodyB->getCollisionShape() : bodyA->getCollisionShape();

    const cbtTransform& trT = terrainIsA ? bodyA->getWorldTransform() : bodyB->getWorldTransform();
    const cbtTransform& trS = terrainIsA ? bodyB->getWorldTransform() : bodyA->getWorldTransform();

    result->setPersistentManifold(m_manifold);
    m_manifold->setBodies(bodyA->getCollisionObject(), bodyB->getCollisionObject());

    const int upAxis = terrain->getUpAxis();
    cbtVector3 upVec(0, 0, 0);
    upVec[upAxis] = 1;


    const cbtVector3 invScale = terrain->getInverseLocalScaling();
    const cbtVector3 localOrg = terrain->getLocalOrigin();
    const cbtScalar tMargin = terrain->getMargin();

    // ─────────────────────────────────────────────────────────────────────────────
    // Sphere vs Height-field – one deepest contact then form a tripod ring around it
    // ─────────────────────────────────────────────────────────────────────────────
    if (shape->getShapeType() == SPHERE_SHAPE_PROXYTYPE) {
        // constants
        const auto* sphere = static_cast<const cbtSphereShape*>(shape);
        const cbtScalar R = sphere->getRadius();  // radius *with* margin
        const cbtVector3 Cw = trS.getOrigin();    // sphere centre (world)
        // tolerance
        const cbtScalar tol = m_manifold->getContactBreakingThreshold() + result->m_closestPointDistanceThreshold;

        // best‐contact initiasations
        cbtScalar bestDepth = 1e30f;
        cbtVector3 bestN, bestP;

        // Build the 13 probe points in sphere‑local space
        cbtAlignedObjectArray<cbtVector3> probes;
        collectSupportPoints(static_cast<const cbtConvexShape*>(shape), probes);

        // support manifold
        struct Cand {
            cbtScalar d;
            cbtVector3 n, p;
        };
        Cand tripod[3] = {{1e30f}, {1e30f}, {1e30f}};  // three best inner‑ring hits

        // convert to heightfield upaxis helper function
        // switches to ensure we get the underside of of the sphere points fromo a default Z‑up offset to upAxis = 0/1/2
        // TODO - we should use the grid helper function above
        auto toHF = [&](const cbtVector3& o) {
            switch (upAxis) {
                case 0:
                    return cbtVector3(o.z(), o.x(), o.y());  // X up
                case 1:
                    return cbtVector3(o.x(), o.z(), o.y());  // Y up
                default:
                    return o;  // Z up
            }
        };

        //  Loop: query height‑field, radial normal and get the deepest contact point
        for (int i = 0; i < probes.size(); ++i) {
            // world probe = centre + (HF basis * permuted offset)
            cbtVector3 Sw = Cw + trT.getBasis() * toHF(probes[i]);

            // sample terrain in world space
            cbtVector3 PwSurf, nWsurf;
            if (!terrain->sampleWorld(trT, Sw, PwSurf, nWsurf))
                continue;

            // slide vertical point onto radial line
            cbtVector3 dirW = Sw - Cw;
            cbtScalar len = dirW.length();
            if (len < SIMD_EPSILON)
                continue;
            dirW /= len;
            cbtScalar t = (PwSurf - Sw).dot(nWsurf) / dirW.dot(nWsurf);
            PwSurf = Sw + dirW * t;

            // radial normal so that its not generating phantom torque
            cbtVector3 radial = Cw - PwSurf;
            if (radial.fuzzyZero())
                continue;
            cbtVector3 nW = radial.normalized();
            if (terrainIsA)
                nW = -nW;

            // penetration depth calc
            cbtScalar depth = radial.length() - R;
            if (depth > tol)
                continue;

            // keep the deepest
            if (depth < bestDepth) {
                bestDepth = depth;
                bestN = nW;
                bestP = terrainIsA ? (Cw + nW * R) : PwSurf;
            }

            // inner‑ring candidates are indices 1‑8
            if (i >= 1 && i <= 8) {
                // insert into tripod[0..2] if it is among the three deepest so far
                for (int s = 0; s < 3; ++s)
                    if (depth < tripod[s].d) {
                        for (int t = 2; t > s; --t)
                            tripod[t] = tripod[t - 1];  // shift down
                        tripod[s] = {depth, nW, terrainIsA ? (Cw + nW * R) : PwSurf};
                        break;
                    }
            }
        }

        // setup for a plane to build a small tripod around the best point
        // this stops he sphere from rolling unaturally/eternally by handing control
        // effectively over to the bullet margin handler process to deal with rather than
        // only a single point contact which can cause neverending movement
        // emit contacts (4 in total, at most)
        if (bestDepth < 1e30f) {
            result->addContactPoint(bestN, bestP, bestDepth);
            for (int s = 0; s < 3; ++s)
                if (tripod[s].d < 1e30f)  // valid?
                    result->addContactPoint(tripod[s].n, tripod[s].p, tripod[s].d);

            result->refreshContactPoints();
        }
        return;  // sphere handled, leave function early
    }


    //-------------------------------------------------------------
    // Generic convex against height-field (sampleWorld + correct transforms)
    //-------------------------------------------------------------
    static thread_local cbtAlignedObjectArray<cbtVector3> verts;
    verts.clear();
    collectSupportPoints(static_cast<const cbtConvexShape*>(shape), verts);

    // margins and vertical band
    const cbtScalar cMargin = static_cast<const cbtConvexShape*>(shape)->getMargin();
    const cbtScalar gap = cMargin + tMargin;

    // compute vMin/vMax and 2D grid bounds from support points
    cbtScalar vMin = 1e30f, vMax = -1e30f;
    cbtScalar gxMin = 1e30f, gxMax = -1e30f;
    cbtScalar gzMin = 1e30f, gzMax = -1e30f;
    for (int i = 0; i < verts.size(); ++i) {
        const cbtVector3& vLocal = verts[i];
        // transform the box‐local point into world‐space
        cbtVector3 Pw = trS * vLocal;
        // then into heightfield local cell coords
        cbtVector3 Pl = (trT.invXform(Pw) * invScale) + localOrg;

        vMin = cbtMin(vMin, Pl[upAxis]);
        vMax = cbtMax(vMax, Pl[upAxis]);

        cbtScalar gx, gz;
        gridCoords(Pl, upAxis, gx, gz);
        gxMin = cbtMin(gxMin, gx);
        gxMax = cbtMax(gxMax, gx);
        gzMin = cbtMin(gzMin, gz);
        gzMax = cbtMax(gzMax, gz);
    }
    vMin -= gap;
    vMax += gap;

    // convert to chunk indices
    constexpr int kChunkShift = 4;
    int cx0 = cbtMax(0, int(gxMin) >> kChunkShift);
    int cx1 = cbtMin((terrain->getWidth() - 2) >> kChunkShift, int(gxMax) >> kChunkShift);
    int cz0 = cbtMax(0, int(gzMin) >> kChunkShift);
    int cz1 = cbtMin((terrain->getLength() - 2) >> kChunkShift, int(gzMax) >> kChunkShift);

    // Bullet’s tolerance
    const cbtScalar tol = m_manifold->getContactBreakingThreshold() + result->m_closestPointDistanceThreshold;

    // per-vertex narrow-phase
    for (int i = 0; i < verts.size(); ++i) {
        const cbtVector3& vLocal = verts[i];
        // box‐space → world‐space
        cbtVector3 Pw = trS * vLocal;
        // world → heightfield‐local
        cbtVector3 Pl = (trT.invXform(Pw) * invScale) + localOrg;

        // 2D grid cull
        cbtScalar gx, gz;
        gridCoords(Pl, upAxis, gx, gz);
        if (gx < gxMin || gx > gxMax || gz < gzMin || gz > gzMax)
            continue;

        // chunk‐box cull
        int cx = int(gx) >> kChunkShift;
        int cz = int(gz) >> kChunkShift;
        if (cx < cx0 || cx > cx1 || cz < cz0 || cz > cz1)
            continue;

        // height‐range cull
        const chrono_hf::Range& R = terrain->GetVBoundsChunk(cx, cz);
        if (R.max < vMin - tol || R.min > vMax + tol)
            continue;

        // robust world sampling
        cbtVector3 surfP, surfN;
        if (!terrain->sampleWorld(trT, Pw, surfP, surfN))
            continue;

        // slide Pw onto the surface along surfN
        cbtScalar slide = (surfP - Pw).dot(surfN);
        cbtVector3 PwSurf = Pw + surfN * slide;

        // compute penetration depth
        cbtScalar pen = -slide - cMargin;
        if (pen > tol)
            continue;

        // add the contact
        cbtVector3 nW = terrainIsA ? -surfN : surfN;
        result->addContactPoint(nW, PwSurf, pen);

    }

    // finalize manifold
    result->refreshContactPoints();
    return;
}
       

cbtScalar cbtConvexHeightfieldAlgo::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                          cbtCollisionObject* body1,
                                                          const cbtDispatcherInfo&,
                                                          cbtManifoldResult*) {
    const bool terrainIsA = (body0->getCollisionShape()->getShapeType() == TERRAIN_SHAPE_PROXYTYPE);

    const auto* terrain = terrainIsA ? static_cast<const chrono_hf*>(body0->getCollisionShape())
                                     : static_cast<const chrono_hf*>(body1->getCollisionShape());
    const auto* convex = terrainIsA ? static_cast<const cbtConvexShape*>(body1->getCollisionShape())
                                    : static_cast<const cbtConvexShape*>(body0->getCollisionShape());

    const cbtTransform& trFrom = terrainIsA ? body1->getWorldTransform() : body0->getWorldTransform();
    const cbtTransform& trTo =
        terrainIsA ? body1->getInterpolationWorldTransform() : body0->getInterpolationWorldTransform();
    const cbtTransform& trTerrain = terrainIsA ? body0->getWorldTransform() : body1->getWorldTransform();

    const cbtVector3 motion = trTo.getOrigin() - trFrom.getOrigin();
    if (motion.length2() < SIMD_EPSILON)
        return 1.f;

    const int up = terrain->getUpAxis();
    const cbtScalar totalMargin = convex->getMargin() + terrain->getMargin();
    const cbtVector3 invScale(1.f / terrain->getLocalScaling().getX(), 1.f / terrain->getLocalScaling().getY(),
                              1.f / terrain->getLocalScaling().getZ());
    const cbtVector3 localOrig = terrain->getLocalOrigin();

    auto grid = [&](const cbtVector3& P, cbtScalar& gx, cbtScalar& gz) {
        switch (up) {
            case 0:
                gx = P.getY();
                gz = P.getZ();
                break;
            case 1:
                gx = P.getX();
                gz = P.getZ();
                break;
            case 2:
                gx = P.getX();
                gz = P.getY();
                break;
        }
    };

    /* simple binary subdivision along linear path of centre point */
    cbtScalar lo = 0.f, hi = 1.f;
    for (int it = 0; it < 25; ++it) {
        const cbtScalar mid = 0.5f * (lo + hi);
        const cbtVector3 Pw = trFrom.getOrigin() + motion * mid;
        const cbtVector3 Pl = (trTerrain.invXform(Pw) * invScale) + localOrig;

        cbtScalar gx, gz;
        grid(Pl, gx, gz);
        bool hit = false;
        if (gx >= 0 && gz >= 0 && gx < terrain->getWidth() - 1 && gz < terrain->getLength() - 1) {
            cbtScalar height;
            cbtVector3 grad;
            terrain->sampleHeight(gx, gz, height, grad);
            cbtScalar radius = 0.f;
            if (convex->getShapeType() == SPHERE_SHAPE_PROXYTYPE)
                radius = static_cast<const cbtSphereShape*>(convex)->getRadius();
            const cbtScalar distSurface = (Pl[up] - localOrig[up])  // centre to plane
                                              -(height - localOrig[up])  // minus terrain height
                                              -radius;                // subtract sphere radius
            hit = distSurface <= 0;
        }
        (hit ? hi : lo) = mid;
    }
    return hi;  // fraction in [0,1]
}

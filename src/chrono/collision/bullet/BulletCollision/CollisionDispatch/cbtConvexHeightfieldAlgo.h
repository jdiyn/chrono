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
//  Triangle‑free convex–heightfield narrow‑phase for Chrono custom bullet
// =============================================================================

#include "BulletCollision/BroadphaseCollision/cbtCollisionAlgorithm.h"
#include "BulletCollision/CollisionShapes/cbtConvexShape.h"
#include "BulletCollision/CollisionShapes/cbtPolyhedralConvexShape.h"
#include "BulletCollision/NarrowPhaseCollision/cbtPersistentManifold.h"
#include "BulletCollision/CollisionShapes/cbtHeightfieldChronoTerrainShape.h"
#include "BulletCollision/CollisionDispatch/cbtCollisionDispatcher.h" 
#include "BulletCollision/CollisionDispatch/cbtCollisionCreateFunc.h"


class cbtConvexHeightfieldAlgo final : public cbtCollisionAlgorithm {
    cbtPersistentManifold* m_manifold;

  public:
    cbtConvexHeightfieldAlgo(const cbtCollisionAlgorithmConstructionInfo& ci,
                             const cbtCollisionObjectWrapper* bodyA,
                             const cbtCollisionObjectWrapper* bodyB);

    ~cbtConvexHeightfieldAlgo() override;

    // -----------------------------------------------------------------
    // discrete contacts
    // -----------------------------------------------------------------
    void processCollision(const cbtCollisionObjectWrapper* bodyA,
                          const cbtCollisionObjectWrapper* bodyB,
                          const cbtDispatcherInfo& info,
                          cbtManifoldResult* result) override;

    // -----------------------------------------------------------------
    // conservative centre‑point sweep (fast CCD)
    // -----------------------------------------------------------------
    cbtScalar calculateTimeOfImpact(cbtCollisionObject* body0,
                                    cbtCollisionObject* body1,
                                    const cbtDispatcherInfo& info,
                                    cbtManifoldResult*) override;

    // manifold access
    void getAllContactManifolds(cbtManifoldArray& array) override { array.push_back(m_manifold); }

        // --- factory helpers expected by the dispatcher loop --------------
    struct CreateFunc : public cbtCollisionAlgorithmCreateFunc {
        cbtCollisionAlgorithm* CreateCollisionAlgorithm(cbtCollisionAlgorithmConstructionInfo& ci,
                                                        const cbtCollisionObjectWrapper* body0,
                                                        const cbtCollisionObjectWrapper* body1) override {
            void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtConvexHeightfieldAlgo));
            return new (mem) cbtConvexHeightfieldAlgo(ci, body0, body1);
        }
    };
    struct SwappedCreateFunc : public cbtCollisionAlgorithmCreateFunc {
        cbtCollisionAlgorithm* CreateCollisionAlgorithm(cbtCollisionAlgorithmConstructionInfo& ci,
                                                        const cbtCollisionObjectWrapper* body0,
                                                        const cbtCollisionObjectWrapper* body1) override {
            void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtConvexHeightfieldAlgo));
            // flip so (convex, terrain) order is preserved
            return new (mem) cbtConvexHeightfieldAlgo(ci, body1, body0);
        }
    };
};

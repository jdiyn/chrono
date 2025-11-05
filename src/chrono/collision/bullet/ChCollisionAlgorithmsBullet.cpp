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
// Authors: Alessandro Tasora, Radu Serban
// =============================================================================

#include "chrono/collision/bullet/ChCollisionAlgorithmsBullet.h"
#include "chrono/collision/bullet/ChCollisionModelBullet.h"
#include "chrono/collision/bullet/ChCollisionUtilsBullet.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtils.h"

#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtSphereShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtCylinderShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtBoxShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtCylindricalShellShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtCapsuleShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbt2DShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtCEtriangleShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtSegmentShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionShapes/cbtHeightfieldChronoTerrainShape.h"
#include "chrono/collision/bullet/BulletCollision/CollisionDispatch/SphereTriangleDetector.h"

namespace chrono {

// ================================================================================================

cbtCapsuleBoxCollisionAlgorithm::cbtCapsuleBoxCollisionAlgorithm(cbtPersistentManifold* mf,
                                                                 const cbtCollisionAlgorithmConstructionInfo& ci,
                                                                 const cbtCollisionObjectWrapper* col0,
                                                                 const cbtCollisionObjectWrapper* col1,
                                                                 bool isSwapped)
    : cbtActivatingCollisionAlgorithm(ci, col0, col1), m_ownManifold(false), m_manifoldPtr(mf), m_isSwapped(isSwapped) {
    const cbtCollisionObjectWrapper* capsuleObjWrap = m_isSwapped ? col1 : col0;
    const cbtCollisionObjectWrapper* boxObjWrap = m_isSwapped ? col0 : col1;

    if (!m_manifoldPtr &&
        m_dispatcher->needsCollision(capsuleObjWrap->getCollisionObject(), boxObjWrap->getCollisionObject())) {
        m_manifoldPtr =
            m_dispatcher->getNewManifold(capsuleObjWrap->getCollisionObject(), boxObjWrap->getCollisionObject());
        m_ownManifold = true;
    }
}

cbtCapsuleBoxCollisionAlgorithm::cbtCapsuleBoxCollisionAlgorithm(const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtCapsuleBoxCollisionAlgorithm ::~cbtCapsuleBoxCollisionAlgorithm() {
    if (m_ownManifold) {
        if (m_manifoldPtr)
            m_dispatcher->releaseManifold(m_manifoldPtr);
    }
}

// Capsule-box intersection test.
void cbtCapsuleBoxCollisionAlgorithm::processCollision(const cbtCollisionObjectWrapper* body0,
                                                       const cbtCollisionObjectWrapper* body1,
                                                       const cbtDispatcherInfo& dispatchInfo,
                                                       cbtManifoldResult* resultOut) {
    (void)dispatchInfo;
    (void)resultOut;
    if (!m_manifoldPtr)
        return;

    const cbtCollisionObjectWrapper* capObjWrap = m_isSwapped ? body1 : body0;
    const cbtCollisionObjectWrapper* boxObjWrap = m_isSwapped ? body0 : body1;

    resultOut->setPersistentManifold(m_manifoldPtr);

    const cbtCapsuleShape* cap = (cbtCapsuleShape*)capObjWrap->getCollisionShape();
    const cbtBoxShape* box = (cbtBoxShape*)boxObjWrap->getCollisionShape();

    // Express capsule in the box frame
    const cbtTransform& abs_X_cap = capObjWrap->getWorldTransform();
    const cbtTransform& abs_X_box = boxObjWrap->getWorldTransform();
    cbtTransform box_X_cap = abs_X_box.inverseTimes(abs_X_cap);

    cbtVector3 a = box_X_cap.getBasis().getColumn(1);  // capsule axis (expressed in box frame)
    cbtVector3 c = box_X_cap.getOrigin();              // capsule center (expressed in box frame)

    // Box dimensions
    cbtVector3 hdims = box->getHalfExtentsWithMargin();

    // Cylinder dimensions
    cbtScalar radius = cap->getRadius();    // capsule radius
    cbtScalar hlen = cap->getHalfHeight();  // capsule half-length

    // Loop over each direction of the box frame (i.e., each of the 3 face normals).
    // In each case, consider two segments on the cylindrical surface that are on a plane defined by the axis and the
    // face normal. (Note that, in principle, we could only consider the segment "closest" to the box, but that is not
    // trivial to define in all configurations). Such segments are parameterized by t in [-H,H].
    //
    // Consider (at most) 4 candidate points on each segment: the 2 ends of the segment and the intersections of the
    // segment with the box (if such an intersection exists).
    //
    // For a capsule, projects these points back onto the capsule axis and check intersection between spheres centered
    // at the points on the axis and the box.  In this case, the projection is done orthogonal to the capsule axis.

    const cbtScalar parallel_tol = cbtScalar(1e-5);           // tolearance for parallelism tests
    const cbtScalar near_tol = cbtScalar(1e-4) * (2 * hlen);  // tolerance for line parameters of near duplicate points

    std::vector<cbtScalar> t_points = {-hlen, +hlen};  // initialize list of candidates with segment ends

    for (int i = 0; i < 3; i++) {
        // "positive" face normal
        cbtVector3 n(0, 0, 0);
        n[i] = 1;

        // If the axis is parallel to the face normal, no additional candidate point.
        cbtVector3 v = n.cross(a);
        if (std::abs(a[i] - 1) < parallel_tol)
            continue;

        // Direction perpendicular to axis
        cbtVector3 r = v.cross(a);

        // Construct center points of the two segments on cylindrical surface
        cbtVector3 c1 = c + radius * r;
        cbtVector3 c2 = c - radius * r;

        // Check if either segment intersects box.
        // If it does, append line parameters for intersection points (clamped to segment limits).
        cbtScalar tMin;
        cbtScalar tMax;
        if (bt_utils::IntersectSegmentBox(hdims, c1, a, hlen, parallel_tol, tMin, tMax)) {
            t_points.push_back(ChClamp(tMin, -hlen, +hlen));
            t_points.push_back(ChClamp(tMax, -hlen, +hlen));
        }
        if (bt_utils::IntersectSegmentBox(hdims, c2, a, hlen, parallel_tol, tMin, tMax)) {
            t_points.push_back(ChClamp(tMin, -hlen, +hlen));
            t_points.push_back(ChClamp(tMax, -hlen, +hlen));
        }
    }

    // Contact distance
    cbtScalar contactDist = radius;
    ////cbtScalar contactDist = radius + m_manifoldPtr->getContactBreakingThreshold();

    // Loop over all candidate points (points on the capsule axis) and perform a sphere-box collision test.
    // In order to eliminate near duplicate points, use a sorted list and keep track of last processed point.
    std::sort(t_points.begin(), t_points.end());
    cbtScalar t_last = -2 * hlen;
    int n_contacts = 0;
    for (auto t : t_points) {
        if (t - t_last < near_tol)
            continue;

        // Update last processed point
        t_last = t;

        // Calculate sphere center (point on axis, expressed in box frame) and snap it to the surface of the box
        cbtVector3 sphPos = c + a * t;
        cbtVector3 boxPos = sphPos;
        bt_utils::SnapPointToBox(hdims, boxPos);

        // If the distance from the sphere center to the closest point is larger than the radius plus the separation
        // value, then there is no contact. Also, ignore contact if the sphere center (almost) coincides with the
        // closest point, in which case we couldn't decide on the proper contact direction.
        cbtVector3 delta = sphPos - boxPos;
        cbtScalar dist2 = delta.length2();
        if (dist2 >= contactDist * contactDist || dist2 <= 1e-12)
            continue;

        // Generate contact information (transform to absolute frame)
        cbtScalar dist = cbtSqrt(dist2);
        cbtScalar penetration = dist - radius;
        cbtVector3 normal = abs_X_box.getBasis() * (delta / dist);
        cbtVector3 point = abs_X_box(boxPos);

        // A new contact point must specify:
        //   normal, pointing from B towards A
        //   point, located on surface of B
        //   distance, negative for penetration
        resultOut->addContactPoint(normal, point, penetration);
        n_contacts++;
    }

    if (m_ownManifold && m_manifoldPtr->getNumContacts()) {
        resultOut->refreshContactPoints();
    }
}

cbtScalar cbtCapsuleBoxCollisionAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                                 cbtCollisionObject* body1,
                                                                 const cbtDispatcherInfo& dispatchInfo,
                                                                 cbtManifoldResult* resultOut) {
    // not yet
    return cbtScalar(1.);
}

void cbtCapsuleBoxCollisionAlgorithm::getAllContactManifolds(cbtManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
        manifoldArray.push_back(m_manifoldPtr);
    }
}

cbtCollisionAlgorithm* cbtCapsuleBoxCollisionAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* body0Wrap,
    const cbtCollisionObjectWrapper* body1Wrap) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtCapsuleBoxCollisionAlgorithm));
    if (!m_swapped) {
        return new (mem) cbtCapsuleBoxCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
    } else {
        return new (mem) cbtCapsuleBoxCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
    }
}

// ================================================================================================

cbtCylshellBoxCollisionAlgorithm::cbtCylshellBoxCollisionAlgorithm(cbtPersistentManifold* mf,
                                                                   const cbtCollisionAlgorithmConstructionInfo& ci,
                                                                   const cbtCollisionObjectWrapper* col0,
                                                                   const cbtCollisionObjectWrapper* col1,
                                                                   bool isSwapped)
    : cbtActivatingCollisionAlgorithm(ci, col0, col1), m_ownManifold(false), m_manifoldPtr(mf), m_isSwapped(isSwapped) {
    const cbtCollisionObjectWrapper* cylshellObjWrap = m_isSwapped ? col1 : col0;
    const cbtCollisionObjectWrapper* boxObjWrap = m_isSwapped ? col0 : col1;

    if (!m_manifoldPtr &&
        m_dispatcher->needsCollision(cylshellObjWrap->getCollisionObject(), boxObjWrap->getCollisionObject())) {
        m_manifoldPtr =
            m_dispatcher->getNewManifold(cylshellObjWrap->getCollisionObject(), boxObjWrap->getCollisionObject());
        m_ownManifold = true;
    }
}

cbtCylshellBoxCollisionAlgorithm::cbtCylshellBoxCollisionAlgorithm(const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtCylshellBoxCollisionAlgorithm::~cbtCylshellBoxCollisionAlgorithm() {
    if (m_ownManifold) {
        if (m_manifoldPtr)
            m_dispatcher->releaseManifold(m_manifoldPtr);
    }
}

// Check and add contact between the given cylshell point and the specified box face.
// 'iface' is +1, +2, +3 for the "positive" x, y, or z box face, respectively.
// 'iface' is -1, -2, -3 for the "negative" x, y, or z box face, respectively.
int addContactPoint(const cbtVector3& pc,
                    int iface,
                    const cbtVector3& hdims,
                    const cbtTransform& X_box,
                    cbtManifoldResult* resultOut) {
    assert(iface >= -3 && iface <= +3 && iface != 0);

    // No contact if point outside box
    if (!bt_utils::PointInsideBox(hdims, pc))
        return 0;  // no contacts added

    // Find point projection on box face and calculate normal and penetration
    // (still working in the box frame)
    cbtVector3 p = pc;
    cbtVector3 n(0, 0, 0);
    cbtScalar penetration;
    if (iface > 0) {
        // "positive" box face
        int i = iface - 1;
        p[i] = hdims[i];
        n[i] = 1;
        penetration = pc[i] - hdims[i];
    } else {
        // "negative" box face
        int i = -iface - 1;
        p[i] = -hdims[i];
        n[i] = -1;
        penetration = -pc[i] - hdims[i];
    }

    // A new contact point must specify (in absolute frame):
    //   normal, pointing from B towards A
    //   point, located on surface of B
    //   distance, negative for penetration
    cbtVector3 normal = X_box.getBasis() * n;
    cbtVector3 point = X_box(p);
    resultOut->addContactPoint(normal, point, penetration);

    ////std::cout << "add contact     nrm = " << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;
    ////std::cout << "                pnt = " << point.x() << " " << point.y() << " " << point.z() << std::endl;
    ////std::cout << "                pen = " << penetration << std::endl;

    return 1;  // one contact added
}

// Add contact between the given box point (assumed to be in or on the cylinder) and the cylshell.
// All input vectors are assumed to be expressed in the box frame.
int addContactPoint(const cbtVector3& p,
                    const cbtVector3& c,
                    const cbtVector3& a,
                    const cbtScalar h,
                    const cbtScalar r,
                    const cbtTransform& X_box,
                    cbtManifoldResult* resultOut) {
    // Find closest point on cylindrical surface to given location
    cbtVector3 q = bt_utils::ProjectPointOnLine(c, a, p);
    cbtVector3 v = p - q;
    cbtScalar dist = v.length();
    cbtVector3 n = v / dist;

    cbtVector3 normal = X_box.getBasis() * (-n);
    cbtVector3 point = X_box(p);

    resultOut->addContactPoint(normal, point, dist - r);

    return 1;
}

// Cylshell-box intersection test:
//   - cylinder caps are ignored
//   - the cylshell is replaced with a capsule on the surface of the cylshell
//   - capsule-box intersection is then reduced to a segment-box intersection
//   - a replacement capsule (one for each direction of the box) may generate 0, 1, or 2 contacts
void cbtCylshellBoxCollisionAlgorithm::processCollision(const cbtCollisionObjectWrapper* body0,
                                                        const cbtCollisionObjectWrapper* body1,
                                                        const cbtDispatcherInfo& dispatchInfo,
                                                        cbtManifoldResult* resultOut) {
    (void)dispatchInfo;
    (void)resultOut;
    if (!m_manifoldPtr)
        return;

    const cbtCollisionObjectWrapper* cylObjWrap = m_isSwapped ? body1 : body0;
    const cbtCollisionObjectWrapper* boxObjWrap = m_isSwapped ? body0 : body1;

    resultOut->setPersistentManifold(m_manifoldPtr);

    const cbtCylindricalShellShape* cyl = (cbtCylindricalShellShape*)cylObjWrap->getCollisionShape();
    const cbtBoxShape* box = (cbtBoxShape*)boxObjWrap->getCollisionShape();

    // Express cylinder in the box frame
    const cbtTransform& abs_X_cyl = cylObjWrap->getWorldTransform();
    const cbtTransform& abs_X_box = boxObjWrap->getWorldTransform();
    cbtTransform box_X_cyl = abs_X_box.inverseTimes(abs_X_cyl);

    cbtVector3 a = box_X_cyl.getBasis().getColumn(1);  // cylinder axis (expressed in box frame)
    cbtVector3 c = box_X_cyl.getOrigin();              // cylinder center (expressed in box frame)

    // Box dimensions
    cbtVector3 hdims = box->getHalfExtentsWithMargin();

    // Cylinder dimensions
    cbtScalar radius = cyl->getRadius();    // cylinder radius
    cbtScalar hlen = cyl->getHalfLength();  // cylinder half-length

    const cbtScalar parallel_tol = cbtScalar(1e-5);  // tolearance for parallelism tests

    int num_contacts = 0;

    // - Loop over each direction of the box frame (i.e., each of the 3 face normals).
    // - For each direction, consider two segments on the cylindrical surface that are on a plane defined by the axis
    //   and the face normal. (Note that, in principle, we could only consider the segment "closest" to the box, but
    //   that is not trivial to define in all configurations). All segments are parameterized by t in [-H,H].
    // - For each segment, if the segment intersects the box, consider 3 candidate contact points: the 2 intersection
    //   points and their midpoint. A contact is added if the segment point is inside the box.
    //   Furthermore, the corresponding box point is located on the box face that is closest to the intersection
    //   midpoint candidate.
    for (int idir = 0; idir < 3; idir++) {
        // current box direction
        cbtVector3 ndir(0, 0, 0);
        ndir[idir] = 1;

        // If the axis is parallel to the current direction, no contact.
        if (std::abs(a[idir] - 1) < parallel_tol || std::abs(a[idir] + 1) < parallel_tol)
            continue;

        // Direction perpendicular to cylinder axis (in direction opposite to ndir)
        cbtVector3 v = ndir.cross(a);
        cbtVector3 r = v.cross(a);
        assert(r.length() > parallel_tol);
        r.normalize();

        // Consider segments in both "negative" and "positive" r direction
        cbtScalar dir[2] = {-1, 1};
        for (int jdir = 0; jdir < 2; jdir++) {
            // Calculate current segment center
            cbtVector3 cs = c + dir[jdir] * radius * r;
            // Check for intersection with box
            cbtScalar tMin, tMax;
            if (bt_utils::IntersectSegmentBox(hdims, cs, a, hlen, parallel_tol, tMin, tMax)) {
                // Consider the intersection points and their midpoint as candidates
                cbtVector3 pMin = cs + a * tMin;
                cbtVector3 pMax = cs + a * tMax;
                cbtVector3 pMid = cs + a * ((tMin + tMax) / 2);

                // Pick box face that is closest to midpoint
                int iface = bt_utils::FindClosestBoxFace(hdims, pMid);

                // Add a contact for any of the candidate points that is inside the box
                num_contacts += addContactPoint(pMin, iface, hdims, abs_X_box, resultOut);  // 1st segment end
                num_contacts += addContactPoint(pMax, iface, hdims, abs_X_box, resultOut);  // 2nd segment end
                num_contacts += addContactPoint(pMid, iface, hdims, abs_X_box, resultOut);  // intersection midpoint
            }
        }
    }

    // If a box face supports the cylinder, do not check box edges.
    if (num_contacts > 0)
        return;

    // - Loop over each direction of the box frame.
    // - For each direction, check intersection with the cylinder for all 4 edges parallel to that direction.
    // - If an edge intersects the cylinder, consider 3 candidate contact points: the 2 intersection points
    //   and their midpoint.
    for (int idir = 0; idir < 3; idir++) {
        // current box edge direction and halflength
        cbtVector3 eD(0, 0, 0);
        eD[idir] = 1;
        cbtScalar eH = hdims[idir];
        // The other two box directions
        int jdir = (idir + 1) % 3;
        int kdir = (idir + 2) % 3;
        for (int j = -1; j <= +1; j += 2) {
            for (int k = -1; k <= +1; k += 2) {
                cbtVector3 eC;
                eC[idir] = 0;
                eC[jdir] = j * hdims[jdir];
                eC[kdir] = k * hdims[kdir];
                // Check for edge intersection with cylinder
                cbtScalar tMin, tMax;
                if (bt_utils::IntersectSegmentCylinder(eC, eD, eH, c, a, hlen, radius, parallel_tol, tMin, tMax)) {
                    // Consider the intersection points and their midpoint as candidates
                    cbtVector3 pMin = eC + eD * tMin;
                    cbtVector3 pMax = eC + eD * tMax;
                    cbtVector3 pMid = eC + eD * ((tMin + tMax) / 2);

                    // Add a contact for any of the candidate points that is inside the cylinder
                    num_contacts += addContactPoint(pMin, c, a, hlen, radius, abs_X_box, resultOut);
                    num_contacts += addContactPoint(pMax, c, a, hlen, radius, abs_X_box, resultOut);
                    num_contacts += addContactPoint(pMid, c, a, hlen, radius, abs_X_box, resultOut);
                }
            }
        }
    }

    ////std::cout << num_contacts << std::endl;

    if (m_ownManifold && m_manifoldPtr->getNumContacts()) {
        resultOut->refreshContactPoints();
    }
}

cbtScalar cbtCylshellBoxCollisionAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                                  cbtCollisionObject* body1,
                                                                  const cbtDispatcherInfo& dispatchInfo,
                                                                  cbtManifoldResult* resultOut) {
    // not yet
    return cbtScalar(1.);
}

void cbtCylshellBoxCollisionAlgorithm::getAllContactManifolds(cbtManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
        manifoldArray.push_back(m_manifoldPtr);
    }
}

cbtCollisionAlgorithm* cbtCylshellBoxCollisionAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* body0Wrap,
    const cbtCollisionObjectWrapper* body1Wrap) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtCylshellBoxCollisionAlgorithm));
    if (!m_swapped) {
        return new (mem) cbtCylshellBoxCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
    } else {
        return new (mem) cbtCylshellBoxCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
    }
}

// ================================================================================================

cbtSphereCylinderCollisionAlgorithm::cbtSphereCylinderCollisionAlgorithm(
    cbtPersistentManifold* mf,
    const cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* col0,
    const cbtCollisionObjectWrapper* col1,
    bool isSwapped)
    : cbtActivatingCollisionAlgorithm(ci, col0, col1), m_ownManifold(false), m_manifoldPtr(mf), m_isSwapped(isSwapped) {
    const cbtCollisionObjectWrapper* sphereObj = m_isSwapped ? col1 : col0;
    const cbtCollisionObjectWrapper* cylObj = m_isSwapped ? col0 : col1;

    if (!m_manifoldPtr && m_dispatcher->needsCollision(sphereObj->getCollisionObject(), cylObj->getCollisionObject())) {
        m_manifoldPtr = m_dispatcher->getNewManifold(sphereObj->getCollisionObject(), cylObj->getCollisionObject());
        m_ownManifold = true;
    }
}

cbtSphereCylinderCollisionAlgorithm::cbtSphereCylinderCollisionAlgorithm(
    const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtSphereCylinderCollisionAlgorithm ::~cbtSphereCylinderCollisionAlgorithm() {
    if (m_ownManifold) {
        if (m_manifoldPtr)
            m_dispatcher->releaseManifold(m_manifoldPtr);
    }
}

void cbtSphereCylinderCollisionAlgorithm::processCollision(const cbtCollisionObjectWrapper* body0,
                                                           const cbtCollisionObjectWrapper* body1,
                                                           const cbtDispatcherInfo& dispatchInfo,
                                                           cbtManifoldResult* resultOut) {
    (void)dispatchInfo;
    (void)resultOut;
    if (!m_manifoldPtr)
        return;

    const cbtCollisionObjectWrapper* sphereObjWrap = m_isSwapped ? body1 : body0;
    const cbtCollisionObjectWrapper* cylObjWrap = m_isSwapped ? body0 : body1;

    resultOut->setPersistentManifold(m_manifoldPtr);

    const cbtSphereShape* sphere0 = (cbtSphereShape*)sphereObjWrap->getCollisionShape();
    const cbtCylinderShape* cylinder = (cbtCylinderShape*)cylObjWrap->getCollisionShape();

    const cbtTransform& m44T = cylObjWrap->getCollisionObject()->getWorldTransform();
    cbtVector3 diff = m44T.invXform(
        sphereObjWrap->getCollisionObject()
            ->getWorldTransform()
            .getOrigin());  // col0->getWorldTransform().getOrigin()-  col1->getWorldTransform().getOrigin();
    cbtScalar radius0 = sphere0->getRadius();
    cbtScalar radius1 = cylinder->getHalfExtentsWithMargin().getX();  // cylinder->getRadius();
    cbtScalar H1 = cylinder->getHalfExtentsWithMargin().getY();

    cbtVector3 r1 = diff;
    r1.setY(0);

    cbtScalar y1 = diff.y();

    cbtScalar r1_len = r1.length();

    cbtVector3 pos1;
    cbtVector3 normalOnSurfaceB(1, 0, 0);
    cbtScalar dist;

    // Case A
    if ((y1 <= H1) && (y1 >= -H1)) {
        /// iff distance positive, don't generate a new contact
        if (r1_len > (radius0 + radius1)) {
            resultOut->refreshContactPoints();
            return;
        }
        /// distance (negative means penetration)
        dist = r1_len - (radius0 + radius1);

        cbtVector3 localnormalOnSurfaceB;
        if (r1_len > SIMD_EPSILON) {
            localnormalOnSurfaceB = r1 / r1_len;
            normalOnSurfaceB = m44T.getBasis() * localnormalOnSurfaceB;
        }
        /// point on B (worldspace)
        pos1 = m44T(cbtVector3(0, y1, 0)) + radius1 * normalOnSurfaceB;
    } else {
        cbtScalar side = 1;
        if (y1 < -H1)
            side = -1;

        if (r1_len > radius1) {
            // case B
            cbtVector3 pos_loc = r1.normalized() * radius1 + cbtVector3(0, H1 * side, 0);
            pos1 = m44T(pos_loc);
            cbtVector3 d = sphereObjWrap->getCollisionObject()->getWorldTransform().getOrigin() - pos1;
            normalOnSurfaceB = d.normalized();
            dist = d.length() - radius0;
        } else {
            // case C
            normalOnSurfaceB = m44T.getBasis() * cbtVector3(0, 1 * side, 0);
            cbtVector3 pos_loc = r1 + cbtVector3(0, H1 * side, 0);
            pos1 = m44T(pos_loc);
            dist = side * (y1 - H1) - radius0;
        }
    }
    /// report a contact. internally this will be kept persistent, and contact reduction is done
    resultOut->addContactPoint(normalOnSurfaceB, pos1, dist);

    resultOut->refreshContactPoints();
}

cbtScalar cbtSphereCylinderCollisionAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                                     cbtCollisionObject* body1,
                                                                     const cbtDispatcherInfo& dispatchInfo,
                                                                     cbtManifoldResult* resultOut) {
    // not yet
    return cbtScalar(1.);
}

void cbtSphereCylinderCollisionAlgorithm::getAllContactManifolds(cbtManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
        manifoldArray.push_back(m_manifoldPtr);
    }
}

cbtCollisionAlgorithm* cbtSphereCylinderCollisionAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* body0Wrap,
    const cbtCollisionObjectWrapper* body1Wrap) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtSphereCylinderCollisionAlgorithm));
    if (!m_swapped) {
        return new (mem) cbtSphereCylinderCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
    } else {
        return new (mem) cbtSphereCylinderCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
    }
}

// ================================================================================================

cbtArcSegmentCollisionAlgorithm::cbtArcSegmentCollisionAlgorithm(cbtPersistentManifold* mf,
                                                                 const cbtCollisionAlgorithmConstructionInfo& ci,
                                                                 const cbtCollisionObjectWrapper* col0,
                                                                 const cbtCollisionObjectWrapper* col1,
                                                                 bool isSwapped)
    : cbtActivatingCollisionAlgorithm(ci, col0, col1), m_ownManifold(false), m_manifoldPtr(mf), m_isSwapped(isSwapped) {
    const cbtCollisionObjectWrapper* arcObjWrap = m_isSwapped ? col1 : col0;
    const cbtCollisionObjectWrapper* segmentObjWrap = m_isSwapped ? col0 : col1;

    if (!m_manifoldPtr &&
        m_dispatcher->needsCollision(arcObjWrap->getCollisionObject(), segmentObjWrap->getCollisionObject())) {
        m_manifoldPtr =
            m_dispatcher->getNewManifold(arcObjWrap->getCollisionObject(), segmentObjWrap->getCollisionObject());
        m_ownManifold = true;
    }
}

cbtArcSegmentCollisionAlgorithm::cbtArcSegmentCollisionAlgorithm(const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtArcSegmentCollisionAlgorithm ::~cbtArcSegmentCollisionAlgorithm() {
    if (m_ownManifold) {
        if (m_manifoldPtr)
            m_dispatcher->releaseManifold(m_manifoldPtr);
    }
}

void cbtArcSegmentCollisionAlgorithm::processCollision(const cbtCollisionObjectWrapper* body0,
                                                       const cbtCollisionObjectWrapper* body1,
                                                       const cbtDispatcherInfo& dispatchInfo,
                                                       cbtManifoldResult* resultOut) {
    (void)dispatchInfo;
    (void)resultOut;
    if (!m_manifoldPtr)
        return;

    const cbtCollisionObjectWrapper* arcObjWrap = m_isSwapped ? body1 : body0;
    const cbtCollisionObjectWrapper* segmentObjWrap = m_isSwapped ? body0 : body1;

    resultOut->setPersistentManifold(m_manifoldPtr);

    // only 1 contact per pair, avoid persistence
    resultOut->getPersistentManifold()->clearManifold();

    const cbt2DarcShape* arc = (cbt2DarcShape*)arcObjWrap->getCollisionShape();
    const cbt2DsegmentShape* segment = (cbt2DsegmentShape*)segmentObjWrap->getCollisionShape();

    // A concave arc (i.e.with outward volume, counterclockwise abscissa) will never collide with segments
    if (arc->get_counterclock())
        return;

    const cbtTransform& m44Tarc = arcObjWrap->getCollisionObject()->getWorldTransform();
    const cbtTransform& m44Tsegment = segmentObjWrap->getCollisionObject()->getWorldTransform();

    // Shapes on two planes that are not so parallel? no collisions!
    cbtVector3 Zarc = m44Tarc.getBasis().getColumn(2);
    cbtVector3 Zsegment = m44Tsegment.getBasis().getColumn(2);
    if (fabs(Zarc.dot(Zsegment)) < 0.99)  //// TODO  threshold as setting
        return;

    // Shapes on two planes that are too far? no collisions!
    cbtVector3 diff = m44Tsegment.invXform(m44Tarc.getOrigin());
    if (fabs(diff.getZ()) > (arc->get_zthickness() + segment->get_zthickness()))
        return;

    // vectors of body 1 in body 2 csys:
    cbtVector3 local_arc_center = m44Tsegment.invXform(m44Tarc * cbtVector3(arc->get_X(), arc->get_Y(), 0));
    cbtVector3 local_arc_X = m44Tsegment.getBasis().transpose() * (m44Tarc.getBasis() * cbtVector3(1, 0, 0));
    double local_arc_rot = atan2(local_arc_X.getY(), local_arc_X.getX());
    double arc1_angle1 = local_arc_rot + arc->get_angle1();
    double arc1_angle2 = local_arc_rot + arc->get_angle2();

    cbtVector3 local_CS1 = local_arc_center - segment->get_P1();
    cbtVector3 local_seg_S2S1 = (segment->get_P2() - segment->get_P1());
    cbtScalar seg_length = local_seg_S2S1.length();
    if (seg_length < 1e-30)
        return;
    cbtVector3 local_seg_D = local_seg_S2S1 / seg_length;
    cbtScalar param = local_CS1.dot(local_seg_D);

    // contact out of segment extrema?
    if (param < 0)
        return;
    if (param > seg_length)
        return;

    cbtVector3 local_P2 = segment->get_P1() + local_seg_D * param;
    cbtVector3 local_CP2 = local_arc_center - local_P2;
    local_CP2.setZ(0);
    cbtVector3 local_R = local_CP2.normalized() * arc->get_radius();
    cbtVector3 local_P1;
    cbtVector3 local_N2;
    if (local_seg_S2S1.cross(local_CP2).getZ() > 0) {
        local_P1 = local_arc_center - local_R;
        local_N2 = local_CP2.normalized();
    } else {
        local_P1 = local_arc_center + local_R;
        local_N2 = -local_CP2.normalized();
    }

    double alpha = atan2(-local_N2.getY(), -local_N2.getX());

    // Discard points out of min-max angles

    // to always positive angles:
    arc1_angle1 = fmod(arc1_angle1 + 1e-30, CH_2PI);
    if (arc1_angle1 < 0)
        arc1_angle1 += CH_2PI;
    arc1_angle2 = fmod(arc1_angle2 + 1e-30, CH_2PI);
    if (arc1_angle2 < 0)
        arc1_angle2 += CH_2PI;
    alpha = fmod(alpha, CH_2PI);
    if (alpha < 0)
        alpha += CH_2PI;

    arc1_angle1 = fmod(arc1_angle1, CH_2PI);
    arc1_angle2 = fmod(arc1_angle2, CH_2PI);

    alpha = fmod(alpha, CH_2PI);

    bool inangle1 = false;

    if (arc1_angle1 < arc1_angle2) {
        if (alpha >= arc1_angle2 || alpha <= arc1_angle1)
            inangle1 = true;
    } else {
        if (alpha >= arc1_angle2 && alpha <= arc1_angle1)
            inangle1 = true;
    }

    if (!inangle1)
        return;

    // transform in absolute coords:
    // cbtVector3 pos1 = m44Tsegment * local_P1; // not needed
    cbtVector3 pos2 = m44Tsegment * local_P2;
    cbtVector3 normal_on_2 = m44Tsegment.getBasis() * local_N2;
    cbtScalar dist = local_N2.dot(local_P1 - local_P2);

    // too far or too interpenetrate? discard.
    if (fabs(dist) > (arc->getMargin() + segment->getMargin()))
        return;

    /// report a contact. internally this will be kept persistent, and contact reduction is done
    resultOut->addContactPoint(normal_on_2, pos2, dist);

    resultOut->refreshContactPoints();
}

cbtScalar cbtArcSegmentCollisionAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                                 cbtCollisionObject* body1,
                                                                 const cbtDispatcherInfo& dispatchInfo,
                                                                 cbtManifoldResult* resultOut) {
    // not yet
    return cbtScalar(1.);
}

void cbtArcSegmentCollisionAlgorithm::getAllContactManifolds(cbtManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
        manifoldArray.push_back(m_manifoldPtr);
    }
}

cbtCollisionAlgorithm* cbtArcSegmentCollisionAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* body0Wrap,
    const cbtCollisionObjectWrapper* body1Wrap) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtArcSegmentCollisionAlgorithm));
    if (!m_swapped) {
        return new (mem) cbtArcSegmentCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
    } else {
        return new (mem) cbtArcSegmentCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
    }
}

// ================================================================================================

cbtArcArcCollisionAlgorithm::cbtArcArcCollisionAlgorithm(cbtPersistentManifold* mf,
                                                         const cbtCollisionAlgorithmConstructionInfo& ci,
                                                         const cbtCollisionObjectWrapper* col0,
                                                         const cbtCollisionObjectWrapper* col1,
                                                         bool isSwapped)
    : cbtActivatingCollisionAlgorithm(ci, col0, col1), m_ownManifold(false), m_manifoldPtr(mf), m_isSwapped(isSwapped) {
    const cbtCollisionObjectWrapper* arcObj1Wrap = m_isSwapped ? col1 : col0;
    const cbtCollisionObjectWrapper* arcObj2Wrap = m_isSwapped ? col0 : col1;

    if (!m_manifoldPtr &&
        m_dispatcher->needsCollision(arcObj1Wrap->getCollisionObject(), arcObj2Wrap->getCollisionObject())) {
        m_manifoldPtr =
            m_dispatcher->getNewManifold(arcObj1Wrap->getCollisionObject(), arcObj2Wrap->getCollisionObject());
        m_ownManifold = true;
    }
}

cbtArcArcCollisionAlgorithm::cbtArcArcCollisionAlgorithm(const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtArcArcCollisionAlgorithm ::~cbtArcArcCollisionAlgorithm() {
    if (m_ownManifold) {
        if (m_manifoldPtr)
            m_dispatcher->releaseManifold(m_manifoldPtr);
    }
}

void cbtArcArcCollisionAlgorithm::processCollision(const cbtCollisionObjectWrapper* body0,
                                                   const cbtCollisionObjectWrapper* body1,
                                                   const cbtDispatcherInfo& dispatchInfo,
                                                   cbtManifoldResult* resultOut) {
    (void)dispatchInfo;
    (void)resultOut;
    if (!m_manifoldPtr)
        return;

    const cbtCollisionObjectWrapper* arcObj1Wrap = m_isSwapped ? body1 : body0;
    const cbtCollisionObjectWrapper* arcObj2Wrap = m_isSwapped ? body0 : body1;

    resultOut->setPersistentManifold(m_manifoldPtr);

    // only 1 contact per pair, avoid persistence
    resultOut->getPersistentManifold()->clearManifold();

    const cbt2DarcShape* arc1 = (cbt2DarcShape*)arcObj1Wrap->getCollisionShape();
    const cbt2DarcShape* arc2 = (cbt2DarcShape*)arcObj2Wrap->getCollisionShape();

    const cbtTransform& m44Tarc1 = arcObj1Wrap->getCollisionObject()->getWorldTransform();
    const cbtTransform& m44Tarc2 = arcObj2Wrap->getCollisionObject()->getWorldTransform();

    // Shapes on two planes that are not so parallel? no collisions!
    cbtVector3 Zarc1 = m44Tarc1.getBasis().getColumn(2);
    cbtVector3 Zarc2 = m44Tarc2.getBasis().getColumn(2);
    if (fabs(Zarc1.dot(Zarc2)) < 0.99)  //// TODO  threshold as setting
        return;

    // Shapes on two planes that are too far? no collisions!
    cbtVector3 diff = m44Tarc2.invXform(m44Tarc1.getOrigin());
    if (fabs(diff.getZ()) > (arc1->get_zthickness() + arc2->get_zthickness()))
        return;

    // vectors and angles of arc 1 in arc 2 csys:
    cbtVector3 local_arc1_center = m44Tarc2.invXform(m44Tarc1 * cbtVector3(arc1->get_X(), arc1->get_Y(), 0));
    cbtVector3 local_arc1_X = m44Tarc2.getBasis().transpose() * (m44Tarc1.getBasis() * cbtVector3(1, 0, 0));
    double local_arc1_rot = atan2(local_arc1_X.getY(), local_arc1_X.getX());
    double arc1_angle1 = local_arc1_rot + arc1->get_angle1();
    double arc1_angle2 = local_arc1_rot + arc1->get_angle2();

    cbtVector3 local_arc2_center = cbtVector3(arc2->get_X(), arc2->get_Y(), 0);
    double arc2_angle1 = arc2->get_angle1();
    double arc2_angle2 = arc2->get_angle2();

    cbtVector3 local_C1C2 = local_arc1_center - local_arc2_center;
    auto len = local_C1C2.length();
    if (len < 1e-8)
        return;
    cbtVector3 local_D12 = local_C1C2 / len;

    cbtVector3 local_P1;
    cbtVector3 local_P2;
    cbtVector3 local_N2;
    double dist = 0;
    bool paired = false;
    double alpha = atan2(local_C1C2.getY(), local_C1C2.getX());
    double alpha1 = 0, alpha2 = 0;

    // convex-convex
    if (arc1->get_counterclock() == false && arc2->get_counterclock() == false) {
        local_P1 = local_arc1_center - local_D12 * arc1->get_radius();
        local_P2 = local_arc2_center + local_D12 * arc2->get_radius();
        local_N2 = local_D12;
        dist = local_C1C2.length() - arc1->get_radius() - arc2->get_radius();
        alpha1 = alpha + CH_PI;
        alpha2 = alpha;
        paired = true;
    }
    // convex-concave
    if (arc1->get_counterclock() == false && arc2->get_counterclock() == true)
        if (arc1->get_radius() <= arc2->get_radius()) {
            local_P1 = local_arc1_center + local_D12 * arc1->get_radius();
            local_P2 = local_arc2_center + local_D12 * arc2->get_radius();
            local_N2 = -local_D12;
            dist = -local_C1C2.length() - arc1->get_radius() + arc2->get_radius();
            alpha1 = alpha;
            alpha2 = alpha;
            paired = true;
        }
    // concave-convex
    if (arc1->get_counterclock() == true && arc2->get_counterclock() == false)
        if (arc1->get_radius() >= arc2->get_radius()) {
            local_P1 = local_arc1_center - local_D12 * arc1->get_radius();
            local_P2 = local_arc2_center - local_D12 * arc2->get_radius();
            local_N2 = -local_D12;
            dist = -local_C1C2.length() + arc1->get_radius() - arc2->get_radius();
            alpha1 = alpha + CH_PI;
            alpha2 = alpha + CH_PI;
            paired = true;
        }

    if (!paired)
        return;

    // Discard points out of min-max angles

    // to always positive angles:
    arc1_angle1 = fmod(arc1_angle1, CH_2PI);
    if (arc1_angle1 < 0)
        arc1_angle1 += CH_2PI;
    arc1_angle2 = fmod(arc1_angle2, CH_2PI);
    if (arc1_angle2 < 0)
        arc1_angle2 += CH_2PI;
    arc2_angle1 = fmod(arc2_angle1, CH_2PI);
    if (arc2_angle1 < 0)
        arc2_angle1 += CH_2PI;
    arc2_angle2 = fmod(arc2_angle2, CH_2PI);
    if (arc2_angle2 < 0)
        arc2_angle2 += CH_2PI;
    alpha1 = fmod(alpha1, CH_2PI);
    if (alpha1 < 0)
        alpha1 += CH_2PI;
    alpha2 = fmod(alpha2, CH_2PI);
    if (alpha2 < 0)
        alpha2 += CH_2PI;

    arc1_angle1 = fmod(arc1_angle1, CH_2PI);
    arc1_angle2 = fmod(arc1_angle2, CH_2PI);
    arc2_angle1 = fmod(arc2_angle1, CH_2PI);
    arc2_angle2 = fmod(arc2_angle2, CH_2PI);
    alpha1 = fmod(alpha1, CH_2PI);
    alpha2 = fmod(alpha2, CH_2PI);

    bool inangle1 = false;
    bool inangle2 = false;

    if (arc1->get_counterclock() == true) {
        if (arc1_angle1 < arc1_angle2) {
            if (alpha1 >= arc1_angle1 && alpha1 <= arc1_angle2)
                inangle1 = true;
        } else {
            if (alpha1 >= arc1_angle1 || alpha1 <= arc1_angle2)
                inangle1 = true;
        }
    } else {
        if (arc1_angle1 < arc1_angle2) {
            if (alpha1 >= arc1_angle2 || alpha1 <= arc1_angle1)
                inangle1 = true;
        } else {
            if (alpha1 >= arc1_angle2 && alpha1 <= arc1_angle1)
                inangle1 = true;
        }
    }

    if (arc2->get_counterclock() == true) {
        if (arc2_angle1 < arc2_angle2) {
            if (alpha2 >= arc2_angle1 && alpha2 <= arc2_angle2)
                inangle2 = true;
        } else {
            if (alpha2 >= arc2_angle1 || alpha2 <= arc2_angle2)
                inangle2 = true;
        }
    } else {
        if (arc2_angle1 < arc2_angle2) {
            if (alpha2 >= arc2_angle2 || alpha2 <= arc2_angle1)
                inangle2 = true;
        } else {
            if (alpha2 >= arc2_angle2 && alpha2 <= arc2_angle1)
                inangle2 = true;
        }
    }

    if (!(inangle1 && inangle2))
        return;

    // transform in absolute coords:
    cbtVector3 pos2 = m44Tarc2 * local_P2;
    cbtVector3 normal_on_2 = m44Tarc2.getBasis() * local_N2;

    // too far or too interpenetrate? discard.
    if (fabs(dist) > (arc1->getMargin() + arc2->getMargin()))
        return;

    /// report a contact.
    resultOut->addContactPoint(normal_on_2, pos2, (cbtScalar)dist);

    resultOut->refreshContactPoints();
}

cbtScalar cbtArcArcCollisionAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                             cbtCollisionObject* body1,
                                                             const cbtDispatcherInfo& dispatchInfo,
                                                             cbtManifoldResult* resultOut) {
    // not yet
    return cbtScalar(1.);
}

void cbtArcArcCollisionAlgorithm::getAllContactManifolds(cbtManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
        manifoldArray.push_back(m_manifoldPtr);
    }
}

cbtCollisionAlgorithm* cbtArcArcCollisionAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* body0Wrap,
    const cbtCollisionObjectWrapper* body1Wrap) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtArcArcCollisionAlgorithm));
    if (!m_swapped) {
        return new (mem) cbtArcArcCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
    } else {
        return new (mem) cbtArcArcCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
    }
}

// ================================================================================================

cbtCEtriangleShapeCollisionAlgorithm::cbtCEtriangleShapeCollisionAlgorithm(
    cbtPersistentManifold* mf,
    const cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* col0,
    const cbtCollisionObjectWrapper* col1,
    bool isSwapped)
    : cbtActivatingCollisionAlgorithm(ci, col0, col1), m_ownManifold(false), m_manifoldPtr(mf), m_isSwapped(isSwapped) {
    const cbtCollisionObjectWrapper* triObj1Wrap = m_isSwapped ? col1 : col0;
    const cbtCollisionObjectWrapper* triObj2Wrap = m_isSwapped ? col0 : col1;

    if (!m_manifoldPtr &&
        m_dispatcher->needsCollision(triObj1Wrap->getCollisionObject(), triObj2Wrap->getCollisionObject())) {
        m_manifoldPtr =
            m_dispatcher->getNewManifold(triObj1Wrap->getCollisionObject(), triObj2Wrap->getCollisionObject());
        m_ownManifold = true;
    }
}

cbtCEtriangleShapeCollisionAlgorithm::cbtCEtriangleShapeCollisionAlgorithm(
    const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtCEtriangleShapeCollisionAlgorithm ::~cbtCEtriangleShapeCollisionAlgorithm() {
    if (m_ownManifold) {
        if (m_manifoldPtr)
            m_dispatcher->releaseManifold(m_manifoldPtr);
    }
}

void cbtCEtriangleShapeCollisionAlgorithm::processCollision(const cbtCollisionObjectWrapper* body0,
                                                            const cbtCollisionObjectWrapper* body1,
                                                            const cbtDispatcherInfo& dispatchInfo,
                                                            cbtManifoldResult* resultOut) {
    (void)dispatchInfo;
    (void)resultOut;
    if (!m_manifoldPtr)
        return;

    const cbtCollisionObjectWrapper* triObj1Wrap = m_isSwapped ? body1 : body0;
    const cbtCollisionObjectWrapper* triObj2Wrap = m_isSwapped ? body0 : body1;

    resultOut->setPersistentManifold(m_manifoldPtr);

    // Avoid persistence of contacts in manifold
    resultOut->getPersistentManifold()->clearManifold();

    const cbtCEtriangleShape* triA = (cbtCEtriangleShape*)triObj1Wrap->getCollisionShape();
    const cbtCEtriangleShape* triB = (cbtCEtriangleShape*)triObj2Wrap->getCollisionShape();
    ChCollisionModelBullet* triModelA = (ChCollisionModelBullet*)triA->getUserPointer();
    ChCollisionModelBullet* triModelB = (ChCollisionModelBullet*)triB->getUserPointer();

    // Discard collisions between connected triangles
    //// TODO: use collision families to bypass during broadphase?
    if (triA->get_p1() == triB->get_p1() || triA->get_p1() == triB->get_p2() || triA->get_p1() == triB->get_p3())
        return;
    if (triA->get_p2() == triB->get_p1() || triA->get_p2() == triB->get_p2() || triA->get_p2() == triB->get_p3())
        return;
    if (triA->get_p3() == triB->get_p1() || triA->get_p3() == triB->get_p2() || triA->get_p3() == triB->get_p3())
        return;

    // Interval boundaries for the distances between vertex-face or edge-edge,
    // these intervals are used to reject distances, where the distance here is assumed for the naked triangles, i.e.
    // WITHOUT the shpereswept_r inflating!
    double max_allowed_dist =
        triModelA->GetEnvelope() + triModelB->GetEnvelope() + triA->sphereswept_r() + triB->sphereswept_r();
    double min_allowed_dist =
        triA->sphereswept_r() + triB->sphereswept_r() - (triModelA->GetSafeMargin() + triModelB->GetSafeMargin());
    double max_edge_dist_earlyout = std::max(max_allowed_dist, std::fabs(min_allowed_dist));

    // Offsets for knowing where the contact points are respect to the points on the naked triangles
    //  - add the sphereswept_r values because one might want to work on the "inflated" triangles for robustness
    //  - TRICK!! offset also by outward 'envelope' because during ReportContacts()
    //    contact points are offset inward by envelope, to cope with GJK method.
    double offset_A = triA->sphereswept_r() + triModelA->GetEnvelope();
    double offset_B = triB->sphereswept_r() + triModelB->GetEnvelope();

    const cbtTransform& m44Ta = triObj1Wrap->getCollisionObject()->getWorldTransform();
    const cbtTransform& m44Tb = triObj2Wrap->getCollisionObject()->getWorldTransform();
    const cbtMatrix3x3& mcbtRa = m44Ta.getBasis();
    const cbtMatrix3x3& mcbtRb = m44Tb.getBasis();
    ChMatrix33<> mRa;
    mRa(0, 0) = mcbtRa[0][0];
    mRa(0, 1) = mcbtRa[0][1];
    mRa(0, 2) = mcbtRa[0][2];
    mRa(1, 0) = mcbtRa[1][0];
    mRa(1, 1) = mcbtRa[1][1];
    mRa(1, 2) = mcbtRa[1][2];
    mRa(2, 0) = mcbtRa[2][0];
    mRa(2, 1) = mcbtRa[2][1];
    mRa(2, 2) = mcbtRa[2][2];
    ChVector3d mOa(m44Ta.getOrigin().x(), m44Ta.getOrigin().y(), m44Ta.getOrigin().z());

    ChMatrix33<> mRb;
    mRb(0, 0) = mcbtRb[0][0];
    mRb(0, 1) = mcbtRb[0][1];
    mRb(0, 2) = mcbtRb[0][2];
    mRb(1, 0) = mcbtRb[1][0];
    mRb(1, 1) = mcbtRb[1][1];
    mRb(1, 2) = mcbtRb[1][2];
    mRb(2, 0) = mcbtRb[2][0];
    mRb(2, 1) = mcbtRb[2][1];
    mRb(2, 2) = mcbtRb[2][2];
    ChVector3d mOb(m44Tb.getOrigin().x(), m44Tb.getOrigin().y(), m44Tb.getOrigin().z());

    // transform points to absolute coords, since models might be roto-translated
    ChVector3d pA1 = mOa + mRa * (*triA->get_p1());
    ChVector3d pA2 = mOa + mRa * (*triA->get_p2());
    ChVector3d pA3 = mOa + mRa * (*triA->get_p3());
    ChVector3d pB1 = mOb + mRb * (*triB->get_p1());
    ChVector3d pB2 = mOb + mRb * (*triB->get_p2());
    ChVector3d pB3 = mOb + mRb * (*triB->get_p3());

    // edges
    ChVector3d eA1 = pA2 - pA1;
    ChVector3d eA2 = pA3 - pA2;
    ChVector3d eA3 = pA1 - pA3;
    ChVector3d eB1 = pB2 - pB1;
    ChVector3d eB2 = pB3 - pB2;
    ChVector3d eB3 = pB1 - pB3;

    // normals
    ChVector3d nA = Vcross(eA1, eA2).GetNormalized();
    ChVector3d nB = Vcross(eB1, eB2).GetNormalized();

    double dist = 1e20;
    bool is_into;
    ChVector3d p_projected;
    double mu, mv;

    // Shortcut: if two degenerate 'skinny' triangles with points 2&3 coincident (ex. used to
    // represent chunks of beams) just do an edge-edge test (as capsule-capsule) and return:
    if ((pA2 == pA3) && (pB2 == pB3) && triA->owns_e1() && triB->owns_e1()) {
        ChVector3d cA, cB, D;
        if (utils::LineLineIntersect(pA1, pA2, pB1, pB2, &cA, &cB, &mu, &mv)) {
            D = cB - cA;
            dist = D.Length();
            if (dist < max_allowed_dist && dist > min_allowed_dist && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                resultOut->refreshContactPoints();
                return;
            }
        }
    }

    // vertex-face tests:

    if (triA->owns_v1()) {
        dist = utils::PointTriangleDistance(pA1, pB1, pB2, pB3, mu, mv, is_into, p_projected);
        if (is_into) {
            if (dist < max_allowed_dist && dist > min_allowed_dist) {
                _add_contact(pA1, p_projected, dist, resultOut, offset_A, offset_B);
            }
        }
    }
    if (triA->owns_v2()) {
        dist = utils::PointTriangleDistance(pA2, pB1, pB2, pB3, mu, mv, is_into, p_projected);
        if (is_into) {
            if (dist < max_allowed_dist && dist > min_allowed_dist) {
                _add_contact(pA2, p_projected, dist, resultOut, offset_A, offset_B);
            }
        }
    }
    if (triA->owns_v3()) {
        dist = utils::PointTriangleDistance(pA3, pB1, pB2, pB3, mu, mv, is_into, p_projected);
        if (is_into) {
            if (dist < max_allowed_dist && dist > min_allowed_dist) {
                _add_contact(pA3, p_projected, dist, resultOut, offset_A, offset_B);
            }
        }
    }

    if (triB->owns_v1()) {
        dist = utils::PointTriangleDistance(pB1, pA1, pA2, pA3, mu, mv, is_into, p_projected);
        if (is_into) {
            if (dist < max_allowed_dist && dist > min_allowed_dist) {
                _add_contact(p_projected, pB1, dist, resultOut, offset_A, offset_B);
            }
        }
    }
    if (triB->owns_v2()) {
        dist = utils::PointTriangleDistance(pB2, pA1, pA2, pA3, mu, mv, is_into, p_projected);
        if (is_into) {
            if (dist < max_allowed_dist && dist > min_allowed_dist) {
                _add_contact(p_projected, pB2, dist, resultOut, offset_A, offset_B);
            }
        }
    }
    if (triB->owns_v3()) {
        dist = utils::PointTriangleDistance(pB3, pA1, pA2, pA3, mu, mv, is_into, p_projected);
        if (is_into) {
            if (dist < max_allowed_dist && dist > min_allowed_dist) {
                _add_contact(p_projected, pB3, dist, resultOut, offset_A, offset_B);
            }
        }
    }
    double beta_A1 = 0, beta_A2 = 0, beta_A3 = 0, beta_B1 = 0, beta_B2 = 0, beta_B3 = 0;  // defaults for free edge
    ChVector3d tA1, tA2, tA3, tB1, tB2, tB3;
    ChVector3d lA1, lA2, lA3, lB1, lB2, lB3;

    //  edges-edges tests

    if (triA->owns_e1()) {
        tA1 = Vcross(eA1, nA).GetNormalized();
        if (triA->get_e1())
            lA1 = (mOa + mRa * (*triA->get_e1())) - pA1;
        else
            lA1 = -tA1;
        beta_A1 = atan2(Vdot(lA1, tA1), Vdot(lA1, nA));
        if (beta_A1 < 0)
            beta_A1 += CH_2PI;
    }
    if (triA->owns_e2()) {
        tA2 = Vcross(eA2, nA).GetNormalized();
        if (triA->get_e2())
            lA2 = (mOa + mRa * (*triA->get_e2())) - pA2;
        else
            lA2 = -tA2;
        beta_A2 = atan2(Vdot(lA2, tA2), Vdot(lA2, nA));
        if (beta_A2 < 0)
            beta_A2 += CH_2PI;
    }
    if (triA->owns_e3()) {
        tA3 = Vcross(eA3, nA).GetNormalized();
        if (triA->get_e3())
            lA3 = (mOa + mRa * (*triA->get_e3())) - pA3;
        else
            lA3 = -tA3;
        beta_A3 = atan2(Vdot(lA3, tA3), Vdot(lA3, nA));
        if (beta_A3 < 0)
            beta_A3 += CH_2PI;
    }
    if (triB->owns_e1()) {
        tB1 = Vcross(eB1, nB).GetNormalized();
        if (triB->get_e1())
            lB1 = (mOb + mRb * (*triB->get_e1())) - pB1;
        else
            lB1 = -tB1;
        beta_B1 = atan2(Vdot(lB1, tB1), Vdot(lB1, nB));
        if (beta_B1 < 0)
            beta_B1 += CH_2PI;
    }
    if (triB->owns_e2()) {
        tB2 = Vcross(eB2, nB).GetNormalized();
        if (triB->get_e2())
            lB2 = (mOb + mRb * (*triB->get_e2())) - pB2;
        else
            lB2 = -tB2;
        beta_B2 = atan2(Vdot(lB2, tB2), Vdot(lB2, nB));
        if (beta_B2 < 0)
            beta_B2 += CH_2PI;
    }
    if (triB->owns_e3()) {
        tB3 = Vcross(eB3, nB).GetNormalized();
        if (triB->get_e3())
            lB3 = (mOb + mRb * (*triB->get_e3())) - pB3;
        else
            lB3 = -tB3;
        beta_B3 = atan2(Vdot(lB3, tB3), Vdot(lB3, nB));
        if (beta_B3 < 0)
            beta_B3 += CH_2PI;
    }

    ChVector3d cA, cB, D;

    double edge_tol = 1e-3;
    //  + edge_tol to discard flat edges with some tolerance:
    double beta_convex_limit = CH_PI_2 + edge_tol;
    //  +/- edge_tol to inflate arc of acceptance of edge vs edge, to cope with singular cases (ex. flat cube vs
    //  flat cube):
    double alpha_lo_limit = -edge_tol;
    double CH_C_PI_mtol = CH_PI - edge_tol;
    double CH_C_PI_2_ptol = CH_PI_2 + edge_tol;

    // edge A1 vs edge B1
    if (triA->owns_e1() && triB->owns_e1())
        if (beta_A1 > beta_convex_limit && beta_B1 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA1, pA2, pB1, pB2, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    double alpha_A = atan2(Vdot(D, tA1), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB1), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A1 - CH_C_PI_2_ptol) && (alpha_B < beta_B1 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A1 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B1 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A1 vs edge B2
    if (triA->owns_e1() && triB->owns_e2())
        if (beta_A1 > beta_convex_limit && beta_B2 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA1, pA2, pB2, pB3, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA1), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB2), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A1 - CH_C_PI_2_ptol) && (alpha_B < beta_B2 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A1 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B2 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A1 vs edge B3
    if (triA->owns_e1() && triB->owns_e3())
        if (beta_A1 > beta_convex_limit && beta_B3 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA1, pA2, pB3, pB1, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA1), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB3), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A1 - CH_C_PI_2_ptol) && (alpha_B < beta_B3 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A1 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B3 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A2 vs edge B1
    if (triA->owns_e2() && triB->owns_e1())
        if (beta_A2 > beta_convex_limit && beta_B1 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA2, pA3, pB1, pB2, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA2), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB1), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A2 - CH_C_PI_2_ptol) && (alpha_B < beta_B1 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A2 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B1 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A2 vs edge B2
    if (triA->owns_e2() && triB->owns_e2())
        if (beta_A2 > beta_convex_limit && beta_B2 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA2, pA3, pB2, pB3, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA2), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB2), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A2 - CH_C_PI_2_ptol) && (alpha_B < beta_B2 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A2 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B2 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A2 vs edge B3
    if (triA->owns_e2() && triB->owns_e3())
        if (beta_A2 > beta_convex_limit && beta_B3 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA2, pA3, pB3, pB1, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA2), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB3), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A2 - CH_C_PI_2_ptol) && (alpha_B < beta_B3 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A2 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B3 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A3 vs edge B1
    if (triA->owns_e3() && triB->owns_e1())
        if (beta_A3 > beta_convex_limit && beta_B1 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA3, pA1, pB1, pB2, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA3), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB1), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A3 - CH_C_PI_2_ptol) && (alpha_B < beta_B1 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A3 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B1 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A3 vs edge B2
    if (triA->owns_e3() && triB->owns_e2())
        if (beta_A3 > beta_convex_limit && beta_B2 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA3, pA1, pB2, pB3, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA3), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB2), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A3 - CH_C_PI_2_ptol) && (alpha_B < beta_B2 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A3 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B2 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }
    // edge A3 vs edge B3
    if (triA->owns_e3() && triB->owns_e3())
        if (beta_A3 > beta_convex_limit && beta_B3 > beta_convex_limit) {
            if (utils::LineLineIntersect(pA3, pA1, pB3, pB1, &cA, &cB, &mu, &mv)) {
                D = cB - cA;
                dist = D.Length();
                if (dist < max_edge_dist_earlyout && mu > 0 && mu < 1 && mv > 0 && mv < 1) {
                    D = cB - cA;
                    double alpha_A = atan2(Vdot(D, tA3), Vdot(D, nA));
                    double alpha_B = atan2(Vdot(-D, tB3), Vdot(-D, nB));
                    if (alpha_A < alpha_lo_limit)
                        alpha_A += CH_2PI;
                    if (alpha_B < alpha_lo_limit)
                        alpha_B += CH_2PI;
                    if ((alpha_A < beta_A3 - CH_C_PI_2_ptol) && (alpha_B < beta_B3 - CH_C_PI_2_ptol)) {
                        if (dist < max_allowed_dist && dist > min_allowed_dist)  // distance interval check - outside
                            _add_contact(cA, cB, dist, resultOut, offset_A, offset_B);
                    } else if (alpha_A > CH_C_PI_mtol && (alpha_A < beta_A3 + CH_PI_2) && alpha_B > CH_C_PI_mtol &&
                               (alpha_B < beta_B3 + CH_C_PI_2_ptol)) {
                        if (-dist < max_allowed_dist && -dist > min_allowed_dist)  // distance interval check - inside
                            _add_contact(cA, cB, -dist, resultOut, offset_A, offset_B);
                    }
                }
            }
        }

    resultOut->refreshContactPoints();
}

cbtScalar cbtCEtriangleShapeCollisionAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                                      cbtCollisionObject* body1,
                                                                      const cbtDispatcherInfo& dispatchInfo,
                                                                      cbtManifoldResult* resultOut) {
    // not yet
    return cbtScalar(1.);
}

void cbtCEtriangleShapeCollisionAlgorithm::getAllContactManifolds(cbtManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
        manifoldArray.push_back(m_manifoldPtr);
    }
}

void cbtCEtriangleShapeCollisionAlgorithm::_add_contact(const ChVector3d& candid_pA,
                                                        const ChVector3d& candid_pB,
                                                        const double dist,
                                                        cbtManifoldResult* resultOut,
                                                        const double offsetA,
                                                        const double offsetB) {
    // convert to Bullet vectors. Note: in absolute csys.
    cbtVector3 absA((cbtScalar)candid_pA.x(), (cbtScalar)candid_pA.y(), (cbtScalar)candid_pA.z());
    cbtVector3 absB((cbtScalar)candid_pB.x(), (cbtScalar)candid_pB.y(), (cbtScalar)candid_pB.z());
    ChVector3d dabsN_onB((candid_pA - candid_pB).GetNormalized());
    cbtVector3 absN_onB((cbtScalar)dabsN_onB.x(), (cbtScalar)dabsN_onB.y(), (cbtScalar)dabsN_onB.z());
    if (dist < 0)
        absN_onB = -absN_onB;  // flip norm to be coherent with dist sign
    resultOut->addContactPoint(absN_onB, absB + absN_onB * (cbtScalar)offsetB, (cbtScalar)(dist - (offsetA + offsetB)));
}

cbtCollisionAlgorithm* cbtCEtriangleShapeCollisionAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* body0Wrap,
    const cbtCollisionObjectWrapper* body1Wrap) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtCEtriangleShapeCollisionAlgorithm));
    if (!m_swapped) {
        return new (mem) cbtCEtriangleShapeCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
    } else {
        return new (mem) cbtCEtriangleShapeCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
    }
}

// ================================================================================================

cbtSegmentSegmentCollisionAlgorithm::cbtSegmentSegmentCollisionAlgorithm(
    cbtPersistentManifold* mf,
    const cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* col0,
    const cbtCollisionObjectWrapper* col1,
    bool isSwapped)
    : cbtActivatingCollisionAlgorithm(ci, col0, col1), m_ownManifold(false), m_manifoldPtr(mf), m_isSwapped(isSwapped) {
    const cbtCollisionObjectWrapper* triObj1Wrap = m_isSwapped ? col1 : col0;
    const cbtCollisionObjectWrapper* triObj2Wrap = m_isSwapped ? col0 : col1;

    if (!m_manifoldPtr &&
        m_dispatcher->needsCollision(triObj1Wrap->getCollisionObject(), triObj2Wrap->getCollisionObject())) {
        m_manifoldPtr =
            m_dispatcher->getNewManifold(triObj1Wrap->getCollisionObject(), triObj2Wrap->getCollisionObject());
        m_ownManifold = true;
    }
}

cbtSegmentSegmentCollisionAlgorithm::cbtSegmentSegmentCollisionAlgorithm(
    const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtSegmentSegmentCollisionAlgorithm ::~cbtSegmentSegmentCollisionAlgorithm() {
    if (m_ownManifold) {
        if (m_manifoldPtr)
            m_dispatcher->releaseManifold(m_manifoldPtr);
    }
}

void cbtSegmentSegmentCollisionAlgorithm::processCollision(const cbtCollisionObjectWrapper* body0,
                                                           const cbtCollisionObjectWrapper* body1,
                                                           const cbtDispatcherInfo& dispatchInfo,
                                                           cbtManifoldResult* resultOut) {
    (void)dispatchInfo;
    (void)resultOut;
    if (!m_manifoldPtr)
        return;

    const cbtCollisionObjectWrapper* triObj1Wrap = m_isSwapped ? body1 : body0;
    const cbtCollisionObjectWrapper* triObj2Wrap = m_isSwapped ? body0 : body1;

    resultOut->setPersistentManifold(m_manifoldPtr);

    // Avoid persistence of contacts in manifold
    resultOut->getPersistentManifold()->clearManifold();

    const cbtSegmentShape* segA = (cbtSegmentShape*)triObj1Wrap->getCollisionShape();
    const cbtSegmentShape* segB = (cbtSegmentShape*)triObj2Wrap->getCollisionShape();
    ChCollisionModelBullet* segModelA = (ChCollisionModelBullet*)segA->getUserPointer();
    ChCollisionModelBullet* segModelB = (ChCollisionModelBullet*)segB->getUserPointer();

    // Discard collisions between connected segments
    //// TODO: use collision families to bypass during broadphase?
    if (segA->get_p1() == segB->get_p1() || segA->get_p1() == segB->get_p2())
        return;
    if (segA->get_p2() == segB->get_p1() || segA->get_p2() == segB->get_p2())
        return;

    //// TODO
    return;
}

cbtScalar cbtSegmentSegmentCollisionAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                                     cbtCollisionObject* body1,
                                                                     const cbtDispatcherInfo& dispatchInfo,
                                                                     cbtManifoldResult* resultOut) {
    // not yet
    return cbtScalar(1.);
}

void cbtSegmentSegmentCollisionAlgorithm::getAllContactManifolds(cbtManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
        manifoldArray.push_back(m_manifoldPtr);
    }
}
void cbtSegmentSegmentCollisionAlgorithm::_add_contact(const ChVector3d& candid_pA,
                                                       const ChVector3d& candid_pB,
                                                       const double dist,
                                                       cbtManifoldResult* resultOut,
                                                       const double offsetA,
                                                       const double offsetB) {
    //// TODO
}

cbtCollisionAlgorithm* cbtSegmentSegmentCollisionAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* body0Wrap,
    const cbtCollisionObjectWrapper* body1Wrap) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtSegmentSegmentCollisionAlgorithm));
    if (!m_swapped) {
        return new (mem) cbtSegmentSegmentCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, false);
    } else {
        return new (mem) cbtSegmentSegmentCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, true);
    }
}




// ---------------------------------------------------------------------------
// Heightfield collision algorithm
cbtCollisionAlgorithm* cbtConvexHeightfieldAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* a,
    const cbtCollisionObjectWrapper* b) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtConvexHeightfieldAlgorithm));
    return new (mem) cbtConvexHeightfieldAlgorithm(ci.m_manifold, ci, a, b, false);
}
cbtCollisionAlgorithm* cbtConvexHeightfieldAlgorithm::SwappedCreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* a,
    const cbtCollisionObjectWrapper* b) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtConvexHeightfieldAlgorithm));
    return new (mem) cbtConvexHeightfieldAlgorithm(ci.m_manifold, ci, a, b, true);
}

cbtConvexHeightfieldAlgorithm::cbtConvexHeightfieldAlgorithm(cbtPersistentManifold* mf,
                                                   const cbtCollisionAlgorithmConstructionInfo& ci,
                                                   const cbtCollisionObjectWrapper* a,
                                                   const cbtCollisionObjectWrapper* b,
                                                   bool swapped)
    : cbtActivatingCollisionAlgorithm(ci, a, b) {
    m_convex = static_cast<const cbtConvexShape*>(swapped ? b->getCollisionShape() : a->getCollisionShape());
    m_ownManifold = (mf == nullptr);
    m_manifoldPtr =
        m_ownManifold ? ci.m_dispatcher1->getNewManifold(a->getCollisionObject(), b->getCollisionObject()) : mf;
}

cbtConvexHeightfieldAlgorithm::cbtConvexHeightfieldAlgorithm(const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtConvexHeightfieldAlgorithm::~cbtConvexHeightfieldAlgorithm() {
    if (m_ownManifold && m_manifoldPtr)
        m_dispatcher->releaseManifold(m_manifoldPtr);
}

void cbtConvexHeightfieldAlgorithm::collectSupportPoints(const cbtConvexShape* s,
                                                              cbtAlignedObjectArray<cbtVector3>& out) {
    if (s->getShapeType() == BOX_SHAPE_PROXYTYPE) {
        // Use support function sampling for boxes to ensure correct geometry
        static const cbtVector3 directions[8] = {cbtVector3(1, 1, 1),   cbtVector3(-1, 1, 1),  cbtVector3(1, -1, 1),
                                                 cbtVector3(-1, -1, 1), cbtVector3(1, 1, -1),  cbtVector3(-1, 1, -1),
                                                 cbtVector3(1, -1, -1), cbtVector3(-1, -1, -1)};
        out.resize(8);
        for (int i = 0; i < 8; ++i) {
            out[i] = s->localGetSupportingVertexWithoutMargin(directions[i]);
        }
        return;
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
    {
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
}



// ============================================================================
// ADAPTIVE CONVEX-HEIGHTFIELD WITH TERRAIN VARIATION CHECK
// ============================================================================

void cbtConvexHeightfieldAlgorithm::processCollision(const cbtCollisionObjectWrapper* wrapperA,
                                                     const cbtCollisionObjectWrapper* wrapperB,
                                                     const cbtDispatcherInfo& /*info*/,
                                                     cbtManifoldResult* resultOut) {
    if (!m_manifoldPtr)
        return;
    resultOut->setPersistentManifold(m_manifoldPtr);

    bool convexIsA = (wrapperA->getCollisionShape() == m_convex);
    const auto* convexWrap = convexIsA ? wrapperA : wrapperB;
    const auto* terrainWrap = convexIsA ? wrapperB : wrapperA;
    const cbtTransform& convexTransform = convexWrap->getWorldTransform();
    const cbtTransform& terrainTransform = terrainWrap->getWorldTransform();
    auto* terrainShape = static_cast<const cbtHeightfieldChronoTerrainShape*>(terrainWrap->getCollisionShape());
    if (!terrainShape)
        return;

    const cbtScalar keepSlop = m_manifoldPtr->getContactBreakingThreshold() * 0.5f;
    auto* convexShape = static_cast<const cbtConvexShape*>(convexWrap->getCollisionShape());
    const cbtScalar convexMargin = convexShape->getMargin();
    const cbtScalar terrainMargin = terrainShape->getMargin();

    const int upAxis = terrainShape->getUpAxis();
    cbtVector3 worldUp(0, 0, 0);
    worldUp[upAxis] = 1.0f;
    cbtVector3 terrainWorldUp = terrainTransform.getBasis() * worldUp;
    terrainWorldUp.normalize();

    // ============================================================================
    // Compute AABB for terrain variation check
    // ============================================================================

    cbtVector3 aabbMin, aabbMax;
    convexShape->getAabb(convexTransform, aabbMin, aabbMax);
    cbtScalar expand = convexMargin + terrainMargin + keepSlop;
    aabbMin -= cbtVector3(expand, expand, expand);
    aabbMax += cbtVector3(expand, expand, expand);

    // Estimate characteristic size of convex shape
    cbtVector3 aabbExtent = aabbMax - aabbMin;
    cbtScalar characteristicSize = (aabbExtent.x() + aabbExtent.y() + aabbExtent.z()) / 3.0f;

    // ============================================================================
    // STEP 2: Check terrain variation in projected AABB
    // Sample 9 points (center + 8 around perimeter) to estimate complexity
    // ============================================================================

    cbtVector3 aabbCenter = (aabbMin + aabbMax) * 0.5f;
    cbtVector3 aabbHalfExtent = (aabbMax - aabbMin) * 0.5f;

    struct TerrainVariationSample {
        cbtScalar height;
        cbtVector3 normal;
        bool valid;
    };

    TerrainVariationSample variationSamples[9];
    int validVariationSamples = 0;

    // Sample positions: center + 8 corners of AABB projection
    cbtVector3 sampleOffsets[9];
    sampleOffsets[0] = cbtVector3(0, 0, 0);  // Center
    sampleOffsets[1] = cbtVector3(1, 0, 1);
    sampleOffsets[2] = cbtVector3(-1, 0, 1);
    sampleOffsets[3] = cbtVector3(1, 0, -1);
    sampleOffsets[4] = cbtVector3(-1, 0, -1);
    sampleOffsets[5] = cbtVector3(1, 0, 0);
    sampleOffsets[6] = cbtVector3(-1, 0, 0);
    sampleOffsets[7] = cbtVector3(0, 0, 1);
    sampleOffsets[8] = cbtVector3(0, 0, -1);

    for (int i = 0; i < 9; ++i) {
        cbtVector3 sampleQuery = aabbCenter + cbtVector3(sampleOffsets[i].x() * aabbHalfExtent.x(),
                                                         sampleOffsets[i].y() * aabbHalfExtent.y(),
                                                         sampleOffsets[i].z() * aabbHalfExtent.z());

        cbtVector3 terrainPt, terrainNormal;
        if (terrainShape->sampleWorld(terrainTransform, sampleQuery, terrainPt, terrainNormal)) {
            if (terrainNormal.length2() >= SIMD_EPSILON) {
                terrainNormal.normalize();
                variationSamples[i].height = terrainPt[upAxis];
                variationSamples[i].normal = terrainNormal;
                variationSamples[i].valid = true;
                validVariationSamples++;
            }
        }
    }

    if (validVariationSamples < 3)
        return;  // Not enough terrain data

    // Calculate variation metrics
    cbtScalar minHeight = variationSamples[0].height;
    cbtScalar maxHeight = variationSamples[0].height;
    cbtScalar avgNormalDot = 0.0f;
    int normalComparisons = 0;

    for (int i = 0; i < 9; ++i) {
        if (!variationSamples[i].valid)
            continue;

        minHeight = cbtMin(minHeight, variationSamples[i].height);
        maxHeight = cbtMax(maxHeight, variationSamples[i].height);

        for (int j = i + 1; j < 9; ++j) {
            if (!variationSamples[j].valid)
                continue;
            avgNormalDot += variationSamples[i].normal.dot(variationSamples[j].normal);
            normalComparisons++;
        }
    }

    if (normalComparisons > 0)
        avgNormalDot /= normalComparisons;

    // Terrain complexity score
    cbtScalar heightVariation = (maxHeight - minHeight) / (characteristicSize * cbtScalar(0.01));
    cbtScalar normalVariation = cbtScalar(1.0) - avgNormalDot;
    cbtScalar terrainComplexity = heightVariation + normalVariation * 10.0f;

    // Decide support sampling density based on complexity
    int numSupportDirs;
    int numFaceSamples;

    if (terrainComplexity < 0.5f) {
        numSupportDirs = 5;  // Flat: minimal (center + 4 cardinal)
        numFaceSamples = 2;  // Only primary faces
    } else if (terrainComplexity < 2.0f) {
        numSupportDirs = 9;  // Moderate (center + 8 radial)
        numFaceSamples = 4;
    } else {
        numSupportDirs = 13;  // Complex (full sampling)
        numFaceSamples = 8;
    }

    // ============================================================================
    // STEP 3: Generate contact candidates with adaptive sampling
    // ============================================================================

    struct ContactCandidate {
        cbtVector3 convexPoint;
        cbtVector3 terrainPoint;
        cbtVector3 normal;
        cbtScalar penetration;
        int priority;
    };

    cbtAlignedObjectArray<ContactCandidate> candidates;
    candidates.reserve(32);

    // STAGE 1: COM
    {
        cbtVector3 worldCOM = convexTransform.getOrigin();
        cbtVector3 terrainPt, terrainNormal;

        if (terrainShape->sampleWorld(terrainTransform, worldCOM, terrainPt, terrainNormal)) {
            if (terrainNormal.length2() >= SIMD_EPSILON) {
                terrainNormal.normalize();
                cbtVector3 comToTerrain = terrainPt - worldCOM;
                cbtScalar signedDist = comToTerrain.dot(terrainNormal);
                cbtScalar penetration = -signedDist - convexMargin - terrainMargin;

                if (penetration < keepSlop) {
                    cbtVector3 convexPt = worldCOM + terrainNormal * signedDist;
                    candidates.push_back({convexPt, terrainPt, terrainNormal, penetration, 0});
                }
            }
        }
    }

    // STAGE 2: Support points (adaptive count)
    {
        cbtVector3 supportDirs[13];
        supportDirs[0] = -terrainWorldUp;

        cbtVector3 perp1 = (cbtFabs(terrainWorldUp.x()) < 0.9f)
                               ? cbtVector3(1, 0, 0).cross(terrainWorldUp).normalized()
                               : cbtVector3(0, 1, 0).cross(terrainWorldUp).normalized();
        cbtVector3 perp2 = terrainWorldUp.cross(perp1).normalized();

        if (numSupportDirs >= 5) {
            // 4 cardinal
            supportDirs[1] = perp1;
            supportDirs[2] = -perp1;
            supportDirs[3] = perp2;
            supportDirs[4] = -perp2;
        }

        if (numSupportDirs >= 9) {
            // 4 diagonal
            for (int i = 0; i < 4; ++i) {
                cbtScalar angle = i * SIMD_PI / 2.0f + SIMD_PI / 4.0f;
                cbtScalar tilt = SIMD_PI / 6.0f;
                supportDirs[5 + i] =
                    (-terrainWorldUp * cbtCos(tilt) + (perp1 * cbtCos(angle) + perp2 * cbtSin(angle)) * cbtSin(tilt))
                        .normalized();
            }
        }

        if (numSupportDirs >= 13) {
            // 4 mid-level
            for (int i = 0; i < 4; ++i) {
                cbtScalar angle = i * SIMD_PI / 2.0f;
                cbtScalar tilt = SIMD_PI / 4.0f;
                supportDirs[9 + i] =
                    (-terrainWorldUp * cbtCos(tilt) + (perp1 * cbtCos(angle) + perp2 * cbtSin(angle)) * cbtSin(tilt))
                        .normalized();
            }
        }

        for (int i = 0; i < numSupportDirs; ++i) {
            cbtVector3 localDir = convexTransform.getBasis().transpose() * supportDirs[i];
            cbtVector3 supportPtLocal = convexShape->localGetSupportingVertexWithoutMargin(localDir);
            cbtVector3 supportPtWorld = convexTransform * supportPtLocal;

            cbtVector3 terrainPt, terrainNormal;
            if (terrainShape->sampleWorld(terrainTransform, supportPtWorld, terrainPt, terrainNormal)) {
                if (terrainNormal.length2() >= SIMD_EPSILON) {
                    terrainNormal.normalize();
                    cbtVector3 supportToTerrain = terrainPt - supportPtWorld;
                    cbtScalar signedDist = supportToTerrain.dot(terrainNormal);
                    cbtScalar penetration = -signedDist - convexMargin - terrainMargin;

                    if (penetration < keepSlop) {
                        candidates.push_back({supportPtWorld, terrainPt, terrainNormal, penetration, 1});
                    }
                }
            }
        }
    }

    // STAGE 3: Face sampling (adaptive count)
    auto* polyShape = dynamic_cast<const cbtPolyhedralConvexShape*>(convexShape);
    if (polyShape) {
        cbtVector3 terrainUpLocal = convexTransform.getBasis().transpose() * terrainWorldUp;
        int numFaces = cbtMin(polyShape->getNumPlanes(), numFaceSamples);

        for (int faceIdx = 0; faceIdx < numFaces; ++faceIdx) {
            cbtVector3 faceNormal, facePt;
            polyShape->getPlane(faceNormal, facePt, faceIdx);

            if (faceNormal.dot(terrainUpLocal) > -0.1f)
                continue;

            cbtVector3 worldFacePt = convexTransform * facePt;

            cbtVector3 terrainPt, terrainNormal;
            if (terrainShape->sampleWorld(terrainTransform, worldFacePt, terrainPt, terrainNormal)) {
                if (terrainNormal.length2() >= SIMD_EPSILON) {
                    terrainNormal.normalize();
                    cbtVector3 faceToTerrain = terrainPt - worldFacePt;
                    cbtScalar signedDist = faceToTerrain.dot(terrainNormal);
                    cbtScalar penetration = -signedDist - convexMargin - terrainMargin;

                    if (penetration < keepSlop) {
                        candidates.push_back({worldFacePt, terrainPt, terrainNormal, penetration, 2});
                    }
                }
            }
        }
    }

    if (candidates.size() == 0)
        return;

    // ============================================================================
    // STEP 4: Contact reduction (same as before)
    // ============================================================================

    struct ContactSorter {
        bool operator()(const ContactCandidate& a, const ContactCandidate& b) const {
            if (a.priority != b.priority)
                return a.priority < b.priority;
            return a.penetration < b.penetration;
        }
    };
    candidates.quickSort(ContactSorter());

    const int MAX_CONTACTS = 4;
    const cbtScalar MIN_SEPARATION = convexMargin * 0.5f;
    const cbtScalar NORMAL_COS_THRESHOLD = cbtCos(15.0f * SIMD_PI / 180.0f);

    cbtAlignedObjectArray<ContactCandidate> finalContacts;
    finalContacts.reserve(MAX_CONTACTS);

    for (int i = 0; i < candidates.size() && finalContacts.size() < MAX_CONTACTS; ++i) {
        if (candidates[i].penetration > 0.0f)
            break;

        bool unique = true;
        for (int j = 0; j < finalContacts.size(); ++j) {
            cbtScalar posDist = (candidates[i].terrainPoint - finalContacts[j].terrainPoint).length();
            cbtScalar normalDot = candidates[i].normal.dot(finalContacts[j].normal);

            if (posDist < MIN_SEPARATION && normalDot > NORMAL_COS_THRESHOLD) {
                if (candidates[i].penetration < finalContacts[j].penetration) {
                    finalContacts[j] = candidates[i];
                }
                unique = false;
                break;
            }
        }

        if (unique) {
            finalContacts.push_back(candidates[i]);
        }
    }

    // ============================================================================
    // STEP 5: Emit contacts
    // ============================================================================

    for (int i = 0; i < finalContacts.size(); ++i) {
        const ContactCandidate& c = finalContacts[i];

        cbtVector3 contactNormal = c.normal;
        cbtVector3 contactPoint = c.convexPoint;
        cbtScalar penetration = c.penetration;

        if (!convexIsA) {
            contactNormal = -contactNormal;
        }

        resultOut->addContactPoint(contactNormal, contactPoint, penetration);
    }

    resultOut->refreshContactPoints();
}


// CCD support: Calculate time of impact (TOI) for convex against heightfield
cbtScalar cbtConvexHeightfieldAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                               cbtCollisionObject* body1,
                                                               const cbtDispatcherInfo& /*dispatchInfo*/,
                                                               cbtManifoldResult* /*resultOut*/) {
    bool terrainIsA = (body0->getCollisionShape()->getShapeType() == TERRAIN_SHAPE_PROXYTYPE);
    auto* terrain = terrainIsA ? static_cast<const cbtHeightfieldChronoTerrainShape*>(body0->getCollisionShape())
                               : static_cast<const cbtHeightfieldChronoTerrainShape*>(body1->getCollisionShape());
    auto* convex = terrainIsA ? static_cast<const cbtConvexShape*>(body1->getCollisionShape())
                              : static_cast<const cbtConvexShape*>(body0->getCollisionShape());

    const cbtTransform& fromTransform = terrainIsA ? body1->getWorldTransform() : body0->getWorldTransform();
    const cbtTransform& toTransform =
        terrainIsA ? body1->getInterpolationWorldTransform() : body0->getInterpolationWorldTransform();
    const cbtTransform& terrainTransform = terrainIsA ? body0->getWorldTransform() : body1->getWorldTransform();

    cbtVector3 motion = toTransform.getOrigin() - fromTransform.getOrigin();
    if (motion.length2() < SIMD_EPSILON)
        return 1.f;  // No motion: No impact

    int upAxis = terrain->getUpAxis();
    cbtScalar totalMargin = convex->getMargin() + terrain->getMargin() + 0.01f;  // Extra tolerance for early detection
    cbtVector3 invScale = terrain->getInverseLocalScaling();
    cbtVector3 localOrigin = terrain->getLocalOrigin();

    // Binary search along motion path
    cbtScalar low = 0.f, high = 1.f;
    for (int iter = 0; iter < 25; ++iter) {  // ~1e-7 precision
        cbtScalar mid = 0.5f * (low + high);
        cbtVector3 interpolatedPos = fromTransform.getOrigin() + motion * mid;
        cbtVector3 localPos = (terrainTransform.invXform(interpolatedPos) * invScale) + localOrigin;

        // For general convex: Use support in -terrainUp (deepest point toward terrain)
        cbtVector3 terrainUp(0, 0, 0);
        terrainUp[upAxis] = 1.f;
        cbtVector3 supportDirLocal = fromTransform.getBasis().transpose() * (-terrainUp);  // To convex local
        cbtVector3 supportVertexLocal = convex->localGetSupportingVertexWithoutMargin(supportDirLocal);
        cbtVector3 supportVertexWorld = fromTransform * supportVertexLocal + motion * mid;  // Interpolated

        // Sample terrain at support point (world)
        cbtVector3 terrainSurface, terrainNormal;
        if (!terrain->sampleWorld(terrainTransform, supportVertexWorld, terrainSurface, terrainNormal))
            continue;  // Outside bounds: No hit

        // Signed distance: Positive above, negative penetrating
        cbtScalar distToSurface = (supportVertexWorld - terrainSurface).dot(terrainNormal) - totalMargin;
        bool penetrating = distToSurface <= 0.f;
        (penetrating ? high : low) = mid;
    }
    return high;  // TOI fraction [0,1]; <1 means impact
}
////////////////////////////////////////////////////



// ---------------------------------------------------------------------------
// Heightfield collision algorithm
cbtCollisionAlgorithm* cbtSphereHeightfieldAlgorithm::CreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* a,
    const cbtCollisionObjectWrapper* b) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtSphereHeightfieldAlgorithm));
    return new (mem) cbtSphereHeightfieldAlgorithm(ci.m_manifold, ci, a, b, false);
}
cbtCollisionAlgorithm* cbtSphereHeightfieldAlgorithm::SwappedCreateFunc::CreateCollisionAlgorithm(
    cbtCollisionAlgorithmConstructionInfo& ci,
    const cbtCollisionObjectWrapper* a,
    const cbtCollisionObjectWrapper* b) {
    void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(cbtSphereHeightfieldAlgorithm));
    return new (mem) cbtSphereHeightfieldAlgorithm(ci.m_manifold, ci, a, b, true);
}

cbtSphereHeightfieldAlgorithm::cbtSphereHeightfieldAlgorithm(cbtPersistentManifold* mf,
                                                             const cbtCollisionAlgorithmConstructionInfo& ci,
                                                             const cbtCollisionObjectWrapper* a,
                                                             const cbtCollisionObjectWrapper* b,
                                                             bool swapped)
    : cbtActivatingCollisionAlgorithm(ci, a, b) {
    m_convex = static_cast<const cbtConvexShape*>(swapped ? b->getCollisionShape() : a->getCollisionShape());
    m_ownManifold = (mf == nullptr);
    m_manifoldPtr =
        m_ownManifold ? ci.m_dispatcher1->getNewManifold(a->getCollisionObject(), b->getCollisionObject()) : mf;
}

cbtSphereHeightfieldAlgorithm::cbtSphereHeightfieldAlgorithm(const cbtCollisionAlgorithmConstructionInfo& ci)
    : cbtActivatingCollisionAlgorithm(ci) {}

cbtSphereHeightfieldAlgorithm::~cbtSphereHeightfieldAlgorithm() {
    if (m_ownManifold && m_manifoldPtr)
        m_dispatcher->releaseManifold(m_manifoldPtr);
}



// ============================================================================
// ADAPTIVE SPHERE-HEIGHTFIELD WITH TERRAIN VARIATION CHECK
// ======================================================
// 
// TODO:- this is messy, complicated and is not efficient in usage of data handling and physcis processing.

void cbtSphereHeightfieldAlgorithm::processCollision(const cbtCollisionObjectWrapper* wrapperA,
                                                     const cbtCollisionObjectWrapper* wrapperB,
                                                     const cbtDispatcherInfo& /*info*/,
                                                     cbtManifoldResult* resultOut) {
    if (!m_manifoldPtr)
        return;
    resultOut->setPersistentManifold(m_manifoldPtr);

    const bool sphereIsA = (wrapperA->getCollisionShape() == m_convex);
    const auto* sphereWrap = sphereIsA ? wrapperA : wrapperB;
    const auto* terrainWrap = sphereIsA ? wrapperB : wrapperA;

    const cbtTransform& sphereTransform = sphereWrap->getWorldTransform();
    const cbtTransform& terrainTransform = terrainWrap->getWorldTransform();
    const auto* terrainShape = static_cast<const cbtHeightfieldChronoTerrainShape*>(terrainWrap->getCollisionShape());
    if (!terrainShape)
        return;

    const cbtScalar contactSlop = m_manifoldPtr->getContactBreakingThreshold() * cbtScalar(0.5);
    const auto* sphereShape = static_cast<const cbtSphereShape*>(m_convex);
    const cbtScalar radius = sphereShape->getRadius();
    const cbtScalar terrainMargin = terrainShape->getMargin();
    const cbtVector3 sphereCenter = sphereTransform.getOrigin();

    const int upAxis = terrainShape->getUpAxis();
    cbtVector3 worldUp(0, 0, 0);
    worldUp[upAxis] = 1.0f;
    cbtVector3 terrainWorldUp = terrainTransform.getBasis() * worldUp;
    terrainWorldUp.normalize();

    // ============================================================================
    // STEP 1: Sample terrain at sphere center + check local variation
    // ============================================================================

    cbtVector3 primaryTerrainPt, primaryNormal;
    if (!terrainShape->sampleWorld(terrainTransform, sphereCenter, primaryTerrainPt, primaryNormal))
        return;

    if (primaryNormal.length2() < SIMD_EPSILON)
        return;

    primaryNormal.normalize();

    cbtVector3 centerToTerrain = primaryTerrainPt - sphereCenter;
    cbtScalar centerDist = centerToTerrain.length();

    if (centerDist > radius + terrainMargin + contactSlop)
        return;

    cbtVector3 primaryContactDir;
    if (centerDist > SIMD_EPSILON) {
        primaryContactDir = centerToTerrain / centerDist;
    } else {
        primaryContactDir = -primaryNormal;
    }

    // ============================================================================
    // STEP 2: ADAPTIVE SAMPLING - Check terrain variation in AABB
    // Sample 5 points (center + 4 cardinal) to estimate terrain complexity
    // ============================================================================

    struct TerrainVariationSample {
        cbtVector3 point;
        cbtVector3 normal;
        cbtScalar height;
        bool valid;
    };

    TerrainVariationSample variationSamples[5];
    int validVariationSamples = 0;

    // Center
    variationSamples[0].point = primaryTerrainPt;
    variationSamples[0].normal = primaryNormal;
    variationSamples[0].height = primaryTerrainPt[upAxis];
    variationSamples[0].valid = true;
    validVariationSamples++;

    // Create basis for cardinal sampling
    cbtVector3 basisX = (cbtFabs(primaryContactDir.x()) < 0.9f)
                            ? cbtVector3(1, 0, 0).cross(primaryContactDir).normalized()
                            : cbtVector3(0, 1, 0).cross(primaryContactDir).normalized();
    cbtVector3 basisY = primaryContactDir.cross(basisX).normalized();

    // Sample at 70% radius in 4 cardinal directions
    const cbtScalar sampleDist = radius * cbtScalar(0.7);
    cbtVector3 cardinalOffsets[4] = {basisX * sampleDist, -basisX * sampleDist, basisY * sampleDist,
                                     -basisY * sampleDist};

    for (int i = 0; i < 4; ++i) {
        cbtVector3 sampleQuery = sphereCenter + primaryContactDir * radius + cardinalOffsets[i];
        cbtVector3 terrainPt, terrainNormal;

        if (terrainShape->sampleWorld(terrainTransform, sampleQuery, terrainPt, terrainNormal)) {
            if (terrainNormal.length2() >= SIMD_EPSILON) {
                terrainNormal.normalize();
                variationSamples[1 + i].point = terrainPt;
                variationSamples[1 + i].normal = terrainNormal;
                variationSamples[1 + i].height = terrainPt[upAxis];
                variationSamples[1 + i].valid = true;
                validVariationSamples++;
            }
        }
    }

    // Calculate terrain variation metrics
    cbtScalar minHeight = variationSamples[0].height;
    cbtScalar maxHeight = variationSamples[0].height;
    cbtScalar avgNormalDot = 0.0f;
    int normalComparisons = 0;

    for (int i = 0; i < 5; ++i) {
        if (!variationSamples[i].valid)
            continue;

        minHeight = cbtMin(minHeight, variationSamples[i].height);
        maxHeight = cbtMax(maxHeight, variationSamples[i].height);

        for (int j = i + 1; j < 5; ++j) {
            if (!variationSamples[j].valid)
                continue;
            avgNormalDot += variationSamples[i].normal.dot(variationSamples[j].normal);
            normalComparisons++;
        }
    }

    if (normalComparisons > 0)
        avgNormalDot /= normalComparisons;

    // Terrain variation score (0 = flat, 1+ = complex)
    cbtScalar heightVariation = (maxHeight - minHeight) / (radius * cbtScalar(0.01));  // Normalized by 1% of radius
    cbtScalar normalVariation = cbtScalar(1.0) - avgNormalDot;                // 0 = identical normals, 1 = opposite
    cbtScalar terrainComplexity = heightVariation + normalVariation * 10.0f;  // Weight normals higher

    // Decide sample count based on complexity
    int numSamples;
    if (terrainComplexity < 0.5f) {
        numSamples = 9;  // Very flat: minimal sampling
    } else if (terrainComplexity < 2.0f) {
        numSamples = 25;  // Moderate: medium sampling
    } else if (terrainComplexity < 5.0f) {
        numSamples = 37;  // Complex: high sampling
    } else {
        numSamples = 49;  // Very complex: maximum sampling
    }

    // ============================================================================
    // STEP 3: Generate contact candidates with adaptive sample count
    // ============================================================================

    struct ContactCandidate {
        cbtVector3 spherePoint;
        cbtVector3 terrainPoint;
        cbtVector3 normal;
        cbtScalar penetration;
    };

    const int MAX_CANDIDATES = 49;
    ContactCandidate candidates[MAX_CANDIDATES];
    int candidateCount = 0;

    // Create hemisphere basis
    cbtVector3 basisZ = primaryContactDir;
    basisX = (cbtFabs(basisZ.x()) < 0.9f) ? cbtVector3(1, 0, 0).cross(basisZ).normalized()
                                          : cbtVector3(0, 1, 0).cross(basisZ).normalized();
    basisY = basisZ.cross(basisX).normalized();

    // Fibonacci hemisphere sampling with adaptive count
    const cbtScalar GOLDEN_ANGLE = cbtScalar(2.39996323f);

    for (int i = 0; i < numSamples && candidateCount < MAX_CANDIDATES; ++i) {
        cbtScalar phi = GOLDEN_ANGLE * i;
        cbtScalar z = cbtScalar(i) / cbtScalar(numSamples - 1);
        cbtScalar r = cbtSqrt(cbtScalar(1.0) - z * z);
        cbtScalar x = cbtCos(phi) * r;
        cbtScalar y = cbtSin(phi) * r;

        cbtVector3 localDir(x, y, z);
        cbtVector3 worldDir = basisX * localDir.x() + basisY * localDir.y() + basisZ * localDir.z();
        worldDir.normalize();

        cbtVector3 sphereSurfacePt = sphereCenter + worldDir * radius;

        cbtVector3 terrainPt, terrainNormal;
        if (!terrainShape->sampleWorld(terrainTransform, sphereSurfacePt, terrainPt, terrainNormal)) {
            continue;
        }

        if (terrainNormal.length2() < SIMD_EPSILON)
            continue;
        terrainNormal.normalize();

        cbtVector3 sphereToTerrain = terrainPt - sphereSurfacePt;
        cbtScalar signedDistance = sphereToTerrain.dot(terrainNormal);
        cbtScalar geometricPenetration = -signedDistance;
        cbtScalar finalPenetration = geometricPenetration + terrainMargin;

        if (finalPenetration < contactSlop) {
            candidates[candidateCount].spherePoint = sphereSurfacePt;
            candidates[candidateCount].terrainPoint = terrainPt;
            candidates[candidateCount].normal = terrainNormal;
            candidates[candidateCount].penetration = finalPenetration;
            candidateCount++;
        }
    }

    if (candidateCount == 0)
        return;

    // ============================================================================
    // STEP 4: Contact selection (same as before)
    // ============================================================================

    const int MAX_CONTACTS = 4;

    for (int i = 0; i < candidateCount - 1; ++i) {
        for (int j = i + 1; j < candidateCount; ++j) {
            if (candidates[j].penetration < candidates[i].penetration) {
                ContactCandidate temp = candidates[i];
                candidates[i] = candidates[j];
                candidates[j] = temp;
            }
        }
    }

    int selectedIndices[MAX_CONTACTS];
    int selectedCount = 0;
    selectedIndices[selectedCount++] = 0;

    const cbtScalar MIN_SEPARATION = radius * cbtScalar(0.35);

    for (int i = 1; i < candidateCount && selectedCount < MAX_CONTACTS; ++i) {
        if (candidates[i].penetration > 0.0f)
            break;

        bool wellSeparated = true;
        for (int j = 0; j < selectedCount; ++j) {
            cbtVector3 delta = candidates[i].spherePoint - candidates[selectedIndices[j]].spherePoint;
            if (delta.length() < MIN_SEPARATION) {
                wellSeparated = false;
                break;
            }
        }

        if (wellSeparated) {
            selectedIndices[selectedCount++] = i;
        }
    }

    // ============================================================================
    // STEP 5: Normal smoothing (adaptive radius based on complexity)
    // ============================================================================

    const cbtScalar SMOOTHING_RADIUS = (terrainComplexity < 1.0f)
                                           ? radius * cbtScalar(0.4)   // More smoothing on flat terrain
                                           : radius * cbtScalar(0.2);  // Less smoothing on complex terrain

    for (int i = 0; i < selectedCount; ++i) {
        int idx = selectedIndices[i];
        ContactCandidate& contact = candidates[idx];

        cbtVector3 smoothedNormal = contact.normal;
        cbtScalar totalWeight = 1.0f;

        for (int j = 0; j < candidateCount; ++j) {
            if (j == idx)
                continue;

            cbtScalar dist = (candidates[j].spherePoint - contact.spherePoint).length();
            if (dist < SMOOTHING_RADIUS) {
                cbtScalar weight = 1.0f - (dist / SMOOTHING_RADIUS);
                smoothedNormal += candidates[j].normal * weight;
                totalWeight += weight;
            }
        }

        if (totalWeight > SIMD_EPSILON) {
            smoothedNormal /= totalWeight;
            if (smoothedNormal.length2() > SIMD_EPSILON) {
                contact.normal = smoothedNormal.normalized();

                cbtVector3 sphereToTerrain = contact.terrainPoint - contact.spherePoint;
                cbtScalar signedDistance = sphereToTerrain.dot(contact.normal);
                contact.penetration = -signedDistance + terrainMargin;
            }
        }
    }

    // ============================================================================
    // STEP 6: Emit contacts
    // ============================================================================

    for (int i = 0; i < selectedCount; ++i) {
        const ContactCandidate& contact = candidates[selectedIndices[i]];

        if (sphereIsA) {
            resultOut->addContactPoint(contact.normal, contact.terrainPoint, contact.penetration);
        } else {
            resultOut->addContactPoint(-contact.normal, contact.spherePoint, contact.penetration);
        }
    }

    resultOut->refreshContactPoints();
}

// void cbtSphereHeightfieldAlgorithm::processCollision(const cbtCollisionObjectWrapper* wrapperA,
//                                                     const cbtCollisionObjectWrapper* wrapperB,
//                                                     const cbtDispatcherInfo& /*info*/,
//                                                     cbtManifoldResult* resultOut) {
//    if (!m_manifoldPtr)
//        return;
//    resultOut->setPersistentManifold(m_manifoldPtr);
//
//    // 1) Identify sphere vs terrain
//    bool sphereIsA = (wrapperA->getCollisionShape() == m_convex);
//    const auto* sphereWrap = sphereIsA ? wrapperA : wrapperB;
//    const auto* terrainWrap = sphereIsA ? wrapperB : wrapperA;
//    const cbtTransform& sphereTransform = sphereWrap->getWorldTransform();
//    const cbtTransform& terrainTransform = terrainWrap->getWorldTransform();
//    auto* terrainShape = static_cast<const cbtHeightfieldChronoTerrainShape*>(terrainWrap->getCollisionShape());
//    if (!terrainShape)
//        return;
//
//    // 2) Basic parameters
//    cbtScalar contactSlop = m_manifoldPtr->getContactBreakingThreshold() * 0.5f;
//    auto* sphereShape = static_cast<const cbtSphereShape*>(m_convex);
//    cbtScalar radius = sphereShape->getRadius();
//    cbtVector3 sphereCenter = sphereTransform.getOrigin();
//
//    // 3) Determine up axis for terrain
//    int upAxis = terrainShape->getUpAxis();
//    cbtVector3 worldUp(0, 0, 0);
//    worldUp[upAxis] = 1.0f;
//    cbtVector3 terrainWorldUp = terrainTransform.getBasis() * worldUp;
//
//    // 4) Set up structures for contact detection
//    struct ContactCandidate {
//        cbtVector3 terrainPoint;   // Point on terrain surface
//        cbtVector3 terrainNormal;  // Pure terrain normal
//        cbtScalar penetration;     // Penetration depth (negative when penetrating)
//    };
//
//    // Adaptive sampling based on radius
//    const int MIN_SAMPLES = 8;
//    const int MAX_SAMPLES = 32;
//    const cbtScalar RADIUS_FACTOR = 8.0f;
//    const int NUM_SAMPLES = cbtMin(MAX_SAMPLES, MIN_SAMPLES + static_cast<int>(radius * RADIUS_FACTOR));
//
//    // Normal clustering threshold
//    const cbtScalar NORMAL_ANGLE_DEG = 10.0f;
//    const cbtScalar NORMAL_ANGLE_RAD = NORMAL_ANGLE_DEG * (SIMD_PI / 180.0f);
//    const cbtScalar COS_NORMAL_ANGLE_THRESHOLD = cbtCos(NORMAL_ANGLE_RAD);
//
//    // Golden angle for Fibonacci distribution
//    const cbtScalar GOLDEN_ANGLE = cbtScalar(2.39996323f);
//
//    // Create orthogonal basis for hemisphere sampling
//    cbtVector3 basis1, basis2, downDirection;
//    downDirection = -terrainWorldUp;  // Point hemisphere downward
//
//    // Find perpendicular vectors to downDirection
//    if (std::abs(downDirection.x()) < 0.57735f) {  // 1/sqrt(3)
//        basis1 = cbtVector3(1, 0, 0).cross(downDirection);
//    } else {
//        basis1 = cbtVector3(0, 1, 0).cross(downDirection);
//    }
//    basis1.normalize();
//    basis2 = downDirection.cross(basis1);
//
//    // 5) Sample points on the SPHERE SURFACE using Fibonacci distribution
//    ContactCandidate samples[MAX_SAMPLES];
//    int validSamples = 0;
//
//    for (int i = 0; i < NUM_SAMPLES && validSamples < MAX_SAMPLES; i++) {
//        // Fibonacci hemisphere distribution on downward hemisphere
//        cbtScalar phi = GOLDEN_ANGLE * i;
//        cbtScalar y = -0.5f * (1.0f + (i / (cbtScalar)(NUM_SAMPLES - 1)));
//
//        // Convert to cartesian coordinates
//        cbtScalar r = cbtSqrt(1.0f - y * y);
//        cbtScalar x = cbtCos(phi) * r;
//        cbtScalar z = cbtSin(phi) * r;
//
//        // Calculate direction and sample point ON SPHERE SURFACE
//        cbtVector3 dir = basis1 * x + basis2 * z + downDirection * y;
//        dir.normalize();
//        cbtVector3 sphereSurfacePoint = sphereCenter + dir * radius;
//
//        // Sample terrain at this sphere surface point
//        cbtVector3 terrainPoint, terrainNormal;
//        if (terrainShape->sampleWorld(terrainTransform, sphereSurfacePoint, terrainPoint, terrainNormal)) {
//            if (terrainNormal.length2() < SIMD_EPSILON)
//                continue;
//
//            terrainNormal.normalize();
//
//            // Calculate distance from sphere center to terrain point
//            cbtVector3 sphereToTerrain = terrainPoint - sphereCenter;
//            cbtScalar distance = sphereToTerrain.length();
//            cbtScalar penetration = distance - radius;  // Negative when penetrating
//
//            // Only accept samples where terrain normal points toward sphere
//            if (sphereToTerrain.dot(terrainNormal) < 0.0f) {
//                samples[validSamples].terrainPoint = terrainPoint;
//                samples[validSamples].terrainNormal = terrainNormal;
//                samples[validSamples].penetration = penetration;
//                validSamples++;
//            }
//        }
//    }
//
//    if (validSamples == 0)
//        return;
//
//    // 6) Sort by penetration depth (deepest first)
//    const int MAX_CONTACTS = 4;
//    int sortLimit = cbtMin(validSamples, MAX_CONTACTS * 2);
//
//    for (int i = 0; i < sortLimit - 1; i++) {
//        for (int j = i + 1; j < validSamples; j++) {
//            if (samples[j].penetration < samples[i].penetration) {
//                ContactCandidate temp = samples[i];
//                samples[i] = samples[j];
//                samples[j] = temp;
//            }
//        }
//    }
//
//    // Only proceed if we have penetrating contacts
//    if (samples[0].penetration > -contactSlop) {
//        return;
//    }
//
//    // 7) Check if terrain is nearly flat
//    bool isTerrainFlat = true;
//    cbtVector3 referenceNormal = samples[0].terrainNormal;
//
//    for (int i = 1; i < cbtMin(validSamples, sortLimit); i++) {
//        cbtScalar dotProduct = referenceNormal.dot(samples[i].terrainNormal);
//        if (dotProduct < COS_NORMAL_ANGLE_THRESHOLD) {
//            isTerrainFlat = false;
//            break;
//        }
//    }
//
//    // 8) Generate final contacts
//    if (isTerrainFlat) {
//        // Single contact for flat terrain
//        if (sphereIsA) {
//            resultOut->addContactPoint(samples[0].terrainNormal, samples[0].terrainPoint, samples[0].penetration);
//        } else {
//            cbtVector3 pointOnSphere = sphereCenter - samples[0].terrainNormal * radius;
//            resultOut->addContactPoint(-samples[0].terrainNormal, pointOnSphere, samples[0].penetration);
//        }
//    } else {
//        // Multiple contacts for varied terrain
//        int selectedContacts[MAX_CONTACTS] = {0, -1, -1, -1};
//        int numSelectedContacts = 1;
//
//        // Add diverse contacts based on terrain normal differences
//        for (int i = 1; i < cbtMin(validSamples, sortLimit) && numSelectedContacts < MAX_CONTACTS; i++) {
//            if (samples[i].penetration > -contactSlop)
//                continue;
//
//            bool isUnique = true;
//            for (int j = 0; j < numSelectedContacts; j++) {
//                cbtScalar dotProduct = samples[i].terrainNormal.dot(samples[selectedContacts[j]].terrainNormal);
//                if (dotProduct > COS_NORMAL_ANGLE_THRESHOLD) {
//                    isUnique = false;
//                    break;
//                }
//            }
//
//            if (isUnique) {
//                selectedContacts[numSelectedContacts++] = i;
//            }
//        }
//
//        // Add contacts to manifold
//        for (int i = 0; i < numSelectedContacts; i++) {
//            int idx = selectedContacts[i];
//            if (sphereIsA) {
//                resultOut->addContactPoint(samples[idx].terrainNormal, samples[idx].terrainPoint,
//                                           samples[idx].penetration);
//            } else {
//                cbtVector3 pointOnSphere = sphereCenter - samples[idx].terrainNormal * radius;
//                resultOut->addContactPoint(-samples[idx].terrainNormal, pointOnSphere, samples[idx].penetration);
//            }
//        }
//    }
//
//    resultOut->refreshContactPoints();
//}
//


// TODO- is this handling the sphere correclty if the sphere comes up against a sharp incline? do we sample about the hemispheres?

cbtScalar cbtSphereHeightfieldAlgorithm::calculateTimeOfImpact(cbtCollisionObject* body0,
                                                               cbtCollisionObject* body1,
                                                               const cbtDispatcherInfo& dispatchInfo,
                                                               cbtManifoldResult* resultOut) {
    bool terrainIsA = (body0->getCollisionShape()->getShapeType() == TERRAIN_SHAPE_PROXYTYPE);
    auto* terrain = terrainIsA ? static_cast<const cbtHeightfieldChronoTerrainShape*>(body0->getCollisionShape())
                               : static_cast<const cbtHeightfieldChronoTerrainShape*>(body1->getCollisionShape());
    auto* sphere = terrainIsA ? static_cast<const cbtSphereShape*>(body1->getCollisionShape())
                              : static_cast<const cbtSphereShape*>(body0->getCollisionShape());

    const cbtTransform& fromTransform = terrainIsA ? body1->getWorldTransform() : body0->getWorldTransform();
    const cbtTransform& toTransform =
        terrainIsA ? body1->getInterpolationWorldTransform() : body0->getInterpolationWorldTransform();
    const cbtTransform& terrainTransform = terrainIsA ? body0->getWorldTransform() : body1->getWorldTransform();

    cbtVector3 motion = toTransform.getOrigin() - fromTransform.getOrigin();
    if (motion.length2() < SIMD_EPSILON)
        return cbtScalar(1.);

    int upAxis = terrain->getUpAxis();
    cbtScalar radius = sphere->getRadius();                         // Sphere radius (with margin)
    cbtScalar totalMargin = radius + terrain->getMargin() + 0.01f;  // Extra tolerance

    // Sample bottom point + cross directions for better edge detection
    cbtScalar minTOI = 1.f;

    // Directions: bottom (-up) + 4 cross (in terrain plane)
    cbtVector3 terrainUp(0, 0, 0);
    terrainUp[upAxis] = 1.f;
    cbtVector3 supportDirLocal = fromTransform.getBasis().transpose() * (-terrainUp);  // To sphere local

    // Sample points: bottom + 4 offsets
    std::vector<cbtVector3> sampleDirs = {
        supportDirLocal,                            // Bottom
        cbtVector3(1, 0, 0), cbtVector3(-1, 0, 0),  // Cross X
        cbtVector3(0, 1, 0), cbtVector3(0, -1, 0)   // Cross Y
    };

    for (const auto& dirLocal : sampleDirs) {
        cbtVector3 supportVertexLocal = sphere->localGetSupportingVertexWithoutMargin(dirLocal.normalized());
        cbtScalar low = 0.f, high = 1.f;
        for (int iter = 0; iter < 25; ++iter) {  // Binary search for TOI
            cbtScalar mid = 0.5f * (low + high);
            cbtVector3 interpolatedPos = fromTransform.getOrigin() + motion * mid;
            cbtVector3 supportVertexWorld = fromTransform * supportVertexLocal + motion * mid;

            // Sample terrain at support point (world)
            cbtVector3 terrainSurface, terrainNormal;
            if (!terrain->sampleWorld(terrainTransform, supportVertexWorld, terrainSurface, terrainNormal))
                continue;  // Outside bounds

            // Signed distance: Positive above, negative penetrating
            cbtScalar distToSurface = (supportVertexWorld - terrainSurface).dot(terrainNormal) - totalMargin;
            bool penetrating = distToSurface <= 0.f;
            (penetrating ? high : low) = mid;
        }
        minTOI = cbtMin(minTOI, high);  // Take min TOI across samples
    }

    return minTOI;
}
}  // namespace chrono

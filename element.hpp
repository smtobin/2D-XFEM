#pragma once

#include "common.hpp"

/** Represents a planar linear triangular element with Neo-Hookean energy.
 * Does not own any simulation data (e.g. the nodal positions), merely does computations on current state.
 */
struct TriangleElement
{
    /** Pointers to the nodal positions that make up the element.
     * Assumes a counter-clockwise ordering.
     */
    const Vec2r* nodes[3];

    /** Pointers to phantom nodes when the element is cut.
     * Phantom node i is associated with node i.
     */
    const Vec2r* phantom_nodes[3];

    /** Undeformed element edge matrix. Used to calculate deformation gradient. */
    Mat2r Dm_inv;

    /** Undeformed element area. */
    Real A0;

    /** Material properties, specified as Lame parameters */
    Real lambda;
    Real mu;

    /** Create an element from 3 nodal pointers.
     * Assumes that the element is in its rest state, and computes the rest element edge matrix.
     */
    TriangleElement(const Vec2r* node1, const Vec2r* node2, const Vec2r* node3, Real lambda_, Real mu_);

    /** Computes and adds internal force to provided force vectors. */
    Vec6r internalForce() const;

    /** Computes the element stiffness matrix */
    Mat6r elementStiffnessMatrix() const;
};
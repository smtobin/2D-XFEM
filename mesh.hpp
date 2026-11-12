
#pragma once

#include "common.hpp"

#include <vector>

/** Represents an element in the mesh */
struct TriangleElement
{
    unsigned nodes[3];
    Mat2r Dm_inv;
    Real A0;
};

/** Stores and does computations on a collection of triangle elements. */
struct Mesh
{
    /** Element storage */
    std::vector<TriangleElement> elements;

    /** Node storage */
    std::vector<Vec2r> nodes;

    /** Phantom node storage */
    std::vector<Vec2r> phantom_nodes;

    /** Material properties - Lame parameters */
    Real lambda;
    Real mu;

    /** Creates a square mesh, with bottom left corner at (0,0). */
    Mesh(Real side_length, unsigned num_elements_per_side, Real lambda, Real mu);

    void _createElement(unsigned n1, unsigned n2, unsigned n3);
    Vec6r _elementInternalForce(unsigned e_idx) const;
    Mat6r _elementStiffnessMatrix(unsigned e_idx) const;

    /** Assembles global stiffness matrix */
    MatXr assembleStiffnessMatrix() const;

    /** Assembles global internal force vector */
    VecXr assembleInternalForces() const;

};
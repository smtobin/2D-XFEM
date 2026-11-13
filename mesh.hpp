
#pragma once

#include "common.hpp"

#include <vector>
#include <unordered_set>

/** Represents an edge in the mesh.
 * n1 < n2
 */
struct Edge
{
    unsigned n1;
    unsigned n2;

    Edge(unsigned n1_, unsigned n2_)
        : n1(std::min(n1_, n2_)), n2(std::max(n1_, n2_))
    {}

    bool operator==(const Edge& other) const
    {
        return (other.n1 == n1) && (other.n2 == n2);
    }
};

struct EdgeHash
{
    size_t operator()(const Edge& e) const {
        auto h1 = std::hash<unsigned>{}(e.n1);
        auto h2 = std::hash<unsigned>{}(e.n2);
        // Better mixing than simple XOR
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

/** Represents a fragment of a triangular element
 * This generalizes both full elements and cut fragments.
 */
struct TriangleElementFragment
{
    unsigned nodes[3];
    std::vector<Vec2r> region_barys;
    Real A0;
};

/** Represents an element in the mesh. */
struct TriangleElement
{
    unsigned nodes[3];
    Mat2r Dm_inv;
    Real A0;
    std::vector<TriangleElementFragment> fragments;
};

/** Represents a cut along a line in the mesh.
 * Taken to be in the undeformed configuration.
 */
struct Cut
{
    // for now, simplest case: just a straight line
    Vec2r p1;
    Vec2r p2;
};


/** Stores and does computations on a collection of triangle elements. */
struct Mesh
{
    /** Element storage */
    std::vector<TriangleElement> elements;

    /** Node storage */
    unsigned initial_num_nodes;
    std::vector<Vec2r> nodes;
    std::vector<Vec2r> initial_nodes;

    /** Phantom node storage */
    std::vector<Vec2r> phantom_nodes;

    /** Cut storage */
    std::vector<Cut> cuts;

    /** List of edges in the mesh
     * These are tested for intersections with cuts.
     */
    std::unordered_set<Edge, EdgeHash> edges;

    /** Material properties - Lame parameters */
    Real lambda;
    Real mu;

    /** Creates a square mesh, with bottom left corner at (0,0). */
    Mesh(Real side_length, unsigned num_elements_per_side, Real lambda, Real mu);

    void _createElement(unsigned n1, unsigned n2, unsigned n3);
    void _createEdgesList();

    Vec6r _elementFragmentInternalForce(const TriangleElement& elem, const TriangleElementFragment& f) const;
    Mat6r _elementFragmentStiffnessMatrix(const TriangleElement& elem, const TriangleElementFragment& f) const;

    /** Assembles global stiffness matrix */
    MatXr assembleStiffnessMatrix() const;

    /** Assembles global internal force vector */
    VecXr assembleInternalForces() const;

    /** Adds a cut to the mesh. */
    void addCut(const Cut& cut);

};
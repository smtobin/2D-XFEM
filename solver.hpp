#pragma once

#include "mesh.hpp"

#include <unordered_set>

/** Takes a mesh and load conditions, and solves the quasistatic problem. */
struct Solver
{
    /** The mesh */
    Mesh* mesh;

    /** External nodal forces */
    std::vector<Vec2r> external_forces;

    /** Tracks which nodes have displacement BCs.
     */
    std::unordered_set<unsigned> nodes_w_displacement_BC;

    /** Initialize the solver with a mesh. */
    Solver(Mesh* mesh);

    /** Sets the external point force for a node. */
    void setNodeExternalForce(unsigned n_idx, const Vec2r& f) { external_forces[n_idx] = f; }

    /** Prescribes a displacement boundary condition. */
    void setNodePosition(unsigned n_idx, const Vec2r& p) { mesh->nodes[n_idx] = p; nodes_w_displacement_BC.insert(n_idx); }

    /** Prescribes a node to be fixed (i.e. zero displacement) */
    void setNodeFixed(unsigned n_idx) { nodes_w_displacement_BC.insert(n_idx); }

    /** Solves the nonlinear equations with Newton-Raphson */
    void solve();
};
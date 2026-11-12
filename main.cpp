#include "common.hpp"
#include "mesh.hpp"
#include "solver.hpp"

#include <iostream>

int main()
{
    // arbitrary material parameters
    Real E = 100;
    Real nu = 0.4;
    Real lambda = E*nu / ( (1+nu)*(1-2*nu) );
    Real mu = E / (2*(1+nu));

    // create mesh
    unsigned elems_per_side = 3;
    Mesh mesh(1, elems_per_side, lambda, mu);

    // create solver
    Solver solver(&mesh);

    // fix bottom edge
    for (unsigned i = 0; i < elems_per_side+1; i++)
    {
        solver.setNodeFixed(i);
    }

    // apply force to top-left node
    solver.setNodeExternalForce(mesh.nodes.size()-1, Vec2r(1,0));

    solver.solve();

    for (unsigned i = 0; i < mesh.nodes.size(); i++)
    {
        std::cout << "Node " << i << ": " << mesh.nodes[i].transpose() << std::endl;
    }
}
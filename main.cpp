#include "common.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "visualization.hpp"

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
    unsigned topright = mesh.nodes.size()-1;

    // fix left edge
    for (unsigned i = 0; i < elems_per_side+1; i++)
    {
        mesh.setNodeFixed(i);
    }

    Mesh initial_mesh = mesh;

    Cut cut;
    cut.p1 = Vec2r(-1, 0.5);
    cut.p2 = Vec2r(5, 0.5);
    mesh.addCut(cut);

    

    // create solver
    Solver solver(&mesh);

    

    // apply force to top-right node
    solver.setNodeExternalForce(topright, Vec2r(0,1));

    solver.solve();

    for (unsigned i = 0; i < mesh.nodes.size(); i++)
    {
        std::cout << "Node " << i << ": " << mesh.nodes[i].transpose() << std::endl;
    }

    std::vector<const Mesh*> meshes = {&initial_mesh, &mesh};
    std::vector<std::string> mesh_names = {"Undeformed", "Deformed"};
    Visualization viz(meshes, mesh_names);
    viz.show();
}
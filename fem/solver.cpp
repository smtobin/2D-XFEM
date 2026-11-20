#include "solver.hpp"

Solver::Solver(Mesh* mesh_)
    : mesh(mesh_)
    , external_forces(mesh->nodes.size(), Vec2r::Zero())    // initialize all external forces to 0
{

}

void Solver::solve()
{
    // assemble external force vector
    VecXr f_ext = VecXr::Zero(mesh->nodes.size()*2);
    for (unsigned n_idx = 0; n_idx < mesh->nodes.size(); n_idx++)
    {
        f_ext.block<2,1>(n_idx*2,0) = external_forces[n_idx];
    }


    unsigned MAX_ITER = 20;
    for (unsigned iter = 0; iter < MAX_ITER; iter++)
    {
        // compute global stiffness matrix and internal force vector
        MatXr K = mesh->assembleStiffnessMatrix();
        VecXr f_int = mesh->assembleInternalForces();

        // zero out rows/columns corresponding to nodes with displacement BCs
        for (unsigned n_idx : mesh->nodes_w_displacement_BCs)
        {
            // zero out rows/cols in stiffness matrix
            K.row(2*n_idx) = VecXr::Zero(2*mesh->nodes.size());
            K.row(2*n_idx+1) = VecXr::Zero(2*mesh->nodes.size());
            K.col(2*n_idx) = VecXr::Zero(2*mesh->nodes.size());
            K.col(2*n_idx+1) = VecXr::Zero(2*mesh->nodes.size());
            K(2*n_idx, 2*n_idx) = 1;
            K(2*n_idx+1, 2*n_idx+1) = 1;
            
            // zero out forces
            f_int[2*n_idx] = 0;
            f_int[2*n_idx+1] = 0;
            f_ext[2*n_idx] = 0;
            f_ext[2*n_idx+1] = 0;
        }

        // compute residual
        VecXr res = f_ext - f_int;

        // compute position update
        // VecXr dx = K.llt().solve(res);

        // basic diagonal preconditioning
        VecXr d = K.diagonal().array().sqrt().inverse();
        // MatXr S = d.asDiagonal();

        // scale system
        MatXr Kscaled = d.asDiagonal() * K * d.asDiagonal();
        VecXr bscaled = d.asDiagonal() * res;

        // solve
        VecXr y = Kscaled.ldlt().solve(bscaled);

        Real rcond_K = K.ldlt().rcond();
        Real rcond_Kscaled = Kscaled.ldlt().rcond();
        std::cout << "K rcond: " << rcond_K << std::endl;
        std::cout << "Kscaled rcond: " << rcond_Kscaled << std::endl;

        // Transform solution back
        VecXr dx = d.asDiagonal() * y;

        // update positions
        for (unsigned n_idx = 0; n_idx < mesh->nodes.size(); n_idx++)
        {
            mesh->nodes[n_idx] += dx.block<2,1>(n_idx*2, 0);
        }
    }
}
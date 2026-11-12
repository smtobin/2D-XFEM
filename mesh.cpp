#include "mesh.hpp"

Mesh::Mesh(Real side_length, unsigned num_elements_per_side, Real lambda_, Real mu_)
    : lambda(lambda_)
    , mu(mu_)
{
    Real elem_side_length = side_length / num_elements_per_side;
    // create square grid of nodes
    nodes.reserve((num_elements_per_side+1)*(num_elements_per_side+1));
    for (unsigned i = 0 ; i < num_elements_per_side+1; i++)
    {
        for (unsigned j = 0; j < num_elements_per_side+1; j++)
        {
            nodes.emplace_back(i*elem_side_length, j*elem_side_length);
        }
    }

    // create elements for square mesh
    // iterate through the grid and create two triangles per square
    elements.reserve(num_elements_per_side*num_elements_per_side);
    for (unsigned i = 0; i < num_elements_per_side; i++)
    {
        for (unsigned j = 0; j < num_elements_per_side; j++)
        {
            unsigned s1 = j*(num_elements_per_side+1) + i;
            unsigned s2 = j*(num_elements_per_side+1) + i + 1;
            unsigned s3 = (j+1)*(num_elements_per_side+1) + i + 1;
            unsigned s4 = (j+1)*(num_elements_per_side+1) + i;

            _createElement(s1, s2, s3);
            _createElement(s1, s3, s4);
        }
    }
}

void Mesh::_createElement(unsigned n1, unsigned n2, unsigned n3)
{
    TriangleElement new_elem;
    new_elem.nodes[0] = n1;
    new_elem.nodes[1] = n2;
    new_elem.nodes[2] = n3;
    
    // compute the undeformed edge matrix for future deformation gradient computation
    Mat2r Dm;
    Dm.col(0) = nodes[n2] - nodes[n1];
    Dm.col(1) = nodes[n3] - nodes[n1];
    new_elem.Dm_inv = Dm.inverse();

    // element undeformed area
    new_elem.A0 = 0.5*std::abs(Dm.determinant());

    elements.push_back(std::move(new_elem));
}


MatXr Mesh::assembleStiffnessMatrix() const
{
    // iterate through elements, compute elemental stiffness matrices, and scatter
    MatXr K = MatXr::Zero(nodes.size()*2, nodes.size()*2);
    for (unsigned e_idx = 0; e_idx < elements.size(); e_idx++)
    {
        const TriangleElement& elem = elements[e_idx];

        Mat6r Ke = _elementStiffnessMatrix(e_idx);
        for (unsigned i = 0; i < 3; i++)
        {
            for (unsigned j = 0; j < 3; j++)
            {
                K.block<2,2>(2*elem.nodes[i], 2*elem.nodes[j]) += Ke.block<2,2>(2*i, 2*j);
            }
        }
    }

    return K;
}

VecXr Mesh::assembleInternalForces() const
{
    // iterate through elements, compute internal forces, and scatter
    VecXr f = VecXr::Zero(nodes.size()*2);
    for (unsigned e_idx = 0; e_idx < elements.size(); e_idx++)
    {
        const TriangleElement& elem = elements[e_idx];

        Vec6r fe = _elementInternalForce(e_idx);
        for (unsigned i = 0; i < 3; i++)
        {
            f.block<2,1>(2*elem.nodes[i], 0) += fe.block<2,1>(2*i, 0);
        }
    }

    return f;
}


Vec6r Mesh::_elementInternalForce(unsigned e_idx) const
{
    const TriangleElement& elem = elements[e_idx];

    // get current nodal positions
    const Vec2r& p1 = nodes[elem.nodes[0]];
    const Vec2r& p2 = nodes[elem.nodes[1]];
    const Vec2r& p3 = nodes[elem.nodes[2]];

    // calculate deformation gradient
    Mat2r Ds;
    Ds.col(0) = p2 - p1;
    Ds.col(1) = p3 - p1;
    Mat2r F = Ds * elem.Dm_inv;

    // calculate first Piola-Kirchhoff stress (dE/dF)
    Real J = F.determinant();
    Mat2r P = lambda*(J - 1)*J*F.transpose().inverse() + mu*F;
    
    // compute internal forces
    Mat2r G = P*elem.Dm_inv.transpose();
    Vec6r f;
    f.block<2,1>(0,0) = elem.A0 * (G.col(0) + G.col(1));
    f.block<2,1>(2,0) = -elem.A0 * G.col(0);
    f.block<2,1>(4,0) = -elem.A0 * G.col(1); 
    return f;
}

Mat6r Mesh::_elementStiffnessMatrix(unsigned e_idx) const
{
    const TriangleElement& elem = elements[e_idx];

    // get current nodal positions
    const Vec2r& p1 = nodes[elem.nodes[0]];
    const Vec2r& p2 = nodes[elem.nodes[1]];
    const Vec2r& p3 = nodes[elem.nodes[2]];

    // calculate deformation gradient
    Mat2r Ds;
    Ds.col(0) = p2 - p1;
    Ds.col(1) = p3 - p1;
    Mat2r F = Ds * elem.Dm_inv;
    Mat2r F_invT = F.transpose().inverse();
    Real J = F.determinant();

    // shape function gradient matrix
    Eigen::Matrix<Real, 3, 2> gradN;
    gradN.block<2,2>(1,0) = elem.Dm_inv;
    gradN.row(0) = -gradN.row(1) - gradN.row(2);

    // strain-displacement matrix
    Eigen::Matrix<Real, 4, 6> B = Eigen::Matrix<Real, 4, 6>::Zero();
    for (unsigned i = 0; i < 3; i++)
    {
        B(0, 2*i) = gradN(i, 0);
        B(1, 2*i) = gradN(i, 1);
        B(2, 2*i+1) = gradN(i, 0);
        B(3, 2*i+1) = gradN(i, 1);
    }

    // material tangent
    Mat4r C;
    for (unsigned i = 0; i < 2; i++)
    {
        for (unsigned A = 0; A < 2; A++)
        {
            unsigned p = 2*i + A;

            for (unsigned j = 0; j < 2; j++)
            {
                for (unsigned B = 0; B < 2; B++)
                {
                    unsigned q = 2*j + B;

                    C(p,q) = lambda * (2*J - 1) * J * F_invT(i, A) * F_invT(j, B)
                            - lambda * (J - 1) * J * F_invT(i, B) * F_invT(j, A)
                            + mu * (i == j) * (A == B);
                }
            }
        }
    }

    Mat6r Ke = elem.A0 * B.transpose() * C * B;
    return Ke;
}
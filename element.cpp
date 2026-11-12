#include "element.hpp"

TriangleElement::TriangleElement(const Vec2r* node1, const Vec2r* node2, const Vec2r* node3, Real lambda_, Real mu_)
    : nodes{node1, node2, node3}
    , phantom_nodes{nullptr, nullptr, nullptr}
    , lambda(lambda_)
    , mu(mu_)
{
    // assume the current state is the undeformed state
    // compute the undeformed edge matrix for future deformation gradient computation
    Mat2r Dm;
    Dm.col(0) = *node2 - *node1;
    Dm.col(1) = *node3 - *node1;
    Dm_inv = Dm.inverse();

    // element undeformed area
    A0 = 0.5*Dm.determinant();
}

Vec6r TriangleElement::internalForce() const
{
    // get current nodal positions
    const Vec2r& p1 = *nodes[0];
    const Vec2r& p2 = *nodes[1];
    const Vec2r& p3 = *nodes[2];

    // calculate deformation gradient
    Mat2r Ds;
    Ds.col(0) = p2 - p1;
    Ds.col(1) = p3 - p1;
    Mat2r F = Ds * Dm_inv;

    // calculate first Piola-Kirchhoff stress (dE/dF)
    Real J = F.determinant();
    Mat2r P = lambda*(J - 1)*J*F.transpose().inverse() + mu*F;
    
    // compute internal forces
    Mat2r G = P*Dm_inv.transpose();
    Vec6r f;
    f.block<2,1>(0,0) = A0 * (G.col(0) + G.col(1));
    f.block<2,1>(2,0) = -A0 * G.col(0);
    f.block<2,1>(4,0) = -A0 * G.col(1); 
    return f;
}

Mat6r TriangleElement::elementStiffnessMatrix() const
{
    // get current nodal positions
    const Vec2r& p1 = *nodes[0];
    const Vec2r& p2 = *nodes[1];
    const Vec2r& p3 = *nodes[2];

    // calculate deformation gradient
    Mat2r Ds;
    Ds.col(0) = p2 - p1;
    Ds.col(1) = p3 - p1;
    Mat2r F = Ds * Dm_inv;
    Mat2r F_invT = F.transpose().inverse();
    Real J = F.determinant();

    // shape function gradient matrix
    Eigen::Matrix<Real, 3, 2> gradN;
    gradN.block<2,2>(1,0) = Dm_inv;
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

    Mat6r Ke = A0 * B.transpose() * C * B;
    return Ke;
}
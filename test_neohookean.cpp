#include "common.hpp"
#include "element.hpp"


/** Test internal forces and stiffness matrix with finite differences */
int main()
{

    Vec2r n1 = Vec2r::Random();
    Vec2r n2 = Vec2r::Random();
    Vec2r n3 = Vec2r::Random();

    Real lambda = 1000;
    Real mu = 500;

    TriangleElement element(&n1, &n2, &n3, lambda, mu);

    // perturb state
    n1 += Vec2r::Random();
    n2 += Vec2r::Random();
    n3 += Vec2r::Random();

    // evaluate internal force
    Vec6r f = element.internalForce();

    Mat6r Ke_FD;
    Real delta = 1e-6;
    for (unsigned i = 0; i < 6; i++)
    {
        if (i/2 == 0)   n1[i%2] += delta;
        else if (i/2 == 1) n2[i%2] += delta;
        else if (i/2 == 2) n3[i%2] += delta;

        Vec6r fp = element.internalForce();

        if (i/2 == 0)   n1[i%2] -= 2*delta;
        else if (i/2 == 1) n2[i%2] -= 2*delta;
        else if (i/2 == 2) n3[i%2] -= 2*delta;

        Vec6r fm = element.internalForce();

        Ke_FD.col(i) = -(fp - fm) / (2*delta);  // use negative because Ke is d^2E/dx^2 and f = -dE/dx

        if (i/2 == 0)   n1[i%2] += delta;
        else if (i/2 == 1) n2[i%2] += delta;
        else if (i/2 == 2) n3[i%2] += delta;
    }

    Mat6r Ke = element.elementStiffnessMatrix();

    Mat6r diff = Ke_FD - Ke;
    std::cout << "Ke:\n" << Ke << std::endl;
    std::cout << "Ke diff:\n" << diff << std::endl;

}
#include "mesh.hpp"

#include <unordered_map>

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

    // save initial mesh state
    initial_nodes = nodes;

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

    // create the edges list
    _createEdgesList();
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

    // create a single fragment for the element
    TriangleElementFragment f;
    f.A0 = new_elem.A0;
    f.nodes[0] = new_elem.nodes[0];
    f.nodes[1] = new_elem.nodes[1];
    f.nodes[2] = new_elem.nodes[2];
    f.region_barys = { Vec2r(0,0), Vec2r(1,0), Vec2r(0,1) };
    new_elem.fragments.push_back(f);

    elements.push_back(std::move(new_elem));
}

void Mesh::_createEdgesList()
{
    edges.clear();
    for (const auto& elem : elements)
    {
        edges.emplace(elem.nodes[0], elem.nodes[1]);
        edges.emplace(elem.nodes[0], elem.nodes[2]);
        edges.emplace(elem.nodes[1], elem.nodes[2]);
    }
}


MatXr Mesh::assembleStiffnessMatrix() const
{
    // iterate through elements, compute elemental stiffness matrices, and scatter
    MatXr K = MatXr::Zero(nodes.size()*2, nodes.size()*2);
    for (unsigned e_idx = 0; e_idx < elements.size(); e_idx++)
    {
        const TriangleElement& elem = elements[e_idx];

        /** TODO: iterate over fragments, not just elements */
        for (const auto& frag : elem.fragments)
        {
            Mat6r Ke = _elementFragmentStiffnessMatrix(elem, frag);
            for (unsigned i = 0; i < 3; i++)
            {
                for (unsigned j = 0; j < 3; j++)
                {
                    K.block<2,2>(2*frag.nodes[i], 2*frag.nodes[j]) += Ke.block<2,2>(2*i, 2*j);
                }
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

        for (const auto& frag : elem.fragments)
        {
            Vec6r fe = _elementFragmentInternalForce(elem, frag);
            for (unsigned i = 0; i < 3; i++)
            {
                f.block<2,1>(2*frag.nodes[i], 0) += fe.block<2,1>(2*i, 0);
            }
        }
        
    }

    return f;
}


Vec6r Mesh::_elementFragmentInternalForce(const TriangleElement& elem, const TriangleElementFragment& frag) const
{
    // get current nodal positions
    const Vec2r& p1 = nodes[frag.nodes[0]];
    const Vec2r& p2 = nodes[frag.nodes[1]];
    const Vec2r& p3 = nodes[frag.nodes[2]];

    // calculate deformation gradient
    Mat2r Ds;
    Ds.col(0) = p2 - p1;
    Ds.col(1) = p3 - p1;
    Mat2r F = Ds * elem.Dm_inv;

    // calculate first Piola-Kirchhoff stress (dE/dF)
    Real J = F.determinant();
    Mat2r P = lambda*(J - 1 - mu/lambda)*J*F.transpose().inverse() + mu*F;
    
    // compute internal forces
    Mat2r G = P*elem.Dm_inv.transpose();
    Vec6r force;
    force.block<2,1>(0,0) = -frag.A0 * (G.col(0) + G.col(1));
    force.block<2,1>(2,0) = frag.A0 * G.col(0);
    force.block<2,1>(4,0) = frag.A0 * G.col(1); 
    return force;
}

Mat6r Mesh::_elementFragmentStiffnessMatrix(const TriangleElement& elem, const TriangleElementFragment& frag) const
{
    // get current nodal positions
    const Vec2r& p1 = nodes[frag.nodes[0]];
    const Vec2r& p2 = nodes[frag.nodes[1]];
    const Vec2r& p3 = nodes[frag.nodes[2]];

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

                    C(p,q) = lambda * (2*J - 1 - mu/lambda) * J * F_invT(i, A) * F_invT(j, B)
                            - lambda * (J - 1 - mu/lambda) * J * F_invT(i, B) * F_invT(j, A)
                            + mu * (i == j) * (A == B);
                }
            }
        }
    }

    Mat6r Ke = frag.A0 * B.transpose() * C * B;
    return Ke;
}


Real _polygonArea(const std::vector<Vec2r>& poly)
{
    Real a = 0;
    for (size_t i = 0; i < poly.size(); i++) 
    {
        const Vec2r& p = poly[i];
        const Vec2r& q = poly[(i + 1) % poly.size()];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;   // positive for CCW
}

void Mesh::addCut(const Cut& cut)
{
    cuts.push_back(cut);

    const Vec2r d = cut.p2 - cut.p1;
    const Real len2 = d.dot(d);
    const Real len = std::sqrt(len2);
    const Vec2r n(d[1] / len, -d[0] / len);   // unit normal
    const Real eps = 1e-9 * len;

    // loop through nodes and compute signed distance of node relative to cut
    std::vector<Real> sd(initial_nodes.size());
    for (unsigned n_idx = 0; n_idx < sd.size(); n_idx++)
    {
        Real sd_i = n.dot(initial_nodes[n_idx] - cut.p1);
        if (std::abs(sd_i) < eps) sd_i = eps;
        sd[n_idx] = sd_i;
    }

    // loop through edges and find cut edges
    // if cut, calculate interpolation parameter where intersection is
    std::unordered_map<Edge, Real, EdgeHash> cut_edges;
    for (const Edge& e : edges)
    {
        Real sd1 = sd[e.n1];
        Real sd2 = sd[e.n2];

        // both nodes on the same "side" of the cut, continue
        if ((sd1 > 0) == (sd2 > 0))
            continue;

        // get location of cut as interpolation parameter in (0,1)
        // t defined as going from lower index (n1) to higher index (n2)
        Real t = sd1 / (sd1 - sd2);

        // ensure that the cut makes it through this edge
        Vec2r p = initial_nodes[e.n1] + t * (initial_nodes[e.n2] - initial_nodes[e.n1]);
        Real tau = (p - cut.p1).dot(d) / len2;
        if (tau < 0 || tau > 1)
            continue;

        cut_edges[e] = t;
    }

    // create phantom nodes lazily, one per original node
    std::vector<int> phantom_of(sd.size(), -1);
    auto get_phantom = [&](unsigned i) -> unsigned {
        if (phantom_of[i] < 0)
        {
            phantom_of[i] = nodes.size();
            nodes.push_back(nodes[i]);
            initial_nodes.push_back(initial_nodes[i]);
        }
        return static_cast<unsigned>(phantom_of[i]);
    };

    // barycentric coordinates of vertices
    const Vec2r V[3] = { Vec2r(0, 0), Vec2r(1, 0), Vec2r(0, 1) };

    // loop through elements and identify which ones have been fully cut (2 edges intersected by the cut)
    // for elements that have been fully cut, compute the element fragments
    for (auto& el : elements)
    {
        bool is_cut[3];
        Real t_edge[3];
        unsigned num_cut = 0;
        for (unsigned k = 0; k < 3; k++)
        {
            // edge k in the element goes from node k to (k+1)%3 
            unsigned a = el.nodes[k];
            unsigned b = el.nodes[(k+1)%3];

            // check if this edge has been cut
            auto it = cut_edges.find(Edge(a,b));
            is_cut[k] = (it != cut_edges.end());

            // if so, store the intersection point
            if (is_cut[k])
            {
                num_cut++;
                t_edge[k] = (a < b) ? it->second : 1 - it->second;
            }
        }
        
        if (num_cut < 2)
            continue;

        assert(num_cut == 2);

        // the uncut edge
        unsigned u = !is_cut[0] ? 0 : (!is_cut[1] ? 1 : 2);
        // the lone vertex
        int L = (u + 2) % 3;
        // the other vertices
        int A = (L + 1) % 3;
        int B = (L + 2) % 3;

        // edge L goes L->A, edge B goes B->L
        Vec2r P_LA = V[L] + t_edge[L] * (V[A] - V[L]);
        Vec2r P_BL = V[B] + t_edge[B] * (V[L] - V[B]);

        // whether or not L is on "positive" side of cut
        bool L_pos = sd[el.nodes[L]] > 0;

        // helper to make element fragment
        auto make_fragment = [&](bool side_pos, std::vector<Vec2r> poly)
        {
            TriangleElementFragment frag;
            for (unsigned v = 0; v < 3; v++)
            {
                unsigned orig = el.nodes[v];
                bool v_pos = sd[orig] > 0;
                // if vertex i in the element is on the same side of the cut, use the original node
                // otherwise, use the phantom node instead of the original node
                frag.nodes[v] = (v_pos == side_pos) ? orig : get_phantom(orig);    // creates a phantom if not created already
            }

            // check for fixed edges - if the newly created fragment overlaps a fixed edge, we must constrain the newly created phantom node
            for (unsigned k = 0; k < 3; k++)
            {
                unsigned a = k;
                unsigned b = (k+1)%3;
                if (nodes_w_displacement_BCs.count(el.nodes[a]) + nodes_w_displacement_BCs.count(el.nodes[b]) < 2)
                    continue;
                bool on_edge = (poly.size() == 3 && k != u)     // for the triangle fragment, it is on the edge unless k is the uncut edge
                                || (poly.size() == 4);          // the quad fragment is on all 3 edges of the original triangle
                
                // the fragment is on the edge, and it is fixed
                if (on_edge)
                {
                    // mark nodes as being fixed
                    nodes_w_displacement_BCs.insert(frag.nodes[a]);
                    nodes_w_displacement_BCs.insert(frag.nodes[b]);
                    // set node positions to their prescribed values
                    nodes[frag.nodes[a]] = nodes[el.nodes[a]];
                    nodes[frag.nodes[b]] = nodes[el.nodes[b]];
                }
            }

            frag.A0 = 2 * el.A0 * _polygonArea(poly);
            frag.region_barys = std::move(poly);
            return frag;
        };

        el.fragments.clear();
        // create the two element fragments created by this cut through the element
        el.fragments.push_back(make_fragment(L_pos, { V[L], P_LA, P_BL } ));
        el.fragments.push_back(make_fragment(!L_pos, { P_LA, V[A], V[B], P_BL }));

        assert(std::abs(el.fragments[0].A0 + el.fragments[1].A0 - el.A0) < 1e-9 * el.A0);
    }
}
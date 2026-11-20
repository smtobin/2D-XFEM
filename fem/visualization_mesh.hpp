#pragma once

#include "common.hpp"
#include "mesh.hpp"

#include <vector>

struct VisualizationMesh
{
    struct VertexSource { unsigned elem, frag; Vec2r bary; };
    std::vector<VertexSource> sources;
    std::vector<std::vector<unsigned>> cells;   // vertex ids per fragment polygon
    std::vector<unsigned> cell_elem;            // for cell data (element id, etc.)

    void rebuildTopology(const Mesh& m)
    {
        sources.clear(); cells.clear(); cell_elem.clear();
        for (unsigned e = 0; e < m.elements.size(); e++)
            for (unsigned f = 0; f < m.elements[e].fragments.size(); f++)
            {
                const auto& frag = m.elements[e].fragments[f];
                std::vector<unsigned> ids;
                for (const Vec2r& b : frag.region_barys) {
                    ids.push_back(sources.size());
                    sources.push_back({e, f, b});
                }
                cells.push_back(ids);
                cell_elem.push_back(e);
            }
    }

    void updatePositions(const Mesh& m, std::vector<Vec2r>& out) const
    {
        out.resize(sources.size());
        for (unsigned i = 0; i < sources.size(); i++) {
            const auto& s = sources[i];
            const auto& frag = m.elements[s.elem].fragments[s.frag];
            Real N[3] = {1 - s.bary[0] - s.bary[1], s.bary[0], s.bary[1]};
            out[i] = N[0]*m.nodes[frag.nodes[0]]
                   + N[1]*m.nodes[frag.nodes[1]]
                   + N[2]*m.nodes[frag.nodes[2]];
        }
    }
};
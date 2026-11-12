#pragma once

#include "common.hpp"
#include "mesh.hpp"

#include <vtkSmartPointer.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkVertex.h>
#include <vtkPolyData.h>
#include <vtkGlyph3DMapper.h>
#include <vtkSphereSource.h>
#include <vtkNamedColors.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkInteractorStyleTrackballCamera.h>

class MeshToggleInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static MeshToggleInteractorStyle* New();

    vtkTypeMacro(MeshToggleInteractorStyle, vtkInteractorStyleTrackballCamera);

    void SetMeshActors(const std::vector<vtkSmartPointer<vtkActor>>& actors) { mesh_actors = actors; }
    void SetNodeActors(const std::vector<vtkSmartPointer<vtkActor>>& actors) { node_actors = actors; }

    void OnKeyPress() override
    {
        std::string key = this->GetInteractor()->GetKeySym();

        // Toggle meshes using number keys.
        if (key.size() == 1 && key[0] >= '0' && key[0] <= '9')
        {
            unsigned index = key[0] - '0';

            if (index < mesh_actors.size())
            {
                auto mesh_actor = mesh_actors[index];
                auto node_actor = node_actors[index];

                mesh_actor->SetVisibility(!mesh_actor->GetVisibility());
                node_actor->SetVisibility(!node_actor->GetVisibility());
                this->GetInteractor()->GetRenderWindow()->Render();
            }

            return;
        }

        // Preserve default keyboard behavior for other keys.
        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

protected:
    MeshToggleInteractorStyle() = default;
    ~MeshToggleInteractorStyle() override = default;

private:
    std::vector<vtkSmartPointer<vtkActor>> mesh_actors;
    std::vector<vtkSmartPointer<vtkActor>> node_actors;
};


//////////////////////////////////////////////////////////////////////////////////////

struct Visualization
{
    std::vector<const Mesh*> meshes;
    std::vector<std::string> mesh_names;

    std::vector<vtkSmartPointer<vtkUnstructuredGrid>> mesh_grids;
    std::vector<vtkSmartPointer<vtkActor>> mesh_actors;
    std::vector<vtkSmartPointer<vtkActor>> node_actors;

    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> render_window;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

    bool show_nodes = true;
    bool show_node_indices = false;
    bool show_element_indices = false;
    bool show_axes = true;

    Visualization(const std::vector<const Mesh*>& meshes, const std::vector<std::string>& mesh_names);

    /** Converts a given mesh into a VTK unstructured grid. */
    void _createVTKGrid(unsigned m_idx);

    /** Sets up the actors for a given mesh */
    void _setupMeshAndNodeActor(unsigned m_idx);

    /** Initializes the rendering pipeline. */
    void _initializeRenderer();

    /** Displays the meshes. */
    void show();
};

#include "visualization.hpp"

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkDataSetMapper.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkLineSource.h>
#include <vtkObjectFactory.h>

Visualization::Visualization(const std::vector<const Mesh*>& meshes_, const std::vector<std::string>& mesh_names_)
    : meshes(meshes_)
    , mesh_names(mesh_names_)
    , mesh_grids(meshes_.size())
    , mesh_actors(meshes_.size())
    , node_actors(meshes_.size())
{
    for (unsigned i = 0; i < meshes.size(); i++)
    {
        mesh_grids[i] = vtkSmartPointer<vtkUnstructuredGrid>::New();
        mesh_actors[i] = vtkSmartPointer<vtkActor>::New();
        node_actors[i] = vtkSmartPointer<vtkActor>::New();

        _createVTKGrid(i);
    }
    

    renderer = vtkSmartPointer<vtkRenderer>::New();

    render_window = vtkSmartPointer<vtkRenderWindow>::New();
    interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    
    _initializeRenderer();
}

void Visualization::_createVTKGrid(unsigned m_idx)
{
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    // Insert mesh nodes.
    const Mesh& mesh = *meshes[m_idx];
    for (const Vec2r& node : mesh.nodes)
    {
        points->InsertNextPoint(
            static_cast<double>(node(0)),
            static_cast<double>(node(1)),
            0.0
        );
    }

    mesh_grids[m_idx]->SetPoints(points);

    // Insert triangular elements.
    for (const TriangleElement& element : mesh.elements)
    {
        vtkIdType ids[3] =
        {
            static_cast<vtkIdType>(element.nodes[0]),
            static_cast<vtkIdType>(element.nodes[1]),
            static_cast<vtkIdType>(element.nodes[2])
        };

        mesh_grids[m_idx]->InsertNextCell(VTK_TRIANGLE, 3, ids);
    }

    mesh_grids[m_idx]->Modified();
}

void Visualization::_setupMeshAndNodeActor(unsigned m_idx)
{
    auto& mesh_actor = mesh_actors[m_idx];
    auto& node_actor = node_actors[m_idx];

    // Mesh mapper.
    vtkSmartPointer<vtkDataSetMapper> mesh_mapper = vtkSmartPointer<vtkDataSetMapper>::New();

    mesh_mapper->SetInputData(mesh_grids[m_idx]);

    mesh_actor->SetMapper(mesh_mapper);

    mesh_actor->GetProperty()->SetColor(0.75, 0.80, 0.90);
    mesh_actor->GetProperty()->SetEdgeColor(0.1, 0.1, 0.1);
    mesh_actor->GetProperty()->SetLineWidth(1.5);
    mesh_actor->GetProperty()->EdgeVisibilityOn();

    // Node visualization.
    vtkSmartPointer<vtkPolyData> node_data =
        vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> node_points =
        vtkSmartPointer<vtkPoints>::New();

    for (const Vec2r& node : meshes[m_idx]->nodes)
    {
        node_points->InsertNextPoint(
            static_cast<double>(node(0)),
            static_cast<double>(node(1)),
            0.01
        );
    }

    node_data->SetPoints(node_points);

    vtkSmartPointer<vtkVertexGlyphFilter> vertex_filter =
        vtkSmartPointer<vtkVertexGlyphFilter>::New();

    vertex_filter->SetInputData(node_data);
    vertex_filter->Update();

    vtkSmartPointer<vtkPolyDataMapper> node_mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();

    node_mapper->SetInputConnection(vertex_filter->GetOutputPort());

    node_actor->SetMapper(node_mapper);

    node_actor->GetProperty()->SetColor(0.85, 0.15, 0.15);
    node_actor->GetProperty()->SetPointSize(8.0);
}

void Visualization::_initializeRenderer()
{
    for (unsigned i = 0; i < meshes.size(); i++)
    {
        _setupMeshAndNodeActor(i);
        renderer->AddActor(mesh_actors[i]);

        if (show_nodes)
            renderer->AddActor(node_actors[i]);
    }

    // Renderer.
    renderer->SetBackground(0.12, 0.14, 0.18);

    // X axis (red).
    auto x_axis = vtkSmartPointer<vtkLineSource>::New();
    x_axis->SetPoint1(0.0, 0.0, 0.00);
    x_axis->SetPoint2(0.2, 0.0, 0.00);

    auto x_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    x_mapper->SetInputConnection(x_axis->GetOutputPort());

    auto x_actor = vtkSmartPointer<vtkActor>::New();
    x_actor->SetMapper(x_mapper);
    x_actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    x_actor->GetProperty()->SetLineWidth(3.0);

    // Y axis (green).
    auto y_axis = vtkSmartPointer<vtkLineSource>::New();
    y_axis->SetPoint1(0.0, 0.0, 0.00);
    y_axis->SetPoint2(0.0, 0.2, 0.00);

    auto y_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    y_mapper->SetInputConnection(y_axis->GetOutputPort());

    auto y_actor = vtkSmartPointer<vtkActor>::New();
    y_actor->SetMapper(y_mapper);
    y_actor->GetProperty()->SetColor(0.0, 1.0, 0.0);
    y_actor->GetProperty()->SetLineWidth(3.0);

    renderer->AddActor(x_actor);
    renderer->AddActor(y_actor);

    // Render window.
    render_window->AddRenderer(renderer);
    render_window->SetSize(1000, 800);
    render_window->SetWindowName("Triangle Mesh Visualization");

    auto style = vtkSmartPointer<MeshToggleInteractorStyle>::New();
    style->SetMeshActors(mesh_actors);
    style->SetNodeActors(node_actors);
    interactor->SetInteractorStyle(style);
    interactor->SetRenderWindow(render_window);

    renderer->ResetCamera();
}

void Visualization::show()
{
    render_window->Render();
    interactor->Initialize();
    interactor->Start();
}

vtkStandardNewMacro(MeshToggleInteractorStyle);

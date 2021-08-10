#include "polyscopersvs.hpp"
#include "RSVSclass.hpp"
#include "RSVSintegration.hpp"
#include "parameters.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "voxel.hpp"
#include <array>
#include <vector>

bool polyscopersvs::POLYSCOPE_INITIALISED = false;

namespace polyscopersvs
{
#ifdef HEADLESS
const std::string polyscopeBackend = "openGL_mock";
#else
const std::string polyscopeBackend = "";
#endif
} // namespace polyscopersvs
/// Initialize polyscope, creating graphics contexts and constructing a window.
/// Should be called exactly once.
polyscopersvs::PolyScopeRSVS::PolyScopeRSVS()
{
    if (!(polyscopersvs::POLYSCOPE_INITIALISED))
    {
        polyscopersvs::POLYSCOPE_INITIALISED = true;
        // TODO: support the mocked backend of polyscope
        polyscope::init(polyscopeBackend);
    }
}

/**
 * @brief Show the polyscope window.
 *
 * Pass control flow to polyscope, displaying the interactive window.
 *  Function will return when user closes the window.
 */
void polyscopersvs::PolyScopeRSVS::show(size_t forFrames)
{
#ifndef HEADLESS
    polyscope::show(forFrames);
#else
    std::cout << "Headless mode, no polyscope window shown." << std::endl;
#endif
}

/**
 * @brief Add a new surface mesh from a mesh object to polyscope.
 *
 * @param name The name of the surface mesh in polyscope. Using an existing name will allow to
 *  overwrite a previous mesh.
 * @param meshIn The mesh object that will be parsed for polyscope. This mesh must have been
 * prepared (Using the `mesh::PrepareForUse` method).
 */
void polyscopersvs::PolyScopeRSVS::addMesh(std::string name, const mesh &meshIn)
{
    if (!meshIn.isready())
    {
        RSVS3D_ERROR("Mesh is not ready, cannot add to polyscope. Call `PrepareForUse` first.");
    }
    std::vector<std::vector<double>> points;
    std::vector<std::vector<int>> faces;
    points.reserve(meshIn.verts.size());
    faces.reserve(meshIn.surfs.size());

    for (size_t i = 0; i < meshIn.verts.size(); ++i)
    {
        if (meshIn.verts(i)->coord.size() != 3)
        {
            RSVS3D_ERROR_RANGE(
                (std::string("Vertex at position ") + std::to_string(i) + " has invalid coordinate.").c_str());
        }

        points.push_back(meshIn.verts(i)->coord);
    }
    for (size_t i = 0; i < meshIn.surfs.size(); ++i)
    {
        auto faceIndices = meshIn.verts.find_list(meshIn.surfs(i)->OrderedVerts(&meshIn));
        if (faceIndices.size() < 3)
        {
            RSVS3D_ERROR_RANGE(
                (std::string("Surface at position") + std::to_string(i) + " has less than 3 vertices.").c_str());
        }
        faces.push_back(faceIndices);
    }

    polyscope::registerSurfaceMesh(name, points, faces);
}

/**
 * @brief
 *
 * @param name
 * @param meshIn
 * @param cellIndices
 * @param isIndex
 * @return float : This is the volume of the last cell that was encountered
 */
float polyscopersvs::PolyScopeRSVS::addCells(std::string name, const mesh &meshIn, const std::vector<int> &&cellIndices,
                                             bool isIndex)
{
    if (!meshIn.isready())
    {
        RSVS3D_ERROR("Mesh is not ready, cannot add to polyscope. Call `PrepareForUse` first.");
    }
    std::vector<std::vector<double>> points;
    std::vector<std::vector<int>> faces;
    points.reserve(meshIn.verts.size());
    faces.reserve(meshIn.surfs.size());
    float outputVolume = 0.0;

    for (size_t i = 0; i < meshIn.verts.size(); ++i)
    {
        if (meshIn.verts(i)->coord.size() != 3)
        {
            RSVS3D_ERROR_RANGE(
                (std::string("Vertex at position ") + std::to_string(i) + " has invalid coordinate.").c_str());
        }

        points.push_back(meshIn.verts(i)->coord);
    }
    std::vector<int> volumePositions;
    if (isIndex)
    {
        volumePositions = meshIn.volus.find_list(cellIndices);
    }
    else
    {
        volumePositions = cellIndices;
    }
    for (auto volumePosition : volumePositions)
    {
        if (volumePosition < 0 || volumePosition >= meshIn.volus.size())
        {
            polyscope::warning("Failed to find cell", "ID : " + std::to_string(volumePosition));
            continue;
        }
        std::cout << "Display cell " << meshIn.volus(volumePosition)->index << " at position [" << volumePosition << "]"
                  << std::endl;
        outputVolume = meshIn.volus(volumePosition)->target;
        for (auto surfaceIndex : meshIn.volus(volumePosition)->surfind)
        {
            auto faceIndices = meshIn.verts.find_list(meshIn.surfs.isearch(surfaceIndex)->OrderedVerts(&meshIn));
            if (faceIndices.size() < 3)
            {
                RSVS3D_ERROR_RANGE(
                    (std::string("Surface ") + std::to_string(surfaceIndex) + " has less than 3 vertices.").c_str());
            }
            faces.push_back(faceIndices);
        }
    }

    polyscope::registerSurfaceMesh(name, points, faces);
    return outputVolume;
}

/**
 * @brief Contains overloads of functions for displaying and configuring the parameter structure.
 *
 * This namespace contains overloads of function `parameterConfigGui` which is used to display a
 *  UI for the configuration of the RSVS parameters using ImGui.
 *
 */
namespace parameter_ui
{
void parameterConfigGui(param::voronoi &paramConf)
{
    paramConf.inputpoints;
    ImGui::InputDouble("Snake coarseness", &paramConf.snakecoarseness, 0.0, 0.0, "%.2f");
    ImGui::InputInt("Voronoi layers", &paramConf.vorosnakelayers);
    ImGui::InputDouble("Bounding box distance", &paramConf.distancebox, 0.0, 0.0, "%.2f");
    paramConf.pointfile.reserve(2048);
    ImGui::InputText("Voronoi points file", paramConf.pointfile.data(), paramConf.pointfile.capacity());
}
void parameterConfigGui(param::voxel &paramConf)
{
    ImGui::InputInt3("Background grid", paramConf.gridsizebackground.data());
    ImGui::InputInt3("Snaking grid", paramConf.gridsizesnake.data());
}

void parameterConfigGui(std::string name, std::array<param::realbounds, 3> &paramConf)
{
    if (ImGui::TreeNode(name.c_str()))
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int i = 0; i < 3; ++i)
            {
                auto fullName = ("##" + name + std::to_string(i) + "," + std::to_string(j)).c_str();
                ImGui::InputDouble(fullName, &(paramConf[i][j]), 0.2, 2.0, "%.2f");
                ImGui::SameLine();
            }
            ImGui::NewLine();
        }
        ImGui::TreePop();
    }
}
void parameterConfigGui(param::grid &paramConf)
{

    if (ImGui::BeginCombo("##Grid Style", paramConf.activegrid.c_str()))
    {
        for (auto key : {"voxel", "voronoi"})
        {

            if (ImGui::Selectable(key, paramConf.activegrid == key))
            {
                paramConf.activegrid = key;
            }
        }
        ImGui::EndCombo();
    }
    parameterConfigGui("Compute domain", paramConf.domain);
    parameterConfigGui("Physical domain", paramConf.physdomain);
    ImGui::SetNextTreeNodeOpen(false, ImGuiCond_FirstUseEver);
    if (ImGui::TreeNode("Voxel"))
    {
        parameterConfigGui(paramConf.voxel);
        ImGui::TreePop();
    }
    ImGui::SetNextTreeNodeOpen(false, ImGuiCond_FirstUseEver);
    if (ImGui::TreeNode("Voronoi"))
    {
        parameterConfigGui(paramConf.voronoi);
        ImGui::TreePop();
    }
}
void parameterConfigGui(param::rsvs &paramConf)
{
    bool change = ImGui::Checkbox("Uniform", &paramConf.cstfill.active);
    if (change && paramConf.cstfill.active)
    {
        paramConf.filefill.active = false;
        paramConf.makefill.active = false;
    }
    ImGui::SameLine();
    ImGui::InputDouble("##Uniform value", &paramConf.cstfill.fill, 0.0, 0.0, "%.2f");

    change = ImGui::Checkbox("File", &paramConf.filefill.active);
    if (change && paramConf.cstfill.active)
    {
        paramConf.makefill.active = false;
        paramConf.cstfill.active = false;
    }
    ImGui::SameLine();
    paramConf.filefill.fill.reserve(2048);
    ImGui::InputText("##File value", paramConf.filefill.fill.data(), paramConf.filefill.fill.capacity());

    change = ImGui::Checkbox("Auto", &paramConf.makefill.active);
    if (change && paramConf.makefill.active)
    {
        paramConf.filefill.active = false;
        paramConf.cstfill.active = false;
    }
    ImGui::SameLine();
    paramConf.makefill.fill.reserve(2048);
    ImGui::InputText("##Auto value", paramConf.makefill.fill.data(), paramConf.makefill.fill.capacity());

    ImGui::InputInt("Linear solver (0-4)", &paramConf.solveralgorithm);
}
void parameterConfigGui(param::snaking &paramConf)
{
}
void parameterConfigGui(param::files &paramConf)
{
}
void parameterConfigGui(param::dev::devparam &paramConf)
{
}
void parameterConfigGui(param::dev::rsvseps &paramConf)
{
}
void parameterConfigGui(param::parameters &paramConf)
{
    paramConf.dev;
    paramConf.files;
    paramConf.grid;
    paramConf.rsvs;
    paramConf.snak;
    if (ImGui::CollapsingHeader("Parameters"))
    {
        ImGui::SetNextTreeNodeOpen(false, ImGuiCond_FirstUseEver);
        if (ImGui::TreeNode("Grid"))
        {
            parameterConfigGui(paramConf.grid);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("RSVS"))
        {
            parameterConfigGui(paramConf.rsvs);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("Snaking"))
        {
            parameterConfigGui(paramConf.snak);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("Files"))
        {
            parameterConfigGui(paramConf.files);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("Development"))
        {
            parameterConfigGui(paramConf.dev);
            ImGui::TreePop();
        }
    }
}
} // namespace parameter_ui

/**
 * @brief A method which sets the polyscope userCallback to a UI controlling the giving
 *   RSVS object.
 *
 * This will make the current object active in the session.
 *
 * @param RSVSobj A fully formed RSVSclass which allows handling of the RSVS iteration process.
 *
 */
void polyscopersvs::PolyScopeRSVS::setInteractiveCallback(integrate::RSVSclass &RSVSobj)
{
    integrate::iteratereturns iterateInfo(0, 0, 0.0);
    int inspectId = 1;
    bool volumeByIndex = true;
    float newVolumeValue = 0.0;
    auto callback = [&] {
        ImGui::PushItemWidth(100);

        if (ImGui::Button("View surfaces"))
        {
            // executes when button is pressed
            RSVSobj.rsvsSnake.PrepareForUse();
            RSVSobj.voluMesh.PrepareForUse();
            integrate::execute::logging::SnakePolyscope(RSVSobj.viewer, RSVSobj.rsvsSnake);
            this->addMesh("Snaking mesh", *RSVSobj.rsvsSnake.snakemesh());
            this->addMesh("Volume mesh", RSVSobj.voluMesh);

            polyscope::getSurfaceMesh("Snaking mesh")->setEnabled(false);
            polyscope::getSurfaceMesh("Volume mesh")->setTransparency(0.5);
            auto snakeSurfaceMesh = polyscope::getSurfaceMesh(integrate::constants::polyscopeSnakeName);
            snakeSurfaceMesh->setEdgeColor({0, 0, 0});
            snakeSurfaceMesh->setEdgeWidth(1);
            snakeSurfaceMesh->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
        }

        if (ImGui::CollapsingHeader("Configuration"))
        {

            if (ImGui::InputInt("Cell ID", &inspectId))
            {
                RSVSobj.voluMesh.PrepareForUse();
                newVolumeValue = this->addCells("Active Cell", RSVSobj.voluMesh, {inspectId}, volumeByIndex);
            }
            ImGui::SameLine();
            ImGui::Checkbox("Use ID (or position)", &volumeByIndex);

            if (ImGui::InputFloat("New volume", &newVolumeValue))
            {
                // executes when button is pressed
                RSVSobj.voluMesh.PrepareForUse();
                int volumePosition = inspectId;
                if (volumeByIndex)
                {
                    volumePosition = RSVSobj.voluMesh.volus.find(inspectId);
                }
                if (volumePosition > -1 && volumePosition < RSVSobj.voluMesh.volus.size())
                {

                    RSVSobj.voluMesh.volus[volumePosition].target =
                        newVolumeValue > 0.0 ? (newVolumeValue <= 1.0 ? newVolumeValue : 1.0) : 0.0;
                    RSVSobj.voluMesh.PrepareForUse();
                }
                else
                {
                    polyscope::warning("Invalid Cell Position", "The cell was not found, it is either out of bounds"
                                                                " or not a valid index.");
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Set All "))
            {
                RSVSobj.voluMesh.PrepareForUse();
                newVolumeValue = newVolumeValue > 0.0 ? (newVolumeValue <= 1.0 ? newVolumeValue : 1.0) : 0.0;
                for (size_t i = 0; i < RSVSobj.voluMesh.volus.size(); ++i)
                {
                    RSVSobj.voluMesh.volus[i].target = newVolumeValue;
                }
            }
        }
        parameter_ui::parameterConfigGui(RSVSobj.paramconf);
        if (ImGui::CollapsingHeader("Execution"))
        {
            if (ImGui::Button("Run RSVS"))
            {
                // executes when button is pressed
                iterateInfo = integrate::execute::RSVSiterate(RSVSobj);
            }
            ImGui::SameLine();
            ImGui::InputInt("Number of steps", &RSVSobj.paramconf.snak.maxsteps);

            if (ImGui::Button("Postprocess RSVS"))
            {
                // executes when button is pressed
                integrate::execute::PostProcessing(RSVSobj, iterateInfo.timeT, iterateInfo.nVoluZone,
                                                   iterateInfo.stepNum);
                integrate::execute::Exporting(RSVSobj);
            }
            ImGui::SameLine();
            if (ImGui::Button("Reset"))
            {
                integrate::Prepare(RSVSobj);
            }
        }
        ImGui::PopItemWidth();
    };
    polyscope::state::userCallback = callback;
}

/**
 * @brief Should initialise polyscope
 *
 * @return int The error code: 0 indicates succes, all other numbers are failures
 */

int polyscopersvs::test::init()
{
    polyscopersvs::PolyScopeRSVS viewer;
    assert(polyscopersvs::POLYSCOPE_INITIALISED);
    return 0;
}

/**
 * @brief Should display the polyscope window
 *
 * @return int The error code: 0 indicates succes, all other numbers are failures
 */
int polyscopersvs::test::show()
{
    polyscopersvs::PolyScopeRSVS viewer;
    viewer.show(30);
    return 0;
}

/**
 * @brief Should display a single small cubic mesh
 *
 * @return int The error code: 0 indicates succes, all other numbers are failures
 */
int polyscopersvs::test::meshShow()
{
    int errorOut = 0;
    polyscopersvs::PolyScopeRSVS viewer;
    mesh mesh1;

    std::array<int, 3> dims = {2, 3, 4};
    errorOut += BuildBlockGrid(dims, mesh1);
    mesh1.PrepareForUse();
    viewer.addMesh("mesh1", mesh1);
    viewer.show(60);
    return 0;
}

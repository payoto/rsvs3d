/**
 * @file polyscopersvs.cpp
 * @author Alexandre Payot
 * @brief Provide the integration of the RSVS with the GUI library polyscope
 * @version 0.1
 * @date 2021-08-10
 *
 * This file provides the interface between the RSVS and polyscope. This interface
 * relies on the class `polyscopersvs::PolyScopeRSVS`.
 *
 * @copyright Copyright (c) 2021
 */
#include "polyscopersvs.hpp"
#include "RSVSclass.hpp"
#include "RSVSintegration.hpp"
#include "meshprocessing.hpp"
#include "parameters.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "rsvsjson.hpp"
#include "voxel.hpp"
#include <array>
#include <sstream>
#include <vector>

bool polyscopersvs::POLYSCOPE_INITIALISED = false;
bool polyscopersvs::test::TEST_HEADLESS = false;
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
#ifndef HEADLESS
    this->isHeadless = true;
#else
    this->isHeadless = false;
#endif
    if (!(polyscopersvs::POLYSCOPE_INITIALISED))
    {
        polyscopersvs::POLYSCOPE_INITIALISED = true;
        // TODO: support the mocked backend of polyscope
        polyscope::init(polyscopersvs::polyscopeBackend);
    }
}

polyscopersvs::PolyScopeRSVS::PolyScopeRSVS(bool isHeadless)
{
    this->isHeadless = isHeadless;
    if (!(polyscopersvs::POLYSCOPE_INITIALISED))
    {
        polyscopersvs::POLYSCOPE_INITIALISED = true;
        // TODO: support the mocked backend of polyscope
        if (isHeadless)
        {
            polyscope::init("openGL_mock");
        }
        else
        {
            polyscope::init();
        }
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
    if (!(this->isHeadless))
    {
        polyscope::show(forFrames);
    }
    else
    {
        std::cout << "Headless mode, no polyscope window shown." << std::endl;
    }
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
    if (this->isHeadless)
        return;

    if (!meshIn.isready())
    {
        RSVS3D_ERROR("Mesh is not ready, cannot add to polyscope. Call `PrepareForUse` first.");
    }
    std::vector<std::vector<double>> points;
    std::vector<std::vector<int>> faces;
    points.reserve(meshIn.verts.size());
    faces.reserve(meshIn.surfs.size());

    for (size_t i = 0; i < size_t(meshIn.verts.size()); ++i)
    {
        if (meshIn.verts(i)->coord.size() != 3)
        {
            RSVS3D_ERROR_RANGE(
                (std::string("Vertex at position ") + std::to_string(i) + " has invalid coordinate.").c_str());
        }

        points.push_back(meshIn.verts(i)->coord);
    }
    for (size_t i = 0; i < size_t(meshIn.surfs.size()); ++i)
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
 * @brief Color a surface mesh using a snake's properties
 *
 * Uses the snake's velocity to color a new surface mesh.
 */
void polyscopersvs::PolyScopeRSVS::addSnake(std::string name, const snake &snakeIn)
{
    if (this->isHeadless)
        return;

    this->addMesh(name, snakeIn.snakeconn);
    std::vector<double> velocities;
    velocities.reserve(snakeIn.snaxs.size());
    for (size_t i = 0; i < size_t(snakeIn.snaxs.size()); ++i)
    {

        velocities.push_back(snakeIn.snaxs(i)->v);
    }
    polyscope::getSurfaceMesh(name)->addVertexScalarQuantity("velocity", velocities);
}

/**
 * @brief Wrap a vector to allow for (i,j) indexing as if it was a 2D
 * array of coordinates.
 *
 * Stride is hard coded at 3. No checks are run.
 *
 */
class _3dVectorWrapper
{
  public:
    _3dVectorWrapper(const std::vector<double> &vec, double multiplier = 1.0) : vec(vec), multiplier(multiplier)
    {
    }
    double operator()(size_t i, size_t j) const
    {
        return vec[i * 3 + j] * this->multiplier;
    }
    size_t size() const
    {
        return vec.size() / 3;
    }

  private:
    const std::vector<double> &vec;
    double multiplier = 1.0;
};

/**
 * @brief Plots all the available surface quantities on a surface mesh.
 *
 * @param name
 * @param surfaceMesh
 */
void polyscopersvs::PolyScopeRSVS::addSurfaceProperties(std::string name, const mesh &surfaceMesh)
{
    if (this->isHeadless)
        return;

    std::vector<double> properties;

    auto psMesh = polyscope::getSurfaceMesh(name);

    psMesh->addEdgeScalarQuantity("edge-curvature", CalculateEdgeCurvature(surfaceMesh));
    psMesh->addEdgeScalarQuantity("edge-length", CalculateEdgeLengths(surfaceMesh));
    psMesh->addVertexScalarQuantity("vertex-curvature", CalculateVertexCurvature(surfaceMesh, 1));
    psMesh->addVertexScalarQuantity("vertex-minimum-length", CalculateVertexMinEdgeLength(surfaceMesh));
    psMesh->addVertexScalarQuantity("vertex-maximum-length", CalculateVertexMaxEdgeLength(surfaceMesh));
    psMesh->addVertexScalarQuantity("vertex-mean-length", CalculateVertexMeanEdgeLength(surfaceMesh));
    psMesh->addVertexVectorQuantity("vertex-normal", _3dVectorWrapper(MeshUnitNormals(surfaceMesh)));
    psMesh->addVertexVectorQuantity("vertex-laplacians-(flipped)", _3dVectorWrapper(MeshLaplacians(surfaceMesh), -1.0));
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

    // Todo: only access vertices that are in the cell.
    for (size_t i = 0; i < size_t(meshIn.verts.size()); ++i)
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
        if (this->isHeadless)
            continue;
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
    // Static variable keeps track of which entry is active in the voronoi config.
    static int activeEntry;

    ImGui::InputDouble("Snake coarseness", &paramConf.snakecoarseness, 0.0, 0.0, "%.2f");
    ImGui::InputInt("Voronoi layers", &paramConf.vorosnakelayers);
    ImGui::InputDouble("Bounding box distance", &paramConf.distancebox, 0.0, 0.0, "%.2f");
    paramConf.pointfile.reserve(2048);
    ImGui::InputText("Voronoi points file", paramConf.pointfile.data(), paramConf.pointfile.capacity());

    if (ImGui::Button("Add Voronoi point"))
    {
        for (size_t i = 0; i < 4; i++)
        {
            paramConf.inputpoints.push_back(0.5);
            // Hacky condition to ensure we finish when we hit a multiple of 4.
            if (paramConf.inputpoints.size() % 4 == 0)
            {
                break;
            }
        }
    }
    int nEntries = paramConf.inputpoints.size() / 4;
    // This is used to initalise the static variable if it is out of bounds
    if (!(activeEntry < nEntries && activeEntry >= 0))
    {
        activeEntry = 0;
    }

    if (nEntries > 0)
    {
        ImGui::InputInt(("Active Point, max ID: " + std::to_string(nEntries)).c_str(), &activeEntry);
        if (!(activeEntry < nEntries && activeEntry >= 0))
        {
            polyscope::warning("ID to access input point was incorrect.", "Maximum value is " +
                                                                              std::to_string(nEntries) + " value was " +
                                                                              std::to_string(activeEntry));
            activeEntry = 0;
        }
        ImGui::Text("Voronoi point definition - Active ID: %i\n   X   ,   Y   ,   Z   , Volume fraction", activeEntry);
        for (int i = 0; i < 4; ++i)
        {
            auto fullName = ("## voropoint setter " + std::to_string(i)).c_str();
            ImGui::InputDouble(fullName, &(paramConf.inputpoints[activeEntry * 4 + i]), 0.0, 0.0, "%.2f");
            ImGui::SameLine();
        }
        ImGui::NewLine();
    }
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
    ImGui::InputDouble("Arrival tolerance", &paramConf.arrivaltolerance, 0.0, 0.0, "%.2e");
    ImGui::InputDouble("Simulataneous arrival tolerance", &paramConf.multiarrivaltolerance, 0.0, 0.0, "%.2e");
    ImGui::InputDouble("Max step distance", &paramConf.snaxdiststep, 0.0, 0.0, "%.2e");
    ImGui::InputDouble("Max step time", &paramConf.snaxtimestep, 0.0, 0.0, "%.2e");
    ImGui::InputDouble("Spawn position", &paramConf.spawnposition, 0.0, 0.0, "%.2e");
    ImGui::InputInt("Initialisation boundary (0/1)", &paramConf.initboundary);
}
void parameterConfigGui(param::files &paramConf)
{
    ImGui::InputInt("Logging level [0-7]", &paramConf.ioout.logginglvl);
    ImGui::InputText("Case name", paramConf.ioin.casename.data(), paramConf.ioin.casename.capacity());
    rsvsjson::json out;
    out = paramConf;
    std::ostringstream stream;
    stream << out.dump(2);
    ImGui::TextUnformatted(stream.str().c_str());
}

void parameterConfigGui(param::dev::devparam &paramConf)
{
    ImGui::InputDouble("limitlagrangian", &paramConf.limitlagrangian, 0.0, 0.0, "%.2e");
    ImGui::InputInt("mindesvarsparse", &paramConf.mindesvarsparse);
    ImGui::InputTextWithHint("smoothstepmethod", "Not sure what appropriate values are",
                             paramConf.smoothstepmethod.data(), paramConf.smoothstepmethod.capacity());
    ImGui::Checkbox("snaxDistanceLimit_conserveShape", &paramConf.snaxDistanceLimit_conserveShape);
    ImGui::Checkbox("surfcentrehessian", &paramConf.surfcentrehessian);
    ImGui::Checkbox("surfcentrejacobian", &paramConf.surfcentrejacobian);

    if (ImGui::TreeNode("Smoothing tolerances"))
    {
        ImGui::InputDouble("rsvsmath_automatic_eps_centre2", &paramConf.rsvsepsilons.rsvsmath_automatic_eps_centre2,
                           0.0, 0.0, "%.2e");
        ImGui::InputDouble("rsvsmath_automatic_eps_centre", &paramConf.rsvsepsilons.rsvsmath_automatic_eps_centre, 0.0,
                           0.0, "%.2e");
        ImGui::InputDouble("rsvsmath_automatic_eps_edge", &paramConf.rsvsepsilons.rsvsmath_automatic_eps_edge, 0.0, 0.0,
                           "%.2e");
        ImGui::InputDouble("rsvsmath_automatic_eps_surf", &paramConf.rsvsepsilons.rsvsmath_automatic_eps_surf, 0.0, 0.0,
                           "%.2e");
        ImGui::InputDouble("rsvsmath_automatic_eps_volu", &paramConf.rsvsepsilons.rsvsmath_automatic_eps_volu, 0.0, 0.0,
                           "%.2e");

        ImGui::TreePop();
    }
}
void parameterExportImportGui(param::parameters &paramConf)
{
    int fileNameCapacity = 1024;
    static std::string importFile;
    static std::string exportFile;
    static bool flatJsonExport = false;
    // Make sure the capacity of importFile is large enough
    if (importFile.capacity() < size_t(fileNameCapacity))
    {
        importFile.reserve(fileNameCapacity);
        importFile = "";
        exportFile.reserve(fileNameCapacity);
        exportFile = "interactive-export.json";
    }

    // Buttons to import and export parameters
    ImGui::PushItemWidth(200);
    ImGui::InputTextWithHint("##Import-parameter", "import.json...", importFile.data(), fileNameCapacity);
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::Button("Import"))
    {
        try
        {
            param::io::read(importFile, paramConf);
        }
        catch (const std::exception &e)
        {
            polyscope::error(std::string("Could not import parameters from file '") + importFile + "': " + e.what());
        }
    }
    ImGui::PushItemWidth(200);
    ImGui::InputTextWithHint("##Export-parameter", "export.json...", exportFile.data(), fileNameCapacity);
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::Button("Export"))
    {
        try
        {
            if (!flatJsonExport)
            {
                param::io::write(exportFile, paramConf);
            }
            else
            {
                param::io::writeflat(exportFile, paramConf);
            }
        }
        catch (const std::exception &e)
        {
            polyscope::error(std::string("Could not export parameters to file '") + exportFile + "': " + e.what());
        }
    }
    ImGui::SameLine();
    ImGui::Checkbox("Flat JSON format", &flatJsonExport);
}
void parameterConfigGui(param::parameters &paramConf)
{
    // Buttons for configuring the parameter structure
    if (ImGui::CollapsingHeader("Parameters"))
    {
        parameterExportImportGui(paramConf);

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
        if (ImGui::TreeNode("File setting (view)"))
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

namespace vos_ui
{
void UpdateVolume(integrate::RSVSclass &RSVSobj, int inspectId, bool volumeByIndex, double newVolumeValue)
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
void UpdateVolumes(integrate::RSVSclass &RSVSobj, double newVolumeValue)
{
    RSVSobj.voluMesh.PrepareForUse();
    newVolumeValue = newVolumeValue > 0.0 ? (newVolumeValue <= 1.0 ? newVolumeValue : 1.0) : 0.0;
    for (size_t i = 0; i < size_t(RSVSobj.voluMesh.volus.size()); ++i)
    {
        RSVSobj.voluMesh.volus[i].target = newVolumeValue;
    }
}
void vosExportImportGui(integrate::RSVSclass &RSVSobj)
{
    static int lineLength = 8;
    static int vosPrecision = 4;
    size_t fileNameCapacity = 1024;
    static std::string importFile;
    static std::string exportFile;
    // Make sure the capacity of importFile is large enough
    if (importFile.capacity() < fileNameCapacity)
    {
        importFile.reserve(fileNameCapacity);
        importFile = "";
        exportFile.reserve(fileNameCapacity);
        exportFile = "interactive-export.fill";
    }
    if (ImGui::TreeNode("Export/Import of Volume fractions"))
    {
        // Buttons to import and export parameters
        ImGui::PushItemWidth(200);
        ImGui::InputTextWithHint("##Import-vosfill", "import.fill...", importFile.data(), fileNameCapacity);
        ImGui::PopItemWidth();
        ImGui::SameLine();
        if (ImGui::Button("Import"))
        {
            try
            {
                RSVSobj.voluMesh.LoadTargetFill(importFile);
            }
            catch (const std::exception &e)
            {
                polyscope::error(std::string("Could not import VOS from file '") + importFile + "': " + e.what());
            }
        }
        ImGui::PushItemWidth(200);
        ImGui::InputTextWithHint("##Export-vosfill", "export.fill...", exportFile.data(), fileNameCapacity);
        ImGui::PopItemWidth();
        ImGui::SameLine();
        if (ImGui::Button("Export"))
        {
            try
            {
                // Write out the target fill of each Volume cell to the export file
                std::fstream file(exportFile, std::ios::out);
                file << std::setprecision(vosPrecision);
                for (int i = 0; i < RSVSobj.voluMesh.volus.size(); ++i)
                {
                    file << RSVSobj.voluMesh.volus[i].target << " ";
                    // every lineLength characters, write a newline
                    if (i % lineLength == lineLength - 1)
                    {
                        file << std::endl;
                    }
                }
            }
            catch (const std::exception &e)
            {
                polyscope::error(std::string("Could not export parameters to file '") + exportFile + "': " + e.what());
            }
        }

        // Input for line length of the output file
        ImGui::InputInt("Items per line", &lineLength);
        ImGui::SameLine();
        ImGui::InputInt("Precision", &vosPrecision);
        ImGui::TreePop();
    }
}
void vosControlGui(integrate::RSVSclass &RSVSobj)
{
    static int inspectId = 1;
    static bool volumeByIndex = true;
    static float newVolumeValue = 0.0;
    if (ImGui::CollapsingHeader("VOS Configuration"))
    {
        // A button to import and export the VOS configuration
        vosExportImportGui(RSVSobj);
        ImGui::Text("\nManual control and viewing Volume fractions");
        bool cellIdChange = ImGui::InputInt("Cell ID", &inspectId);
        ImGui::SameLine();
        ImGui::Checkbox("Use ID (or position)", &volumeByIndex);
        bool volumeHasChanged = ImGui::InputFloat("New volume", &newVolumeValue);
        ImGui::SameLine();
        bool volumesHaveChanged = ImGui::Button("Set All");
        if (volumeHasChanged)
        {
            UpdateVolume(RSVSobj, inspectId, volumeByIndex, newVolumeValue);
        }
        if (volumesHaveChanged)
        {
            UpdateVolumes(RSVSobj, newVolumeValue);
        }
        if (cellIdChange || volumeHasChanged || volumesHaveChanged)
        {
            RSVSobj.voluMesh.PrepareForUse();
            newVolumeValue = RSVSobj.viewer.addCells("Active Cell", RSVSobj.voluMesh, {inspectId}, volumeByIndex);
            polyscope::getSurfaceMesh("Active Cell")
                ->setSurfaceColor({1.0 - newVolumeValue, 1.0 - newVolumeValue, 1.0 - newVolumeValue});
        }
    }
}

} // namespace vos_ui

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

    /**
     * @brief A ImGui callback for configuring and controlling the RSVS execution.
     *
     * This lambda function is only meant to be used as the user callback of polyscope.
     */
    auto callback = [&] {
        bool runPreparation = false;
        ImGui::PushItemWidth(100);

        bool viewSurfaces = ImGui::Button("View surfaces");

        if (ImGui::CollapsingHeader("RSVS Execution"))
        {
            ImGui::InputInt("Number of steps", &RSVSobj.paramconf.snak.maxsteps);
            if (ImGui::Button("Run"))
            {
                // executes when button is pressed
                iterateInfo = integrate::execute::RSVSiterate(RSVSobj);
            }
            ImGui::SameLine();

            if (ImGui::Button("Postprocess"))
            {
                // executes when button is pressed
                integrate::execute::PostProcessing(RSVSobj, iterateInfo.timeT, iterateInfo.nVoluZone,
                                                   iterateInfo.stepNum);
                integrate::execute::Exporting(RSVSobj);
            }
            ImGui::SameLine();
            runPreparation = ImGui::Button("Reset");
            if (runPreparation)
            {
                integrate::Prepare(RSVSobj);
            }

            RSVSobj.paramconf.snak.engine.reserve(2048);
            ImGui::InputText("Velocity Engine", RSVSobj.paramconf.snak.engine.data(),
                             RSVSobj.paramconf.snak.engine.capacity());
        }
        if (viewSurfaces || runPreparation)
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
        vos_ui::vosControlGui(RSVSobj);
        parameter_ui::parameterConfigGui(RSVSobj.paramconf);

        ImGui::PopItemWidth();
    };
    polyscope::state::userCallback = callback;
}

/**
 * @brief Should initialise polyscope
 *
 * @return int The error code: 0 indicates success, all other numbers are failures
 */
int polyscopersvs::test::init()
{
    polyscopersvs::PolyScopeRSVS viewer(polyscopersvs::test::TEST_HEADLESS);
    assert(polyscopersvs::POLYSCOPE_INITIALISED);
    return 0;
}

/**
 * @brief Should display the polyscope window
 *
 * @return int The error code: 0 indicates success, all other numbers are failures
 */
int polyscopersvs::test::show()
{
    polyscopersvs::PolyScopeRSVS viewer(polyscopersvs::test::TEST_HEADLESS);
    viewer.show(30);
    return 0;
}

/**
 * @brief Should display a single small cubic mesh
 *
 * @return int The error code: 0 indicates success, all other numbers are failures
 */
int polyscopersvs::test::meshShow()
{
    int errorOut = 0;
    polyscopersvs::PolyScopeRSVS viewer(polyscopersvs::test::TEST_HEADLESS);
    mesh mesh1;

    std::array<int, 3> dims = {2, 3, 4};
    errorOut += BuildBlockGrid(dims, mesh1);
    mesh1.PrepareForUse();
    viewer.addMesh("mesh1", mesh1);
    viewer.show(60);
    return 0;
}

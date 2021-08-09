#include "polyscopersvs.hpp"
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

#include "polyscopersvs.hpp"
#include "polyscope/polyscope.h"

bool polyscopersvs::POLYSCOPE_INITIALISED = false;

/// Initialize polyscope, creating graphics contexts and constructing a window.
/// Should be called exactly once.
polyscopersvs::PolyScopeRSVS::PolyScopeRSVS()
{
    if (!(polyscopersvs::POLYSCOPE_INITIALISED))
    {
        polyscopersvs::POLYSCOPE_INITIALISED = true;
        // TODO: support the mocked backend of polyscope
        polyscope::init();
    }
}

/**
 * @brief Show the polyscope window.
 *
 * Pass control flow to polyscope, displaying the interactive window.
 *  Function will return when user closes the window.
 */
void polyscopersvs::PolyScopeRSVS::show()
{
    polyscope::show();
}

int polyscopersvs::test::init()
{
    polyscopersvs::PolyScopeRSVS viewer;
    return 0;
}

int polyscopersvs::test::show()
{
    polyscopersvs::PolyScopeRSVS viewer;
    viewer.show();
    return 0;
}

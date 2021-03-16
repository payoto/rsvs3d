#include "rsvsjson.hpp"

#include "warning.hpp"
void rsvsjson::flatupdate(rsvsjson::json &jfin, rsvsjson::json &jnew, bool isFlatFin, bool isFlatNew)
{
    /*
    Allows recursing update into sub fields
    */
    if (!isFlatNew)
    {
        jnew = jnew.flatten();
    }
    if (!isFlatFin)
    {
        jfin = jfin.flatten();
    }

    // Insert values read into the parameter structure
    jfin.update(jnew);

    try
    {
        jnew = jnew.unflatten();
    }
    catch (std::exception const &ex)
    {
        std::cerr << jnew.dump(2) << std::endl;
        std::cerr << ex.what() << std::endl;
        RSVS3D_ERROR("Could not unflatten jnew");
    }
    try
    {
        jfin = jfin.unflatten();
    }
    catch (std::exception const &ex)
    {
        std::cerr << jfin.dump(2) << std::endl;
        std::cerr << ex.what() << std::endl;
        RSVS3D_ERROR("Could not unflatten jfin");
    }
}

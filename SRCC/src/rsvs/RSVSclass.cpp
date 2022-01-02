
#include "RSVSclass.hpp"

#include <memory>

#include "RSVScalc.hpp"
#include "SnakeVelocities.hpp"
#include "triangulate.hpp"

void integrate::RSVSclass::setVelocityEngine(const std::string &engine)
{
    if (engine.compare("rsvs") == 0)
        this->calcObj = std::make_shared<RSVScalc>();
    else if (engine.compare("random-slow") == 0)
        this->calcObj = std::make_shared<rsvs3d::VelocityFunction>(CalculateSnakeVel);
    else if (engine.compare("random") == 0)
        this->calcObj = std::make_shared<rsvs3d::VelocityFunction>(CalculateSnakeVelRand);
    else if (engine.compare("unit") == 0)
        this->calcObj = std::make_shared<rsvs3d::VelocityFunction>(CalculateSnakeVelUnit);
    else if (engine.compare("unit-reflect") == 0)
        this->calcObj = std::make_shared<rsvs3d::VelocityFunction>(CalculateSnakeVelUnitReflect);
    else if (engine.compare("fast") == 0)
        this->calcObj = std::make_shared<rsvs3d::VelocityFunction>(CalculateSnakeVelFast);
    else
    {
        std::cerr << "Engine " << engine << std::endl;
        RSVS3D_ERROR_ARGUMENT("Unknown engine option passed to the RSVS");
    }
}

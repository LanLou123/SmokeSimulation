

#include "constants.h"

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {50, 30, 50};
#endif
const double dt =0.085;
const double theCellSize = 0.5;
const double fluidDensity = 2;
const double theAirDensity = 1;
const double theBoundConstant = (fluidDensity * theCellSize) / dt;
const double theBuoyancyAlpha = 0.08; // Gravity's effect on the smoke particles.
const double theBuoyancyBeta = 0.97; // Buoyancy's effect due to temperature difference.
const double theBuoyancyAmbientTemperature = 0.0; // Ambient temperature.
const double theVorticityEpsilon = 0.1;


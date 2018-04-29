

#include "constants.h"

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {40, 30, 40};
#endif
const double dt =0.008;
const double theCellSize = 0.5;
const double fluidDensity = 1;
const double theAirDensity = 1;
const double theBoundConstant = (fluidDensity * theCellSize) / dt;
const double theBuoyancyAlpha = 0.08; // Gravity's effect on the smoke particles.
const double theBuoyancyBeta = 5.37; // Buoyancy's effect due to temperature difference.
const double theBuoyancyAmbientTemperature = 0.0; // Ambient temperature.
const double theVorticityEpsilon = 0.2;




#include "constants.h"

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {20, 20, 20};
#endif
const double dt =0.004;
const double theCellSize = 0.5;
double fluidDensity = 1;
const double theAirDensity = 0.01;
const double theBoundConstant = (fluidDensity * theCellSize) / dt;
const double theBuoyancyAlpha = 0.08; // Gravity's effect on the smoke particles.
const double theBuoyancyBeta = 0.037; // Buoyancy's effect due to temperature difference.
const double theBuoyancyAmbientTemperature = 0.0; // Ambient temperature.
const double theVorticityEpsilon = 4;


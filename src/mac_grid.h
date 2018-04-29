#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)
# define OMParallelize

# ifdef OMParallelize
# define TOTALThreads 8
# endif
#include "open_gl_headers.h" 
#include "vec.h"
#include "grid_data.h"
#include "grid_data_matrix.h" 
#include <Partio.h>
class Camera;

class MACGrid
{

public:
	MACGrid();
	~MACGrid();
	MACGrid(const MACGrid& orig);
	MACGrid& operator=(const MACGrid& orig);

	void reset();

	void draw(const Camera& c);
	void updateSources();
	void advectVelocity(double dt);
	void addExternalForces(double dt);
	void project(double dt);
	void advectTemperature(double dt);
	void advectDensity(double dt);
	void advectRenderingParticles(double dt);

protected:

	// Setup
	void initialize();

	// Simulation
	void computeBouyancy(double dt);
	void computeVorticityConfinement(double dt);
    void computeCentralVel();
    void computeOmega();
    void computeOmegaGradient();
    void computeDivergence();
    void computeBound();
	// Rendering
	struct Cube { vec3 pos; vec4 color; double dist; };
	void drawWireGrid();
	void drawSmokeCubes(const Camera& c);
	void drawSmoke(const Camera& c);
	void drawCube(const MACGrid::Cube& c);
	void drawFace(const MACGrid::Cube& c);
	void drawVelocities();
	vec4 getRenderColor(int i, int j, int k);
	vec4 getRenderColor(const vec3& pt);
	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);

	// GridData accessors
	enum Direction { X, Y, Z };
	vec3 getVelocity(const vec3& pt);
	double getVelocityX(const vec3& pt);
	double getVelocityY(const vec3& pt);
	double getVelocityZ(const vec3& pt);
	double getTemperature(const vec3& pt);
	double getDensity(const vec3& pt);
	vec3 getCenter(int i, int j, int k);

	void checkPressure( int& i, int& j, int& k, const GridData& p, vec3& pLow, vec3& pHigh );
	vec3 getRewoundPosition(const vec3 & currentPosition, const double dt);
	vec3 clipToGrid(const vec3& outsidePoint, const vec3& insidePoint);
	double getSize(int dimension);
	int getCellIndex(int i, int j, int k);
    void getCellIndexReverse(int idx,int & i,int & j,int & k );
	int getNumberOfCells();
	bool isValidCell(int i, int j, int k);
	bool isValidFace(int dimension, int i, int j, int k);
	vec3 getFacePosition(int dimension, int i, int j, int k);
	void calculateAMatrix();
	bool preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
    bool conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	void calculatePreconditioner(const GridDataMatrix & A);
	void applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z);
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);


	GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
	GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
	GridDataZ mW; // W component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
	GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
	GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
	GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ
    GridData centeral_vel_x;
    GridData centeral_vel_y;
    GridData centeral_vel_z;

    GridData omegaX;
    GridData omegaY;
    GridData omegaZ;
    GridData omegaN;

    GridData omegaGX;
    GridData omegaGY;
    GridData omegaGZ;

    GridData vorConfFX;
    GridData vorConfFY;
    GridData vorConfFZ;

    GridData diverGence;

	GridDataMatrix AMatrix;
	GridData precon;
    vec3 sphereCentr;

public:

	// rendering particles
	std::vector<vec3> rendering_particles;
	std::vector<vec3> rendering_particles_vel;

	enum RenderMode { CUBES, SHEETS };
	static RenderMode theRenderMode;
	static bool theDisplayVel;
	
	
	void saveSmoke(const char* fileName);
	void saveParticle(std::string filename);
	void saveDensity(std::string filename);
};

#endif
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "smoke_sim.h"
#include "constants.h" 
#include "open_gl_headers.h" 
#include "stb_image_write.h" 
#include "custom_output.h" 
#include "basic_math.h"
#include <fstream>

SmokeSim::SmokeSim() : mFrameNum(0), mTotalFrameNum(0), mRecordEnabled(true)
{
   reset();
}

SmokeSim::~SmokeSim()
{
}

void SmokeSim::reset()
{
   mGrid.reset();
	mTotalFrameNum = 0;
}

/*
void SmokeSim::setGridDimensions(int x, int y, int z)
{
   //extern int theDim[3]; // Naughty globals...
   theDim[0] = x;
   theDim[1] = y;
   theDim[2] = z;
   reset();
}
*/

void SmokeSim::step()
{

    std::cout<<"============================================================="<<mTotalFrameNum<<std::endl;
    std::cout<<"============================================================="<<mTotalFrameNum<<std::endl;
    std::cout<<"============================================================="<<mTotalFrameNum<<std::endl;
    std::cout<<"============================================================="<<mTotalFrameNum<<std::endl;
    std::cout<<"============================================================="<<mTotalFrameNum<<std::endl;
    std::cout<<"============================================================="<<mTotalFrameNum<<std::endl;
    std::cout<<"============================================================="<<mTotalFrameNum<<std::endl;
    std::cout<<mGrid.sphereC;
	double dt = 0.085;//0.1;
	mGrid.sphereC+=vec3(0,0,0.06);
   // Step0: Gather user forces
    if(mTotalFrameNum<=130)
   mGrid.updateSources();

   // Step1: Calculate new velocities
   mGrid.advectVelocity(dt);
   mGrid.addExternalForces(dt);
   mGrid.project(dt);

   // Step2: Calculate new temperature
   mGrid.advectTemperature(dt);

   // Step3: Calculate new density 
   mGrid.advectDensity(dt);


	// Step4: Advect rendering particles
	mGrid.advectRenderingParticles(dt);


	mTotalFrameNum++;
}

void SmokeSim::setRecording(bool on, int width, int height)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
	
	
	recordWidth = width;
	recordHeight = height;
}

bool SmokeSim::isRecording()
{
   return mRecordEnabled;
}

void SmokeSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
   if (mRecordEnabled) grabScreen();
}

void SmokeSim::drawAxes()
{
	glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
		glDisable(GL_LIGHTING);

		glLineWidth(2.0); 
		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(1.0, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 1.0, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 1.0);
		glEnd();
	glPopAttrib();
}

void SmokeSim::grabScreen()  // Code adapted from asst#1 . USING STB_IMAGE_WRITE INSTEAD OF DEVIL.
{
	if (mFrameNum > 9999) exit(0);

	// Save density field to a .bgeo file
	std::string densityFile = "../records/DensityFrame" + std::to_string(mFrameNum) + ".bgeo";
	mGrid.saveDensity(densityFile);

	// Save an image:
	unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];
	for (int i=0; i<recordHeight; i++) 
	{
		glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
			bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
	}
	char anim_filename[2048];
	snprintf(anim_filename, 2048, "../records/smoke_%04d.png", mFrameNum); 
	stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);
	delete [] bitmapData;

	// Dump out rendering particle data in .bgeo file
	std::string particleFile = "../records/frame" + std::to_string(mFrameNum) + ".bgeo";
	mGrid.saveParticle(particleFile);

	mFrameNum++;
}

int SmokeSim::getTotalFrames() {
	return mTotalFrameNum;
}
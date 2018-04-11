
#ifndef smokeSim_H_
#define smokeSim_H_

#include "mac_grid.h"
#include <Partio.h>

class Camera;
class SmokeSim
{
public:
   SmokeSim();
   virtual ~SmokeSim();

   virtual void reset();
   virtual void step();
   virtual void draw(const Camera& c);
   //virtual void setGridDimensions(int x, int y, int z); 
   virtual void setRecording(bool on, int width, int height);
   virtual bool isRecording();
	
	
	int getTotalFrames();

protected:
   virtual void drawAxes();
   virtual void grabScreen();

protected:
	MACGrid mGrid;
	bool mRecordEnabled;
	int mFrameNum;
	int mTotalFrameNum; 
	
	
	int recordWidth;
	int recordHeight;
};

#endif
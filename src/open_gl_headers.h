

#pragma once

// WINDOWS:
// #include <windows.h>
// #include "../../GL/GL.h"
// #include "../../GL/GLU.h"
// #include "../../GL/glut.h"

#ifdef __linux__
    #include <GL/gl.h>
    #include <GL/glu.h>
    #include <GL/glut.h>
#elif _WIN32
    #include "../../GL/GL.h"
    #include "../../GL/GLU.h"
    #include "../../GL/glut.h"
#elif __APPLE__
    #include <OpenGL/gl.h>
    #include <OpenGL/glu.h>
    #include <GLUT/glut.h>
#endif
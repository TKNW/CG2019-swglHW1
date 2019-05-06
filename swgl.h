#ifndef __swgl_h__
#define __swgl_h__

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//implement the following function call for the opengl pipeline


//hw1
void swTranslated(GLdouble x, GLdouble y, GLdouble z);//done
void swScaled(GLdouble x, GLdouble y, GLdouble z);//done
void swRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z);//done

void swMatrixMode(GLenum mode);//done
void swLoadIdentity(void);//done
void swLoadMatrixd(const GLdouble * m);//done
void swMultMatrixd(const GLdouble * m);//done

void swPushMatrix(void);
void swPopMatrix(void);

void swuLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
 	           GLdouble centerX, GLdouble centerY, GLdouble centerZ,
 	           GLdouble upX, GLdouble upY, GLdouble upZ);//done
void swFrustum(	GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble nearVal, GLdouble farVal);
void swuPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);//done

void swViewport(GLint x, GLint y, GLsizei width, GLsizei height);//done

bool swTransformation(const GLdouble h[4], GLdouble w[4]);//done


#endif                  /* __swgl_h__ */

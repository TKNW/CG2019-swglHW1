#include "swgl.h"
#include <cmath>
#include <stack>
using namespace std;
GLdouble CTM_MV[16];	//Current Transformation Matrix: ModelView
GLdouble CTM_P[16];		//Current Transformation Matrix: Projection
GLdouble *CTM;			//Pointer to Current Transformation Matrix
int viewport[4];
struct Matrix{
    GLdouble P[16];
};
stack <Matrix> MatrixStack;
//h: input the vertex in world space(h)
//w: vertex in windows space(w)
GLdouble* matrix4X4(GLdouble *M1,const GLdouble *M2){
    GLdouble *R = new GLdouble[16];
    for(int i = 0;i< 16;i++){
        for(int j = 0; j < 4; j++){
            R[i] += M1[j + (i / 4)*4] * M2[j*4 + i % 4];
        }
    }
    return R;
}
//step1
void swMultMatrixd(const GLdouble * m){
    GLdouble *R = new GLdouble[16];
    for(int i = 0;i < 16;i++){
        for(int j = 0; j < 4; j++){
            R[i] += CTM[j + (i / 4)*4] * m[j*4 + i % 4];
        }
    }
    for(int i = 0;i < 16;i++)
        CTM[i] = R[i];
}
bool swTransformation(const GLdouble h[4], GLdouble w[4]){
	//p = CTM_P*CTM_MV*h
	GLdouble *temp = matrix4X4(CTM_P,CTM_MV);
	GLdouble p[4];
	p[0] = temp[0]*h[0] + temp[1]*h[1] + temp[2]*h[2] + temp[3]*h[3];
    p[1] = temp[4]*h[0] + temp[5]*h[1] + temp[6]*h[2] + temp[7]*h[3];
    p[2] = temp[8]*h[0] + temp[9]*h[1] + temp[10]*h[2] + temp[11]*h[3];
    p[3] = temp[12]*h[0] + temp[13]*h[1] + temp[14]*h[2] + temp[15]*h[3];
	//prespective division
	for(int i = 0; i< 4; i++)
        p[i] = p[i]/p[3];
	//viewport transformation
	w[0] = (p[0] + 1)*(viewport[2]/2) + viewport[0];
	w[1] = (p[1] + 1)*(viewport[3]/2) + viewport[1];
	w[2] = p[2];
	w[3] = p[3];
	return true;
}
void swViewport(GLint x, GLint y, GLsizei width, GLsizei height){
    viewport[0] = x;
    viewport[1] = y;
    viewport[2] = width;
    viewport[3] = height;
}
void swMatrixMode(GLenum mode){
    if(mode == GL_PROJECTION)
        CTM = CTM_P;
    else if(mode == GL_MODELVIEW)
        CTM = CTM_MV;

}
void swLoadIdentity(void){
    for(int i = 0;i < 16;i++){
        if(i == 0 || i == 5 || i == 10 || i == 15)
            CTM[i] = 1;
        else
            CTM[i] = 0;
    }
}
//step2

void swuPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar){
    GLdouble P[16];
    for(int i = 0; i < 16;i++)
        P[i] = 0;
    P[0] = 1/tan((fovy/2) * (2*acos(0)/180)) / aspect;
    P[5] = 1/tan((fovy/2) * (2*acos(0)/180));
    P[10] = (zFar + zNear)/(zNear - zFar);
    P[11] = (2 * zFar * zNear)/(zNear - zFar);
    P[14] = -1;
    swMultMatrixd(P);
}
GLdouble* cross(GLdouble* v1, GLdouble* v2) {
    GLdouble *R = new GLdouble[3];
	R[0] = v1[1]*v2[2] - v1[2]*v2[1];
	R[1] = v1[2]*v2[0] - v1[0]*v2[2];
	R[2] = v1[0]*v2[1] - v1[1]*v2[0];
	return R;
}
void swuLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
 	           GLdouble centerX, GLdouble centerY, GLdouble centerZ,
 	           GLdouble upX, GLdouble upY, GLdouble upZ){
    GLdouble F[3] = {centerX - eyeX,centerY - eyeY,centerZ - eyeZ};
    GLdouble UP[3] = {upX,upY,upZ};
    GLdouble P[16];
    for(int i = 0;i < 16; i++){
        P[i] = 0;
    }
    double Fnor = sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);
    double UPnor = sqrt(UP[0]*UP[0] + UP[1]*UP[1] + UP[2]*UP[2]);
    for(int i = 0; i < 3; i++){
        F[i] = F[i]/Fnor;
        UP[i] = UP[i]/UPnor;
    }
    GLdouble* s = cross(F,UP);
    GLdouble* SN = s;
    GLdouble Snor = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
    for(int i = 0; i < 3; i++){
        SN[i] = SN[i]/Snor;
    }
    GLdouble* u = cross(SN,F);
    for(int i = 0; i < 3;i++){
        P[i] = s[i];
        P[i+4] = u[i];
        P[i+8] = -F[i];
    }
    P[15] = 1;
    swMultMatrixd(P);
    swTranslated(-eyeX, -eyeY, -eyeZ);
}
void swTranslated(GLdouble x, GLdouble y, GLdouble z){
    GLdouble P[16];
    for(int i = 0;i < 16; i++){
        P[i] = 0;
    }
    P[0] = 1;
    P[3] = x;
    P[5] = 1;
    P[7] = y;
    P[10] = 1;
    P[11] = z;
    P[15] = 1;
    swMultMatrixd(P);
}
//step3
void swScaled(GLdouble x, GLdouble y, GLdouble z){
    GLdouble P[16];
    for(int i = 0;i < 16; i++){
        P[i] = 0;
    }
    P[0] = x;
    P[5] = y;
    P[10] = z;
    P[15] = 1;
    swMultMatrixd(P);
}
void swRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z){
    GLdouble P[16];
    for(int i = 0; i < 16;i++)
        P[i] = 0;
    GLdouble c = cos(angle * 2*acos(0)/180);
    GLdouble s = sin(angle * 2*acos(0)/180);
    GLdouble nor = sqrt(x*x + y*y + z*z);
    x = x/nor;
    y = y/nor;
    z = z/nor;
    P[0] = x*x*(1 - c) + c;
    P[1] = x*y*(1 - c) - z*s;
    P[2] = x*z*(1 - c) + y*s;
    P[4] = y*x*(1 - c) + z*s;
    P[5] = y*y*(1 - c) + c;
    P[6] = y*z*(1 - c) - x*s;
    P[8] = x*z*(1 - c) - y*s;
    P[9] = y*z*(1 - c) + x*s;
    P[10] = z*z*(1 - c) + c;
    P[15] = 1;
    swMultMatrixd(P);
}
//step4
void swLoadMatrixd(const GLdouble * m){
    for(int i = 0; i < 16;i++){
        CTM[i] = m[i];
    }
}
void swPushMatrix(void){
    Matrix temp;
    for(int i = 0; i < 16;i++){
        temp.P[i] = CTM[i];
    }
    MatrixStack.push(temp);
}
void swPopMatrix(void){
    Matrix temp = MatrixStack.top();
    MatrixStack.pop();
    swLoadMatrixd(temp.P);
}



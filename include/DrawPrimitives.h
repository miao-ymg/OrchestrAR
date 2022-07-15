#include <GLFW/glfw3.h>


#include <math.h>


/* PI */
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


void drawSphere(double r, int lats, int longs) {
	int i, j;
	for(i = 0; i <= lats; i++) {
		double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
		double z0  = r * sin(lat0);
		double zr0 = r *  cos(lat0);

		double lat1 = M_PI * (-0.5 + (double) i / lats);
		double z1  = r * sin(lat1);
		double zr1 = r * cos(lat1);

		glBegin(GL_QUAD_STRIP);
		for(j = 0; j <= longs; j++) {
			double lng = 2 * M_PI * (double) (j - 1) / longs;
			double x = cos(lng);
			double y = sin(lng);

			glNormal3f(x * zr0, y * zr0, z0);
			glVertex3f(x * zr0, y * zr0, z0);
			glNormal3f(x * zr1, y * zr1, z1);
			glVertex3f(x * zr1, y * zr1, z1);
		}
		glEnd();
	}
}


void drawCone(GLdouble base, GLdouble height, GLint slices, GLint stacks) {

	// draw the upper part of the cone
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0, 0, height);
	for (int angle = 0; angle < 360; angle++) {
		glVertex3f(sin((double)angle) * base, cos((double)angle) * base, 0.f);
	}
	glEnd();

	// draw the base of the cone
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0, 0, 0);
	for (int angle = 0; angle < 360; angle++) {
		// normal is just pointing down
		glNormal3f(0, -1, 0);
		glVertex3f(sin((double)angle) * base, cos((double)angle) * base, 0.f);
	}
	glEnd();
}


void drawCircle(float r, int segments) {
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0, 0, 0);
	for (size_t n = 0; n <= segments; ++n) {
		const float t = 2 * M_PI * (float) n / (float) segments;
		glVertex3f(cos(t) * r, sin(t) * r, 0);
	}
	glEnd();
}


void drawRectangle(float x1, float z1, float x2, float z2) {
	glBegin(GL_QUADS);
	glVertex3f(x1, 0, z1);
	glVertex3f(x1, 0, z2);
	glVertex3f(x2, 0, z2);
	glVertex3f(x2, 0, z1);
	glEnd();
}


void drawCylinder(GLfloat radius, GLfloat height, GLubyte R, GLubyte G, GLubyte B) {
    GLfloat x = 0.0;
    GLfloat y = 0.0;
    GLfloat angle = 30;
    GLfloat angle_stepsize = 0.1;

    /** Draw the tube */
    glColor4f(R - 0.1, G - 0.1, B - 0.1, 1);
    glBegin(GL_QUAD_STRIP);
    angle = 0.0;
        while( angle < 2*M_PI ) {
            x = radius * cos(angle);
            y = radius * sin(angle);
            glVertex3f(x, y , height);
            glVertex3f(x, y , 0.0);
            angle = angle + angle_stepsize;
        }
        glVertex3f(radius, 0.0, height);
        glVertex3f(radius, 0.0, 0.0);
    glEnd();

    /** Draw the circle on top of cylinder */
    glColor4f(R, G, B, 1);
    glBegin(GL_POLYGON);
    angle = 0.0;
        while( angle < 2*M_PI ) {
            x = radius * cos(angle);
            y = radius * sin(angle);
            glVertex3f(x, y , height);
            angle = angle + angle_stepsize;
        }
        glVertex3f(radius, 0.0, height);
    glEnd();
}

void drawEllipse(float rx, float ry, int num_segments) { 
    float theta = 2 * M_PI / float(num_segments); 
    float c = cosf(theta);//precalculate the sine and cosine
    float s = sinf(theta);
    float t;

    float x = 1;//we start at angle = 0 
    float y = 0; 

    glBegin(GL_POLYGON); 
    for(int ii = 0; ii < num_segments; ii++) 
    { 
        //apply radius and offset
        glVertex3f(x * rx, 0, y * ry);//output vertex 

        //apply the rotation matrix
        t = x;
        x = c * x - s * y;
        y = s * t + c * y;
    } 
    glEnd(); 
}

void drawMusicNote() {
	glRotatef(-25, 0, 1, 0);
    drawEllipse(0.14, 0.1, 32);
    glRotatef(25, 0, 1, 0);
    drawRectangle(0.1, 0, 0.13, 0.7);
}
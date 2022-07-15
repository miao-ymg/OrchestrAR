#include <GLFW/glfw3.h>


#include <math.h>


/* PI */
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

enum alignment{XY, XZ, YZ};


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


void drawCircle(GLfloat r, int segments) {
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0, 0, 0);
	for (size_t n = 0; n <= segments; ++n) {
		const float t = 2 * M_PI * (float) n / (float) segments;
		glVertex3f(cos(t) * r, sin(t) * r, 0);
	}
	glEnd();
}


void drawAlignedRectangle(GLfloat a1, GLfloat b1, GLfloat a2, GLfloat b2, alignment al) {
    glBegin(GL_QUADS);
	switch (al) {
		case XY:
			glVertex3f(a1, b1, 0);
			glVertex3f(a1, b2, 0);
			glVertex3f(a2, b2, 0);
			glVertex3f(a2, b1, 0);
			break;
		case XZ:
			glVertex3f(a1, 0, b1);
			glVertex3f(a1, 0, b2);
			glVertex3f(a2, 0, b2);
			glVertex3f(a2, 0, b1);
			break;
		case YZ:
			glVertex3f(0, a1, b1);
			glVertex3f(0, a1, b2);
			glVertex3f(0, a2, b2);
			glVertex3f(0, a2, b1);
			break;
		default:
			throw std::invalid_argument("ERROR: Invalid alignment!");
			break;
	}
	glEnd();
}


void drawCylinder(GLfloat radius, GLfloat height, GLfloat R, GLfloat G, GLfloat B) {
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
			glVertex3f(x, y , 0);
            glVertex3f(x, y , height);
            angle = angle + angle_stepsize;
        }
        glVertex3f(radius, 0.0, height);
    glEnd();
}

void drawEllipse(GLfloat rx, GLfloat ry, int num_segments) { 
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
    drawAlignedRectangle(0.09, 0.0, 0.135, 0.7, XZ);
}

void drawCuboid(GLfloat l, GLfloat w, GLfloat h, GLfloat R, GLfloat G, GLfloat B) {
	glColor4f(R, G, B, 1);
    glBegin(GL_QUADS);
	//bottom
	drawAlignedRectangle(-l/2, -w/2, l/2, w/2, XY);
	glPushMatrix();
	//sides
	glColor4f(R - 0.1, G - 0.1, B - 0.1, 1);
	glTranslatef(0, w/2, 0);
	drawAlignedRectangle(-l/2, 0, l/2, h, XZ);
	glTranslatef(0, -w, 0);
	drawAlignedRectangle(l/2, 0, -l/2, h, XZ);
	glTranslatef(0, w/2, 0);

	glRotatef(90, 0, 0, 1);
	glTranslatef(0, l/2, 0);
	drawAlignedRectangle(w/2, 0, -w/2, h, XZ);
	glTranslatef(0, -l, 0);
	drawAlignedRectangle(w/2, 0, -w/2, h, XZ);

	//top
	glColor4f(R, G, B, 1);
	glPopMatrix();
	glTranslatef(0, 0, h);
	drawAlignedRectangle(l/2, -w/2, -l/2, w/2, XY);
	glEnd();
}

void drawSpeakers(GLfloat l, GLfloat w, GLfloat h, GLfloat R, GLfloat G, GLfloat B){
	glPushMatrix();
	drawCuboid(l, w, h, R, G, B);
	glPopMatrix();
	//membranes 1 and 2, 1 half the size of 2
	GLfloat r1 = fmin(w/2, h/4) * 0.9;
	GLfloat r2 = fmin((r1/2) * 1.2, h/4);

	//membrane radius is periodic
	float sizeFactor = 0.95 + (sin(12*glfwGetTime()) * 0.05);
	float r1_ = r1 * sizeFactor;
	float r2_ = r2 * sizeFactor;

	GLfloat z1 = r1;
	GLfloat z2 = h/2 + r2;

	glColor4f(0.0, 0.0, 0.0, 1.0);
	glTranslatef(0, l/2, 0);
	glTranslatef(0, 0, z1);
	glRotatef(90, -1, 0, 0);
	drawCylinder(r1_, 0.02, 1, 1, 1);
	glTranslatef(0, -1*(2*r1 + r2), 0);
	drawCylinder(r2_, 0.02, 1, 1, 1);
}

void drawPianoKeys(GLfloat l, GLfloat w, GLfloat h){
	float r = 1.0;
	float g = 1.0;
	float b = 1.0;

	//black keys are smaller
	GLfloat bLen = l * 0.8;
	GLfloat bWid = w * 0.5;

	//lower white row
	glPushMatrix();
	drawCuboid(l, w, h, r, g, b);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, w+0.01, 0);
	drawCuboid(l, w, h, r, g, b);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, -w-0.01, 0);
	drawCuboid(l, w, h, r, g, b);
	glPopMatrix();

	r = 0.1;
	g = 0.1;
	b = 0.1;

	//upper black row
	glTranslatef(l * 0.2, 0, h*0.5);
	glPushMatrix();
	glTranslatef(0, w*0.5, 0);
	drawCuboid(bLen, bWid, h, r, g, b);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, -w*0.5, 0);
	drawCuboid(bLen, bWid, h, r, g, b);
	glPopMatrix();

}
#include <GLFW/glfw3.h>
#include <math.h>
#include <algorithm>

#include "Instrument.h"


/* PI */
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


enum plane{XY, XZ, YZ};


float colorIntensity(int pitch) {
    float x = pitch % 24;
    return (x < 12 ? fmin(1.0/4 * x, 1.0) : fmax(-1.0/4 * x + 4, 0.0));
}


// Paint an object depending on the sample's pitch
std::array<float, 3> paintObject(Pitch pitch) {
    // Sample has no pitch -> Custom color
    if (pitch == 24)
        return {0.7, 0.7, 0.7};

    float r = colorIntensity(pitch + 8);
    float g = colorIntensity(pitch);
    float b = colorIntensity(pitch - 8);

    return {r, g, b};
}


// --- BASIC SHAPES ---

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


/**
 * Draws an axis aligned rectangle which is parallel to the plane given in pl, enclosed by two vertices v, w
 *
 * @param v1, v2, w1, w2 Coordinates of the two vertices v = (v1, v2) and w = (w1, w2)
 * @param pl Plane which is parallel to the rectangle
 */
void drawAlignedRectangle(GLfloat v1, GLfloat v2, GLfloat w1, GLfloat w2, plane pl) {
    glBegin(GL_QUADS);
	switch (pl) {
		case XY:
			glVertex3f(v1, v2, 0);
			glVertex3f(v1, w2, 0);
			glVertex3f(w1, w2, 0);
			glVertex3f(w1, v2, 0);
			break;
		case XZ:
			glVertex3f(v1, 0, v2);
			glVertex3f(v1, 0, w2);
			glVertex3f(w1, 0, w2);
			glVertex3f(w1, 0, v2);
			break;
		case YZ:
			glVertex3f(0, v1, v2);
			glVertex3f(0, v1, w2);
			glVertex3f(0, v2, w2);
			glVertex3f(0, v2, v2);
			break;
		default:
			throw std::invalid_argument("ERROR: Invalid alignment!");
	}
	glEnd();
}


/**
 * Draws a cylinder
 *
 * @param radius Radius of the cylinder
 * @param height Height of the cylinder
 * @param r, g, b Values defining the color code
 */
void drawCylinder(GLfloat radius, GLfloat height, GLfloat r, GLfloat g, GLfloat b) {
    GLfloat x, y, angle;
    GLfloat angle_stepsize = 0.1;

    // Draw the tube
    glColor4f(r - 0.1, g - 0.1, b - 0.1, 1);
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

    // Draw the circles at the top and bottom of the cylinder
    glColor4f(r, g, b, 1);
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


/**
 * Draws an ellipse
 *
 * @param rx, ry Length of the radiuses for both x and y direction
 * @param num_segments Number of straight lines the outer ellipse line consists of
 */
void drawEllipse(GLfloat rx, GLfloat ry, int num_segments) { 
    float theta = 2 * M_PI / float(num_segments); 
    float c = cosf(theta);//precalculate the sine and cosine
    float s = sinf(theta);
    float t;

    float x = 1;//we start at angle = 0 
    float y = 0; 

    glBegin(GL_POLYGON); 
    for(int i = 0; i < num_segments; ++i) {
        // Apply radius and offset
        glVertex3f(x * rx, 0, y * ry);//output vertex 

        // Apply the rotation matrix
        t = x;
        x = c * x - s * y;
        y = s * t + c * y;
    } 
    glEnd(); 
}


/**
 * Draws a cuboid
 *
 * @param l, w, h Length, width and height of the cuboid
 * @param r, g, b Values defining the color code
 */
void drawCuboid(GLfloat l, GLfloat w, GLfloat h, GLfloat r, GLfloat g, GLfloat b) {
    glColor4f(r, g, b, 1);
    glBegin(GL_QUADS);
    // Bottom
    drawAlignedRectangle(-l/2, -w/2, l/2, w/2, XY);
    glPushMatrix();
    // Sides
    glColor4f(r - 0.1, g - 0.1, b - 0.1, 1);
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
    // Top
    glColor4f(r, g, b, 1);
    glPopMatrix();
    glTranslatef(0, 0, h);
    drawAlignedRectangle(l/2, -w/2, -l/2, w/2, XY);
    glEnd();
}



// --- COMPLEX OBJECTS ---

// Draws a speaker membrane with given radius r
void drawSpeakerMembrane(GLfloat r) {
    drawCylinder(r, 0.02, 1, 1, 1);
    glTranslatef(0, 0, 0.01);
    drawCylinder(r / 2, 0.01, 0.3, 0.3, 0.3);
    glTranslatef(0, 0, -0.01);
}


/**
 * Draws a speaker with a vibrating membrane for BASS samples
 *
 * @param l, w, h Length, width and height of the speaker
 * @param r, g, b Values defining the color code
 */
void drawSpeaker(GLfloat l, GLfloat w, GLfloat h, GLfloat r, GLfloat g, GLfloat b) {
    glScalef(1.4, 1.4, 1.4);
    // Draw the box
	glPushMatrix();
	drawCuboid(l, w, h, r, g, b);
	glPopMatrix();

    // For the membrane vibration
	GLfloat r1 = fmin(w/2, h/4);
	GLfloat r2 = fmin(0.75 * r1, h/4);
	float sizeFactor = 0.75 + (sin(12*glfwGetTime()) * 0.05);
	float r1_vib = r1 * sizeFactor;
	float r2_vib = r2 * sizeFactor;

    // Draw the individual speakers
	glTranslatef(0, l/2, 0);
	glTranslatef(0, 0, 1.2 * r1);
	glRotatef(90, -1, 0, 0);
    drawSpeakerMembrane(r1_vib);
	glTranslatef(0, -(r1 + 1.1 * r2), 0);
    drawSpeakerMembrane(r2_vib);
}


/**
 * Draws a drum for BEATS samples
 *
 * @param height, width Height and width of the drum
 * @param r, g, b Values defining the color code
 */
void drawDrum(GLfloat height, GLfloat width, GLfloat r, GLfloat g, GLfloat b) {
    for (size_t i = 0; i < 6; i++) {
        glTranslatef(width * 2 * M_PI * sin(i) / 6, width * 2 * M_PI * cos(i) / 6, 0);
        drawCylinder(0.03, 0.7, 1, 1, 1);
        glTranslatef(width * -2 * M_PI * sin(i) / 6, width * -2 * M_PI * cos(i) / 6, 0);
    }
    // Drum shell
    glTranslatef(0, 0, 0.2);
    drawCylinder(width, height, r, g, b);
    // Drum skin
    glTranslatef(0, 0, height - 0.0);
    drawCylinder(width + 0.05, 0.05, 1, 1, 1);
}


// Draws a musical note for MELODY samples
void drawMusicalNote() {
    glRotatef(-25, 0, 1, 0);
    drawEllipse(0.14, 0.1, 32);
    glRotatef(25, 0, 1, 0);
    drawAlignedRectangle(0.09, 0.0, 0.135, 0.7, XZ);
}


/**
 * Draws a set of piano keys for KEYS samples
 *
 * @param l, w, h Length, width and height of the individual white keys
 * @param r, g, b Values defining the color code
 */
void drawPianoKeys(GLfloat l, GLfloat w, GLfloat h, GLfloat r, GLfloat g, GLfloat b) {
    // Dimensions of black keys
    GLfloat bLen = l * 0.65;
    GLfloat bWid = w * 0.5;

    // "White" keys
    glPushMatrix();
    drawCuboid(l, w, h, r, g, b);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0, w + 0.01, 0);
    drawCuboid(l, w, h, r, g, b);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0, -w - 0.01, 0);
    drawCuboid(l, w, h, r, g, b);
    glPopMatrix();

    // Black keys
    glTranslatef(0.15 * l, 0, h / 2);
    glPushMatrix();
    glTranslatef(0, w / 2, 0);
    drawCuboid(bLen, bWid, h, 0.1, 0.1, 0.1);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0, -w / 2, 0);
    drawCuboid(bLen, bWid, h, 0.1, 0.1, 0.1);
    glPopMatrix();
}


/**
 * Draws a microphone for VOCAL samples
 *
 * @param r, g, b Values defining the color code
 */
void drawMicrophone(GLfloat r, GLfloat g, GLfloat b) {
    glRotatef(15, 1, 0, 0);
    drawCylinder(0.15, 1.5, r, g, b);
    glTranslatef(0, 0, 1.4);
    drawCylinder(0.3, 0.04, r, g, b);
    glColor4f(0.2, 0.2, 0.2, 1);
    drawSphere(0.3, 10, 10);
    glTranslatef(0, 0, -2);
}

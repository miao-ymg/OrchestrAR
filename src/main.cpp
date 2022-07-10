#include <iostream>
#include <array>
#include <unordered_map>
#include <cmath>
#include <iomanip>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h>

#include <GLFW/glfw3.h>
#include <SDL.h>
#include <SDL_mixer.h>
// #include <cpct_dll.h>

#include "glext.h"
#include "DrawPrimitives.h"

//#include "../opencv2/opencv.hpp"
//#include "../include/ownlib.h"

#include "Instrument.cpp"
#include "MarkerTracker.cpp"
#include "PoseEstimation.cpp"


using namespace cv;
using namespace std;

#define TEST 0
#define EX5 1
#define EX5_RAW 0
#define DRAW_CONTOUR 0
#define DRAW_RECTANGLE 0

#define THICKNESS_VALUE 4


// --- Audio properties ---
const int audio_rate = 22050;
const Uint16 audio_format = AUDIO_S16SYS;
const int audio_channels = 4;
const int audio_buffers = 4096;

cv::VideoCapture cap;

// Camera settings
int camera_width  = 1280; 
int camera_height = 720; 
const int virtual_camera_angle = 30;

// queue<vector<int> > identifierHistory;
unordered_map<int, Instrument> instruments;

/*
int initSDL() {
	//initialize SDL
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) != 0) {
		fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
		return 1;
	}
	if (Mix_OpenAudio(audio_rate, audio_format, audio_channels, audio_buffers) != 0) {
		fprintf(stderr, "Unable to initialize audio: %s\n", Mix_GetError());
		return 1;
	}

	// --- Load sound files ---
	for (auto& instr : instruments) {
		instr.second.loadSound();
	}

	return 0;
}
*/

void initVideoStream(cv::VideoCapture &cap) {
	if(cap.isOpened())
		cap.release();

	cap.open(0);
	camera_width  = cap.get(CV_CAP_PROP_FRAME_WIDTH); 
	camera_height = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
	cout << "Width: " << camera_width << endl;
	cout << "Height: " << camera_height << endl;

	if (cap.isOpened() == false) {
		std::cout << "No webcam found, using a video file" << std::endl;
		cap.open("MarkerMovie.mpg");
		if (cap.isOpened() == false) {
			std::cout << "No video file found. Exiting." << std::endl;
			exit(0);
		}
	} else {
		//get camera resolution
		camera_height = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
		camera_width = cap.get(CV_CAP_PROP_FRAME_WIDTH);

		cout << "res width: " << camera_width << endl;
		cout << "res height: " << camera_height << endl;
	}
}


/* Program & OpenGL initialization */
void initGL(int argc, char *argv[]) {

// Added in Exercise 8 - End *****************************************************************

	// For our connection between OpenCV/OpenGL
    // Pixel storage/packing stuff -> how to handle the pixel on the graphics card
	// For glReadPixels​ -> Pixel representation in the frame buffer
    glPixelStorei(GL_PACK_ALIGNMENT,   1);
	// For glTexImage2D​ -> Define the texture image representation
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	// Turn the texture coordinates from OpenCV to the texture coordinates OpenGL
    glPixelZoom(1.0, -1.0);

// Added in Exercise 8 - End *****************************************************************

    // Enable and set colors
    glEnable(GL_COLOR_MATERIAL);
    glClearColor(0, 0, 0, 1.0);

	// Enable and set depth parameters
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0);

	// Light parameters
	GLfloat light_amb[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat light_pos[] = {1.0, 1.0, 1.0, 0.0};
    GLfloat light_dif[] = {0.7, 0.7, 0.7, 1.0};

	// Enable lighting
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_amb);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_dif);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}

void display(GLFWwindow* window, const cv::Mat &img_bgr, vector<Instrument> &visibleInstruments) {

// Added in Exercise 8 - Start *****************************************************************

	const float degreePerSec = 90.0f;

	const float angle = (float) glfwGetTime() * degreePerSec;

	unsigned char bkgnd[camera_width * camera_height * 3];

	// Copy picture data into bkgnd array
	memcpy(bkgnd, img_bgr.data, sizeof(bkgnd));

// Added in Exercise 8 - End *****************************************************************

	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Needed for rendering the real camera image
	glMatrixMode(GL_MODELVIEW);
	// No position changes
	glLoadIdentity();
	
// Added in Exercise 8 - Start *****************************************************************

	glDisable(GL_DEPTH_TEST);

	glMatrixMode(GL_PROJECTION);
	// Push the projection matrix (frustum) -> frustum will be saved on the stack
	glPushMatrix();
	glLoadIdentity();
	// In the ortho view all objects stay the same size at every distance
	glOrtho(0.0, camera_width, 0.0, camera_height,-1,1);

	// -> Render the camera picture as background texture
	// Making a raster of the image -> -1 otherwise overflow
	glRasterPos2i(0, camera_height-1);
	// Load and render the camera image -> unsigned byte because of bkgnd as unsigned char array
	// bkgnd 3 channels -> pixelwise rendering
	glDrawPixels(camera_width, camera_height, GL_BGR_EXT, GL_UNSIGNED_BYTE, bkgnd);

	// Go back to the previous projection -> frustum
	glPopMatrix();

	// Activate depth -> that snowman can be scaled with depth
	glEnable(GL_DEPTH_TEST);

	// Move to marker-position
	glMatrixMode(GL_MODELVIEW);

	for (auto& instr : visibleInstruments) {
		
		// Sadly doesn't work for Windows -> so we made own solution!
		//glLoadTransposeMatrixf(resultMatrix);

		std::array<float, 16> poseMatrix = instr.getPoseMatrix();

		// -> Transpose the Modelview Matrix
		float resultTransposedMatrix[16];
		for (int x=0; x<4; ++x) {
			for (int y=0; y<4; ++y) {
				// Change columns to rows
				resultTransposedMatrix[x*4+y] = poseMatrix[y*4+x];
			}
		}

		// Load the transpose matrix
		glLoadMatrixf(resultTransposedMatrix);
		
		// Rotate 90 desgress in x-direction
		glRotatef(-90, 0, 0, 1);
		// Scale down!
		glScalef(0.03, 0.03, 0.03);

	// Added in Exercise 8 - End *****************************************************************

		glRotatef(angle, 0, 0, 1);

		// Make alpha value depending on volume (between 0 and 1)
		// float alpha = (float) instr.getVolume() / 128;

		switch (instr.getRole()) {
			case BASS:
				glColor4f(1.0, 0.0, 0.0, 1);
				break;
			case BEAT:
				glColor4f(0.0, 1.0, 0.0, 1);
				break;
			case KEYS:
				glColor4f(0.0, 0.5, 0.5, 1);
				break;
			case MELODY:
				glColor4f(0.0, 0.0, 1.0, 1);
				break;
			case VOCAL:
				glColor4f(1.0, 1.0, 0.0, 1);
				break;
			default:
				glColor4f(1.0, 1.0, 1.0, 1);
				break;
		}

		drawCircle(1, 6);
		drawSphere(0.5, 10, 10);
		// Draw 3 white spheres
		/*
		glColor4f(1.0, 1.0, 1.0, 1.0);
		drawSphere(0.8, 10, 10);
		glTranslatef(0.0, 0.8, 0.0);
		drawSphere(0.6, 10, 10);
		glTranslatef(0.0, 0.6, 0.0);
		drawSphere(0.4, 10, 10);

		// Draw the eyes
		glPushMatrix();
		glColor4f(0.0, 0.0, 0.0, 1.0);
		glTranslatef(0.2, 0.2, 0.2);
		drawSphere(0.066, 10, 10);
		glTranslatef(0, 0, -0.4);
		drawSphere(0.066, 10, 10);
		glPopMatrix();

		// Draw a nose
		glColor4f(1.0, 0.5, 0.0, 1.0);
		glTranslatef(0.3, 0.0, 0.0);
		glRotatef(90, 0, 1, 0);
		drawCone(0.1, 0.3, 10, 10);
		*/
	}
}

void reshape( GLFWwindow* window, int width, int height ) {
	// Set a whole-window viewport
	glViewport( 0, 0, (GLsizei)width, (GLsizei)height );

	// Create a perspective projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// The camera should be calibrated -> a calibration results in the projection matrix -> then load the matrix
	// -> into GL_PROJECTION
	// -> adjustment of FOV is needed for each camera
	float ratio = (GLfloat)width / (GLfloat)height;
	
	float near = 0.01f, far = 100.f;
	float top = tan((double)(virtual_camera_angle*M_PI / 360.0f)) * near;
	float bottom = -top;
	float left = ratio * bottom;
	float right = ratio * top;
	glFrustum(left, right, bottom, top, near, far);
}


int main(int argc, char* argv[]) {

	// --- Add instrument markers ---
	
	// --- Orchestra 1 ---
	instruments.insert(std::make_pair(0x0690, Instrument(0x0690, BASS, "../Sounds/Bass/130_bass_Bm_Ridem-Cowgirl-Mid-Bass.wav")));
	instruments.insert(std::make_pair(0x0272, Instrument(0x0272, BEAT, "../Sounds/Beat/130_beat_Trap-Drum-130bpm.wav")));
	instruments.insert(std::make_pair(0x1c44, Instrument(0x1c44, MELODY, "../Sounds/Melody/130_melody_D_Paris-Emotional-Piano-Loop.wav")));
	// instruments.insert(std::make_pair(0x005a, Instrument(0x005a, VOCAL, "../Sounds/Vocal/130_vocal_D_Emotions.wav")));
	
	/*
	// --- Orchestra 2 ---
	instruments.insert(std::make_pair(0x0690, Instrument(0x0690, BASS, "../Sounds/Bass/130_bass_Em_LilTecca-LilMosey-Type-Melody-Part-4.wav")));
	instruments.insert(std::make_pair(0x0272, Instrument(0x0272, BEAT, "../Sounds/Beat/130_beat_Trap-Drum-130bpm.wav")));
	instruments.insert(std::make_pair(0x1c44, Instrument(0x1c44, MELODY, "../Sounds/Melody/130_melody_Em_Piano-YXNG-SXN.wav")));
	*/

	GLFWwindow* window;

	// Initialize the library
	if (!glfwInit())
		return -1;

	// --- Initialize SDL ---
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) != 0) {
		fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
		return 1;
	}
	if (Mix_OpenAudio(audio_rate, audio_format, audio_channels, audio_buffers) != 0) {
		fprintf(stderr, "Unable to initialize audio: %s\n", Mix_GetError());
		return 1;
	}

	// --- Load sound files ---
	for (auto& instr : instruments) {
		instr.second.loadSound();
	}

	// const GLFWvidmode* mode = glfwGetVideoMode(NULL);
	
	// glfwWindowHint(GLFW_R);

	// Initialize the window system
	// Create a windowed mode window and its OpenGL context
	window = glfwCreateWindow(camera_width, camera_height, "OrchestrAR", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}

	//glfwSetWindowSize(window, camera_width, camera_height);
	glfwSetWindowAspectRatio(window, camera_width, camera_height);
	//glfwSetWindowSizeLimits(window, camera_width, camera_height, camera_width, camera_height);
	
	// Set callback functions for GLFW
	glfwSetFramebufferSizeCallback(window, reshape);

	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);
	
	int window_width, window_height;
	glfwGetFramebufferSize(window, &window_width, &window_height);
	reshape(window, window_width, window_height);

	// Initialize the GL library
	initGL(argc, argv);

    // Setup OpenCV
	cv::Mat img_bgr;
	// Get video stream
	initVideoStream(cap);
	// [m]
	const double kMarkerSize = 0.045;
	// Constructor with the marker size (similar to Exercise 5)
	MarkerTracker markerTracker(kMarkerSize);
	
	float resultMatrix[16];
	// Loop until the user closes the window
	while (!glfwWindowShouldClose(window)) {
		// Capture here
		cap >> img_bgr;
		
		// Make sure that we got a frame -> otherwise crash
		if(img_bgr.empty()) {
			std::cout << "Could not query frame. Trying to reinitialize." << std::endl;
			initVideoStream(cap);
			// Wait for one sec.
			cv::waitKey(1000);
			continue;
		}

		vector<int> identifiers;

		// Track a marker and get the pose of the marker
		markerTracker.findMarker(img_bgr, resultMatrix, instruments, identifiers);

		vector<Instrument> visibleInstruments;

		// ---Play sound if markers are present ---
		for (auto& instr : instruments) {
			instr.second.toggleSound(identifiers);
			if (find(identifiers.begin(), identifiers.end(), instr.second.getID()) != identifiers.end()) {
				visibleInstruments.push_back(instr.second);
			}
		}

		// glViewport( 0, 0, camera_width, camera_height);

		// Render here
		display(window, img_bgr, visibleInstruments);

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();
	}

	for (auto& instr : instruments) {
		instr.second.freeChunk();
	}

	Mix_CloseAudio();
	SDL_Quit();

	// Important -> Avoid memory leaks!
	glfwTerminate();

    return 0;
}

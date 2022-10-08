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

#include "DrawObjects.h"
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
unordered_map<int, Instrument > instruments;


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
	/*
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_amb);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_dif);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	*/
}

void display(GLFWwindow* window, const cv::Mat &img_bgr, vector<Instrument> &visibleInstruments) {

// Added in Exercise 8 - Start *****************************************************************

	const float degreePerSec = 60.0f;

	//models "jump"
	float height = abs(sin((float)glfwGetTime() * 3)) * 0.3;


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
	glDrawPixels(camera_width, camera_height, 0x80E0, GL_UNSIGNED_BYTE, bkgnd);

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
		
		glRotatef(-90, 0, 0, 1);
		glRotatef(180, 1, 0, 0);
		// Scale down!
		glScalef(0.05, 0.05, 0.05);

	// Added in Exercise 8 - End *****************************************************************

		glRotatef(angle, 0, 0, 1);

		glTranslatef(0,0,height);

		instr.drawObject();
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
	bool toggle = true;

	if (toggle) {
		// --- Sound set 1 ---
		instruments.insert(std::make_pair(0x005a, Instrument(0x005a, BEAT, None, "../samples/beat/130-beat-N-Trap_Drum_130bpm.wav")));
		instruments.insert(std::make_pair(0x0690, Instrument(0x0690, BASS, Em, "../samples/bass/130-bass-Em-LilTecca_LilMosey_Type_Melody_Part_4.wav")));
		instruments.insert(std::make_pair(0x0272, Instrument(0x0272, MELODY, Em, "../samples/melody/130-melody-Em-Piano_YXNG_SXN.wav")));
		//instruments.insert(std::make_pair(0x1c44, Instrument(0x1c44, KEYS, Em, "../samples/keys/130-keys-Em-Duel_Of_The_Fates_I_String_Staccato.wav")));
		//instruments.insert(std::make_pair(0x0B44, Instrument(0x0B44, MELODY, Em, "../samples/melody/130-melody-Em-Gunna_Money_Man_BROKEN_By_Danil040.wav")));
		instruments.insert(std::make_pair(0x1228, Instrument(0x1228, VOCAL, Am, "../samples/vocal/130-vocal-Am-voc.wav")));

	} else {
		// --- Sound set 2 ---
		instruments.insert(std::make_pair(0x005a, Instrument(0x005a, BEAT, None, "../samples/beat/130-beat-N-Trap_Drum_130bpm.wav")));
		instruments.insert(std::make_pair(0x0690, Instrument(0x0690, BASS, Bm, "../samples/bass/130-bass-Bm-Ridem_Cowgirl_Mid_Bass.wav")));
		instruments.insert(std::make_pair(0x0272, Instrument(0x0272, MELODY, D, "../samples/melody/130-melody-D-Paris_Emotional_Piano_Loop.wav")));
		//instruments.insert(std::make_pair(0x1228, Instrument(0x1228, VOCAL, D, "../samples/vocal/130-vocal-D-Emotions.wav")));
	}
	
	
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
		instr.second.loadRandomSound(instr.second.getRole());
	}

	// const GLFWvidmode* mode = glfwGetVideoMode(NULL);
	
	// glfwWindowHint(GLFW_R);

	// Initialize the window system
	// Create a windowed mode window and its OpenGL context
	window = glfwCreateWindow(640, 360, "OrchestrAR", NULL, NULL);
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

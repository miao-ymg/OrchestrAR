#include "ownlib.h"

#include "../opencv2/opencv.hpp"
#include <iostream>
#include <array>

#include <SDL.h>
#include <SDL_mixer.h>
#include "Instrument.h"
#include "Instrument.cpp"

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



// Struct holding all infos about each strip, e.g. length
struct MyStrip {
	int stripeLength;
	int nStop;
	int nStart;
	Point2f stripeVecX;
	Point2f stripeVecY;
};

// List of points
typedef vector<Point> contour_t;
// List of contours
typedef vector<contour_t> contour_vector_t;


const int threshold_slider_max = 255;
int threshold_slider = 0;

const int fps = 30;

Mat videoStreamFrameGray;
Mat videoStreamFrameOutput;

const string stripWindow = "Strip Window";
const string kWinName4 = "Exercise 4 - Marker";

// Added in sheet 4 - Start *****************************************************************

int bw_thresh = 55;

// Added in sheet 4 - End *******************************************************************

// Pos is from UI, dereferencing of the pointer
static void on_trackbar(int pos, void* slider_value) {
	*((int*)slider_value) = pos;
	// C++ >= 11 -> Standard
	//*(static_cast<int*>(slider_value)) = pos;
}

// Added in sheet 4 - Start *****************************************************************

void bw_trackbarHandler(int pos, void* slider_value) {
	*((int*)slider_value) = pos;
}

// Added in sheet 4 - End *******************************************************************

int subpixSampleSafe(const Mat& pSrc, const Point2f& p) {
	// Point is float, slide 14
	int fx = int(floorf(p.x));
	int fy = int(floorf(p.y));

	if (fx < 0 || fx >= pSrc.cols - 1 ||
		fy < 0 || fy >= pSrc.rows - 1)
		return 127;

	// Slides 15
	int px = int(256 * (p.x - floorf(p.x)));
	int py = int(256 * (p.y - floorf(p.y)));

	// Here we get the pixel of the starting point
	unsigned char* i = (unsigned char*)((pSrc.data + fy * pSrc.step) + fx);

	// Internsity, shift 3
	int a = i[0] + ((px * (i[1] - i[0])) >> 8);
	i += pSrc.step;
	int b = i[0] + ((px * (i[1] - i[0])) >> 8);

	// We want to return Intensity for the subpixel
	return a + ((py * (b - a)) >> 8);
}

// Added in Sheet 3 - Ex7 (a) Start *****************************************************************
Mat calculate_Stripe(double dx, double dy, MyStrip& st) {
	// Norm (euclidean distance) from the direction vector is the length (derived from the Pythagoras Theorem)
	double diffLength = sqrt(dx * dx + dy * dy);

	// Length proportional to the marker size
	st.stripeLength = (int)(0.8 * diffLength);

	if (st.stripeLength < 5)
		st.stripeLength = 5;

	// Make stripeLength odd (because of the shift in nStop), Example 6: both sides of the strip must have the same length XXXOXXX
	st.stripeLength |= 1;
	/*if (st.stripeLength % 2 == 0)
		st.stripeLength++;*/

		// E.g. stripeLength = 5 --> from -2 to 2: Shift -> half top, the other half bottom
	st.nStop = st.stripeLength >> 1;
	//st.nStop = st.stripeLength / 2;
	st.nStart = -st.nStop;

	Size stripeSize;

	// Sample a strip of width 3 pixels
	stripeSize.width = 3;
	stripeSize.height = st.stripeLength;

	// Normalized direction vector
	st.stripeVecX.x = dx / diffLength;
	st.stripeVecX.y = dy / diffLength;

	// Normalized perpendicular vector
	st.stripeVecY.x = st.stripeVecX.y;
	st.stripeVecY.y = -st.stripeVecX.x;

	// 8 bit unsigned char with 1 channel, gray
	return Mat(stripeSize, CV_8UC1);
}

array<int, 4> getCanonicalId(vector<array<int, 4> > rotations) {
	if (rotations.size() == 0) return { 0,0,0,0 };

	array<int, 4> result = rotations.at(0);
	for (int i = 0; i < rotations.size(); i++) {
		if (rotations[i] < result) result = rotations[i];
	}
	return result;
}

int main(int argc, char **args) {
	Mat frame;
	VideoCapture cap(0);

	queue<vector<int> > identifierHistory;


	vector<Instrument> instruments;

	// --- Add instrument markers
	instruments.push_back(Instrument(1680, "../../Sounds/Melody/130_latinTrap_guitar.wav"));
	instruments.push_back(Instrument(626, "../../Sounds/Drums/130_basicTrap_drums.wav"));


	const string streamWindow = "Stream";

	int capResHeight = 240;
	int capResWidth = 320;

	if (!cap.isOpened()) {
		cout << "No webcam, using video file" << endl;
		cap.open("MarkerMovie.MP4");
		if (cap.isOpened() == false) {
			cout << "No video!" << endl;
			exit(0);
			return -1;
		}
	}
	else {
		//get camera resolution
		capResHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
		capResWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH);

		cout << "res width: " << capResWidth << endl;
		cout << "res height: " << capResHeight << endl;
	}
	
	//initialize SDL
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) != 0) {
		fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
		return 1;
	}
	if (Mix_OpenAudio(audio_rate, audio_format, audio_channels, audio_buffers) != 0) {
		fprintf(stderr, "Unable to initialize audio: %s\n", Mix_GetError());
		exit(1);
	}


	// --- Load sound files ---
	for (Instrument& instr : instruments) {
		instr.loadSound();
	}


	// Added in Sheet 3 - Start *****************************************************************

	bool isFirstStripe = true;

	// Added in Sheet 3 - End *******************************************************************

	// Added in sheet 4 - Start *****************************************************************

	bool isFirstMarker = true;

	// Added in sheet 4 - End *******************************************************************

	const string contoursWindow = "Webcam";
	const string UI = "Threshold";
	namedWindow(contoursWindow, CV_WINDOW_AUTOSIZE);
	// namedWindow(stripWindow, CV_WINDOW_AUTOSIZE);
	const string markerWindow = "Marker";
	// namedWindow(markerWindow, CV_WINDOW_NORMAL);


	int slider_value = 140;
	createTrackbar(UI, contoursWindow, &slider_value, 255, on_trackbar, &slider_value);

	// Added in sheet 4 - Start *****************************************************************

	//int bw_sileder_value = bw_thresh;

	//createTrackbar("BW Threshold", contoursWindow, &bw_sileder_value, 255, bw_trackbarHandler, &bw_sileder_value);

	namedWindow(kWinName4, CV_WINDOW_NORMAL);

	resizeWindow(kWinName4, 120, 120);

	// Added in sheet 4 - End *******************************************************************

	Mat imgFiltered;

		while (cap.read(frame)) {

		vector<int> identifiers;

		// --- Process Frame ---

		Mat grayScale;
		imgFiltered = frame.clone();
		cvtColor(imgFiltered, grayScale, COLOR_BGR2GRAY);

		// Threshold to reduce the noise
		threshold(grayScale, grayScale, slider_value, 255, THRESH_BINARY);

		contour_vector_t contours;

		// RETR_LIST is a list of all found contour, SIMPLE is to just save the begin and ending of each edge which belongs to the contour
		findContours(grayScale, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);

		//drawContours(imgFiltered, contours, -1, Scalar(0, 255, 0), 4, 1);

		// size is always positive, so unsigned int -> size_t; if you have not initialized the vector it is -1, hence crash
		for (size_t k = 0; k < contours.size(); k++) {

			// -------------------------------------------------

			// --- Process Contour ---

			contour_t approx_contour;

			// Simplifying of the contour with the Ramer-Douglas-Peuker Algorithm
			// true -> Only closed contours
			// Approximation of old curve, the difference (epsilon) should not be bigger than: perimeter(->arcLength)*0.02
			approxPolyDP(contours[k], approx_contour, arcLength(contours[k], true) * 0.02, true);

#if DRAW_CONTOUR
			contour_vector_t cov, aprox;
			cov.emplace_back(contours[k]);
			aprox.emplace_back(approx_contour);
			if (approx_contour.size() > 1) {
				drawContours(imgFiltered, aprox, -1, Scalar(0, 255, 0), 4, 1);
				drawContours(imgFiltered, cov, -1, Scalar(255, 0, 0), 4, 1);
				continue;
			}
#endif // DRAW_CONTOUR

			Scalar QUADRILATERAL_COLOR(0, 0, 255);
			Scalar colour;
			// Convert to a usable rectangle
			Rect r = boundingRect(approx_contour);

#if DRAW_RECTANGLE
			rectangle(imgFiltered, r, Scalar(0, 0, 255), 4);
			continue;
#endif //DRAW_RECTANGLE

			// 4 Corners -> We color them
			if (approx_contour.size() == 4) {
				colour = QUADRILATERAL_COLOR;
			}
			else {
				continue;
			}

			// --- Filter tiny ones --- If the found contour is too small (20 -> pixels, frame.cols - 10 to prevent extreme big contours)
			if (r.height < 20 || r.width < 20 || r.width > imgFiltered.cols - 10 || r.height > imgFiltered.rows - 10) {
				continue;
			}

			// -> Cleaning done

			// 1 -> 1 contour, we have a closed contour, true -> closed, 4 -> thickness
			polylines(imgFiltered, approx_contour, true, colour, THICKNESS_VALUE);

			// Added in sheet 4 Ex10 - Start *******************************************************************

			// Direction vector (x0,y0) and contained point (x1,y1) -> For each line -> 4x4 = 16
			float lineParams[16];
			// lineParams is shared, CV_32F -> Same data type like lineParams
			Mat lineParamsMat(Size(4, 4), CV_32F, lineParams);

			// Added in sheet 4 Ex10 - End *******************************************************************

			// -----------------------------

			// --- Process Corners ---

			for (size_t i = 0; i < approx_contour.size(); ++i) {
				// Render the corners, 3 -> Radius, -1 filled circle
				circle(imgFiltered, approx_contour[i], 3, CV_RGB(0, 255, 0), -1);

				// Euclidic distance, 7 -> parts, both directions dx and dy
				double dx = ((double)approx_contour[(i + 1) % 4].x - (double)approx_contour[i].x) / 7.0;
				double dy = ((double)approx_contour[(i + 1) % 4].y - (double)approx_contour[i].y) / 7.0;

				// Added in Sheet 3 - Start *****************************************************************

				MyStrip strip;

				// A simple array of unsigned char cv::Mat
				Mat imagePixelStripe = calculate_Stripe(dx, dy, strip);

				// Added in Sheet 3 - End *******************************************************************

				// Added in sheet 4 - Start *****************************************************************

				// Array for edge point centers
				Point2f edgePointCenters[6];

				// Added in sheet 4 - End *****************************************************************

				// First point already rendered, now the other 6 points
				for (int j = 1; j < 7; ++j) {
					// Position calculation
					double px = (double)approx_contour[i].x + (double)j * dx;
					double py = (double)approx_contour[i].y + (double)j * dy;

					Point p;
					p.x = (int)px;
					p.y = (int)py;
					circle(imgFiltered, p, 2, CV_RGB(0, 0, 255), -1);

					//------------------------------------------- EX 3 ---------------------------------------------------------

					// Columns: Loop over 3 pixels
					for (int m = -1; m <= 1; ++m) {
						// Rows: From bottom to top of the stripe, e.g. -3 to 3
						for (int n = strip.nStart; n <= strip.nStop; ++n) {
							Point2f subPixel;

							// m -> going over the 3 pixel thickness of the stripe, n -> over the length of the stripe, direction comes from the orthogonal vector in st
							// Going from bottom to top and defining the pixel coordinate for each pixel belonging to the stripe
							subPixel.x = (double)p.x + ((double)m * strip.stripeVecX.x) + ((double)n * strip.stripeVecY.x);
							subPixel.y = (double)p.y + ((double)m * strip.stripeVecX.y) + ((double)n * strip.stripeVecY.y);

							Point p2;
							p2.x = (int)subPixel.x;
							p2.y = (int)subPixel.y;

							// The one (purple color) which is shown in the stripe window
							if (isFirstStripe)
								circle(imgFiltered, p2, 1, CV_RGB(255, 0, 255), -1);
							else
								circle(imgFiltered, p2, 1, CV_RGB(0, 255, 255), -1);

							// Combined Intensity of the subpixel
							int pixelIntensity = subpixSampleSafe(grayScale, subPixel);

							// Converte from index to pixel coordinate
							// m (Column, real) -> -1,0,1 but we need to map to 0,1,2 -> add 1 to 0..2
							int w = m + 1;

							// n (Row, real) -> add stripeLenght >> 1 to shift to 0..stripeLength
							// n=0 -> -length/2, n=length/2 -> 0 ........ + length/2
							int h = n + (strip.stripeLength >> 1);

							// Set pointer to correct position and safe subpixel intensity
							imagePixelStripe.at<uchar>(h, w) = (uchar)pixelIntensity;

							// Added in Sheet 3 - Ex7 (a) End *****************************************************************
						}
					}

					// Added in Sheet 3 - Ex7 (b) Start *****************************************************************
					// Use sobel operator on stripe

					// ( -1 , -2, -1 )
					// (  0 ,  0,  0 )
					// (  1 ,  2,  1 )

					// The first and last row must be excluded from the sobel calculation because they have no top or bottom neighbors
					vector<double> sobelValues(strip.stripeLength - 2.);

					// To use the kernel we start with the second row (n) and stop before the last one
					for (int n = 1; n < (strip.stripeLength - 1); n++) {
						// Take the intensity value from the stripe 
						unsigned char* stripePtr = &(imagePixelStripe.at<uchar>(n - 1, 0));

						// Calculation of the gradient with the sobel for the first row
						double r1 = -stripePtr[0] - 2. * stripePtr[1] - stripePtr[2];

						// r2 -> Is equal to 0 because of sobel

						// Go two lines for the third line of the sobel, step = size of the data type, here uchar
						stripePtr += 2 * imagePixelStripe.step;

						// Calculation of the gradient with the sobel for the third row
						double r3 = stripePtr[0] + 2. * stripePtr[1] + stripePtr[2];

						// Writing the result into our sobel value vector
						unsigned int ti = n - 1;
						sobelValues[ti] = r1 + r3;
					}

					double maxIntensity = -1;
					int maxIntensityIndex = 0;

					// Finding the max value
					for (int n = 0; n < strip.stripeLength - 2; ++n) {
						if (sobelValues[n] > maxIntensity) {
							maxIntensity = sobelValues[n];
							maxIntensityIndex = n;
						}
					}

					// Added in Sheet 3 - Ex7 (b) End *****************************************************************

					// Added in Sheet 3 - Ex7 (d) Start *****************************************************************

					// f(x) slide 7 -> y0 .. y1 .. y2
					double y0, y1, y2;

					// Point before and after
					unsigned int max1 = maxIntensityIndex - 1, max2 = maxIntensityIndex + 1;

					// If the index is at the border we are out of the stripe, then we will take 0
					y0 = (maxIntensityIndex <= 0) ? 0 : sobelValues[max1];
					y1 = sobelValues[maxIntensityIndex];
					// If we are going out of the array of the sobel values
					y2 = (maxIntensityIndex >= strip.stripeLength - 3) ? 0 : sobelValues[max2];

					// Formula for calculating the x-coordinate of the vertex of a parabola, given 3 points with equal distances 
					// (xv means the x value of the vertex, d the distance between the points): 
					// xv = x1 + (d / 2) * (y2 - y0)/(2*y1 - y0 - y2)

					// Equation system
					// d = 1 because of the normalization and x1 will be added later
					double pos = (y2 - y0) / (4 * y1 - 2 * y0 - 2 * y2);

					// What happens when there is no solution? -> /0 or Number = other Number
					// If the found pos is not a number -> there is no solution
					if (isnan(pos)) {
						continue;
					}
					// Check if there is a solution to the calculation, cool trick
					/*if (pos != pos) {
						// Value is not a number -> NAN
						continue;
					}*/

					// Exact point with subpixel accuracy
					Point2d edgeCenter;

					// Back to Index positioning, Where is the edge (max gradient) in the picture?
					int maxIndexShift = maxIntensityIndex - (strip.stripeLength >> 1);

					// Shift the original edgepoint accordingly -> Is the pixel point at the top or bottom?
					edgeCenter.x = (double)p.x + (((double)maxIndexShift + pos) * strip.stripeVecY.x);
					edgeCenter.y = (double)p.y + (((double)maxIndexShift + pos) * strip.stripeVecY.y);

					// Highlight the subpixel with blue color
					circle(imgFiltered, edgeCenter, 2, CV_RGB(0, 0, 255), -1);

					// Added in Sheet 3 - Ex7 (d) End *****************************************************************

					// Added in sheet 4 -  Start *****************************************************************

					edgePointCenters[j - 1].x = edgeCenter.x;
					edgePointCenters[j - 1].y = edgeCenter.y;

					// Added in sheet 4 - End *****************************************************************

					// Added in Sheet 3 - Ex7 (c) Start *****************************************************************

					// Draw the stripe in the image
					if (isFirstStripe) {
						Mat iplTmp;
						// The intensity differences on the stripe
						resize(imagePixelStripe, iplTmp, Size(100, 300));

						// imshow(stripWindow, iplTmp);
						isFirstStripe = false;
					}

					// Added in Sheet 3 - Ex7 (c) End *****************************************************************
				}

				// Added in sheet 4 Ex9(a) - Start * ****************************************************************

				// We now have the array of exact edge centers stored in "points", every row has two values -> 2 channels!
				Mat highIntensityPoints(Size(1, 6), CV_32FC2, edgePointCenters);

				// fitLine stores the calculated line in lineParams per column in the following way:
				// vec.x, vec.y, point.x, point.y
				// Norm 2, 0 and 0.01 -> Optimal parameters
				// i -> Edge points
				fitLine(highIntensityPoints, lineParamsMat.col(i), CV_DIST_L2, 0, 0.01, 0.01);
				// We need two points to draw the line
				Point p1;
				// We have to jump through the 4x4 matrix, meaning the next value for the wanted line is in the next row -> +4
				// d = -50 is the scalar -> Length of the line, g: Point + d*Vector
				// p1<----Middle---->p2
				//   <-----100----->
				p1.x = (int)lineParams[8 + i] - (int)(50.0 * lineParams[i]);
				p1.y = (int)lineParams[12 + i] - (int)(50.0 * lineParams[4 + i]);

				Point p2;
				p2.x = (int)lineParams[8 + i] + (int)(50.0 * lineParams[i]);
				p2.y = (int)lineParams[12 + i] + (int)(50.0 * lineParams[4 + i]);

				// Draw line
				line(imgFiltered, p1, p2, CV_RGB(0, 255, 255), 1, 8, 0);

				// Added in sheet 4 Ex9 (a)- End *******************************************************************
			}

			// Added in sheet 4 Ex9 (b)- Start *****************************************************************

			// So far we stored the exact line parameters and show the lines in the image now we have to calculate the exact corners
			Point2f corners[4];

			// Calculate the intersection points of both lines
			for (int i = 0; i < 4; ++i) {
				// Go through the corners of the rectangle, 3 -> 0
				int j = (i + 1) % 4;

				double x0, x1, y0, y1, u0, u1, v0, v1;

				// We have to jump through the 4x4 matrix, meaning the next value for the wanted line is in the next row -> +4
				// g: Point + d*Vector
				// g1 = (x0,y0) + scalar0*(u0,v0) == g2 = (x1,y1) + scalar1*(u1,v1)
				x0 = lineParams[i + 8]; y0 = lineParams[i + 12];
				x1 = lineParams[j + 8]; y1 = lineParams[j + 12];

				// Direction vector
				u0 = lineParams[i]; v0 = lineParams[i + 4];
				u1 = lineParams[j]; v1 = lineParams[j + 4];

				// (x|y) = p + s * vec --> Vector Equation

				// (x|y) = p + (Ds / D) * vec

				// p0.x = x0; p0.y = y0; vec0.x= u0; vec0.y=v0;
					// p0 + s0 * vec0 = p1 + s1 * vec1
					// p0-p1 = vec(-vec0 vec1) * vec(s0 s1)	

					// s0 = Ds0 / D (see cramer's rule)
					// s1 = Ds1 / D (see cramer's rule)   
					// Ds0 = -(x0-x1)v1 + (y0-y1)u1 --> You need to just calculate one, here Ds0

				// (x|y) = (p * D / D) + (Ds * vec / D)
				// (x|y) = (p * D + Ds * vec) / D

					// x0 * D + Ds0 * u0 / D    or   x1 * D + Ds1 * u1 / D     --> a / D
					// y0 * D + Ds0 * v0 / D    or   y1 * D + Ds1 * v1 / D     --> b / D						   		

				// (x|y) = a / c;

				// Cramer's rule
				// 2 unknown a,b -> Equation system
				double a = x1 * u0 * v1 - y1 * u0 * u1 - x0 * u1 * v0 + y0 * u0 * u1;
				double b = -x0 * v0 * v1 + y0 * u0 * v1 + x1 * v0 * v1 - y1 * v0 * u1;

				// Calculate the cross product to check if both direction vectors are parallel -> = 0
				// c -> Determinant = 0 -> linear dependent -> the direction vectors are parallel -> No division with 0
				double c = v1 * u0 - v0 * u1;
				if (fabs(c) < 0.001) {
					std::cout << "lines parallel" << std::endl;
					continue;
				}

				// We have checked for parallelism of the direction vectors
				// -> Cramer's rule, now divide through the main determinant
				a /= c;
				b /= c;

				// Exact corner
				corners[i].x = a;
				corners[i].y = b;

				// Added in sheet 4 Ex9 (b)- End *******************************************************************

				Point p;
				p.x = (int)corners[i].x;
				p.y = (int)corners[i].y;

				// Added in sheet 4 Ex9 (c)- Start *****************************************************************

				circle(imgFiltered, p, 5, CV_RGB(255, 255, 0), -1);

				// Added in sheet 4 Ex9 (c)- End *******************************************************************

			} // End of the loop to extract the exact corners

			// Added in sheet 4 Ex10 (a)- Start *****************************************************************

			// Coordinates on the original marker images to go to the actual center of the first pixel -> 6x6
			Point2f targetCorners[4];
			targetCorners[0].x = -0.5; targetCorners[0].y = -0.5;
			targetCorners[1].x = 5.5; targetCorners[1].y = -0.5;
			targetCorners[2].x = 5.5; targetCorners[2].y = 5.5;
			targetCorners[3].x = -0.5; targetCorners[3].y = 5.5;

			// Create and calculate the matrix of perspective transform -> non-affine -> parallel stays not parallel
			// Homography is a matrix to describe the transformation from an image region to the 2D projected image
			Mat homographyMatrix(Size(3, 3), CV_32FC1);
			// Corner which we calculated and our target Mat, find the transformation
			homographyMatrix = getPerspectiveTransform(corners, targetCorners);

			// Added in sheet 4 Ex10 (a)- End *******************************************************************

			// Added in sheet 4 Ex10 (b)- Start *****************************************************************

			// Create image for the marker
			Mat imageMarker(Size(6, 6), CV_8UC1);

			// Change the perspective in the marker image using the previously calculated Homography Matrix
			// In the Homography Matrix there is also the position in the image saved
			warpPerspective(grayScale, imageMarker, homographyMatrix, Size(6, 6));

			// Added in sheet 4 Ex10 (b)- End *******************************************************************

			// Now we have a B/W image of a supposed Marker
			threshold(imageMarker, imageMarker, bw_thresh, 255, CV_THRESH_BINARY);

			// Added in sheet 4 Ex10 (c)- Start *****************************************************************

			int code = 0;
			for (int i = 0; i < 6; ++i) {
				// Check if border is black
				int pixel1 = imageMarker.at<uchar>(0, i); //top
				int pixel2 = imageMarker.at<uchar>(5, i); //bottom
				int pixel3 = imageMarker.at<uchar>(i, 0); //left
				int pixel4 = imageMarker.at<uchar>(i, 5); //right

				// 0 -> black
				if ((pixel1 > 0) || (pixel2 > 0) || (pixel3 > 0) || (pixel4 > 0)) {
					code = -1;
					break;
				}
			}

			if (code < 0) {
				continue;
			}

			// Copy the BW values into cP -> codePixel on the marker 4x4 (inner part of the marker, no black border)
			int cP[4][4];
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < 4; ++j) {
					// +1 -> no borders!
					cP[i][j] = imageMarker.at<uchar>(i + 1, j + 1);
					// If black then 1 else 0
					cP[i][j] = (cP[i][j] == 0) ? 1 : 0;
				}
			}

			// Save the ID of the marker, for each side
			int codes[4];
			codes[0] = codes[1] = codes[2] = codes[3] = 0;

			// Calculate the code from all sides in one loop
			for (int i = 0; i < 16; i++) {
				// /4 to go through the rows
				int row = i >> 2;
				int col = i % 4;

				// Multiplied by 2 to check for black values -> 0*2 = 0
				codes[0] <<= 1;
				codes[0] |= cP[row][col]; // 0°

				// 4x4 structure -> Each column represents one side 
				codes[1] <<= 1;
				codes[1] |= cP[3 - col][row]; // 90°

				codes[2] <<= 1;
				codes[2] |= cP[3 - row][3 - col]; // 180°

				codes[3] <<= 1;
				codes[3] |= cP[col][3 - row]; // 270°

				/*cout << "iteration: " << dec << i << endl;
				cout << "Code 0: " << hex << codes[0] << endl;
				cout << "Code 1: " << hex << codes[1] << endl;
				cout << "Code 2: " << hex << codes[2] << endl;
				cout << "Code 3: " << hex << codes[3] << endl;*/
			}

			// Account for symmetry -> One side complete white or black
			if ((codes[0] == 0) || (codes[0] == 0xffff)) {
				continue;
			}

			// Added in Exercise 5 - Start *****************************************************************

			int angle = 0;

			// Added in Exercise 5 - End *******************************************************************

			// Search for the smallest marker ID
			code = codes[0];
			for (int i = 1; i < 4; ++i) {
				if (codes[i] < code) {
					code = codes[i];

					// Added in Exercise 5 - Start *****************************************************************

					angle = i;

					// Added in Exercise 5 - End *******************************************************************
				}
			}

			// Print ID
			//printf("Found: %04x\n", code);

			identifiers.push_back(code);

			// Show the first detected marker in the image
			if (isFirstMarker) {
				imshow(kWinName4, imageMarker);
				isFirstMarker = false;
			}

			// Added in sheet 4 Ex10 (c)- End *****************************************************************

			// Added in sheet 4 - End *******************************************************************

			// Added in Exercise 5 - Start *****************************************************************

			// Correct the order of the corners, if 0 -> already have the 0 degree position
			/* TASK: How to make sure, that we always have the correct order? */
			if (angle != 0) {
				Point2f corrected_corners[4];
				// Smallest id represents the x-axis, we put the values in the corrected_corners array
				for (int i = 0; i < 4; i++)	corrected_corners[(i + angle) % 4] = corners[i];
				// We put the values back in the array in the sorted order
				for (int i = 0; i < 4; i++)	corners[i] = corrected_corners[i];
			}

			// Normally we should do a camera calibration to get the camera paramters such as focal length
			// Two ways: Inner parameters, e.g. focal length (intrinsic parameters); camera with 6 dof (R|T) (extrinsic parameters)
			// Transfer screen coords to camera coords -> To get to the principal point
			for (int i = 0; i < 4; i++) {
				// TODO Here you have to use your own camera resolution (x) * 0.5
				corners[i].x -= capResWidth * 0.5;
				// -(corners.y) -> is needed because y is inverted
				// TODO Here you have to use your own camera resolution (y) * 0.5
				corners[i].y = -corners[i].y + capResHeight * 0.5;
			}

			// 4x4 -> Rotation | Translation
			//        0  0  0  | 1 -> (Homogene coordinates to combine rotation, translation and scaling)
			float resultMatrix[16];
			// TODO Marker size in meters!
			estimateSquarePose(resultMatrix, (Point2f*)corners, (double) 0.04346);

			// This part is only for printing
			/*
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < 4; ++j) {
					cout << setw(6); // Total 6
					cout << setprecision(4); // Numbers of decimal places = 4 (of the 6)
					// Again we are going through a matrix which is saved as an array
					cout << resultMatrix[4 * i + j] << " ";
				}
				cout << "\n";
			}
			cout << "\n";

			*/
			float x, y, z;
			// Translation values in the transformation matrix to calculate the distance between the marker and the camera
			x = resultMatrix[3];
			y = resultMatrix[7];
			z = resultMatrix[11];
			// Euclidian distance
			//TODO check distance, will likely be too short
			//cout << "length: " << sqrt(x * x + y * y + z * z) << "\n";
			//cout << "\n";

			// Added in Exercise 5 - End *****************************************************************

			// -----------------------------

			// -----------------------------
		}

		// ---Play sound if markers are present ---
		for (Instrument& instr : instruments) {
			instr.toggleSound(identifiers);
		}


		imshow(contoursWindow, imgFiltered);
		isFirstStripe = true;

		int pressedKey = waitKey(10);
		if (pressedKey == 27) {
			//end program
			break;
		}
		/*
		else if (pressedKey == 32) {
			//switch synth volume on or off
			if (Mix_Volume(channel_drums, 0) == 0) Mix_Volume(channel_drums, 128);
		}
		*/

	}

	destroyWindow(contoursWindow);

	for (Instrument& instr : instruments) {
		instr.freeChunk();
	}

	Mix_CloseAudio();
	SDL_Quit();

	return 0;
}

#include<opencv2/opencv.hpp>
#include<iostream>
#include<array>
#include "SDL2/SDL.h"
#include "SDL2/SDL_mixer.h"


using namespace cv;
using namespace std;

#define TEST 0
#define EX5 1
#define EX5_RAW 0
#define DRAW_CONTOUR 0
#define DRAW_RECTANGLE 0

#define THICKNESS_VALUE 4

#define SYNTH_PATH "" //TODO

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

const std::string kWinName4 = "Exercise 4 - Marker";

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
	int audio_rate = 22050;
	Uint16 audio_format = AUDIO_S16SYS;
	int audio_channels = 4;
	int audio_buffers = 4096;

	Mat frame;
	VideoCapture cap(0);

	vector<array<int, 4> > identifiers;

	const string streamWindow = "Stream";

	if (!cap.isOpened()) {
		cout << "No webcam, using video file" << endl;
		cap.open("MarkerMovie.MP4");
		if (cap.isOpened() == false) {
			cout << "No video!" << endl;
			exit(0);
			return -1;
		}
	}

	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) != 0) {
		fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
		return 1;
	}
	if (Mix_OpenAudio(audio_rate, audio_format, audio_channels, audio_buffers) != 0) {
		fprintf(stderr, "Unable to initialize audio: %s\n", Mix_GetError());
		exit(1);
	}

	Mix_Chunk* sound_synth = Mix_LoadWAV("/Users/sam/Desktop/arwork/ex01/synth130.wav");
	if (sound_synth == NULL) {
		fprintf(stderr, "Unable to load WAV file: %s\n", Mix_GetError());
	}
	int channel_synth = Mix_PlayChannel(-1, sound_synth, 0);
	if (channel_synth == -1) {
		fprintf(stderr, "Unable to play WAV file: %s\n", Mix_GetError());
	}
	Mix_Volume(channel_synth, 0);


	Mix_Chunk* sound_drums = Mix_LoadWAV("/Users/sam/Desktop/arwork/ex01/drums130.wav");
	if (sound_drums == NULL) {
		fprintf(stderr, "Unable to load WAV file: %s\n", Mix_GetError());
	}
	int channel_drums = Mix_PlayChannel(-1, sound_drums, 0);
	if (channel_drums == -1) {
		fprintf(stderr, "Unable to play WAV file: %s\n", Mix_GetError());
	}
	Mix_Volume(channel_drums, 0);


	// Added in Sheet 3 - Start *****************************************************************

	bool isFirstStripe = true;

	// Added in Sheet 3 - End *******************************************************************

	// Added in sheet 4 - Start *****************************************************************

	bool isFirstMarker = true;

	// Added in sheet 4 - End *******************************************************************

	const string contoursWindow = "Contours";
	const string UI = "Threshold";
	namedWindow(contoursWindow, CV_WINDOW_AUTOSIZE);
	//namedWindow(stripWindow, CV_WINDOW_AUTOSIZE);
	const string markerWindow = "Marker";
	namedWindow(markerWindow, CV_WINDOW_NORMAL);


	int slider_value = 100;
	createTrackbar(UI, contoursWindow, &slider_value, 255, on_trackbar, &slider_value);

	// Added in sheet 4 - Start *****************************************************************

	//int bw_sileder_value = bw_thresh;

	//createTrackbar("BW Threshold", contoursWindow, &bw_sileder_value, 255, bw_trackbarHandler, &bw_sileder_value);

	namedWindow(kWinName4, CV_WINDOW_NORMAL);

	resizeWindow(kWinName4, 120, 120);

	// Added in sheet 4 - End *******************************************************************

	Mat imgFiltered;

	while (cap.read(frame)) {

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

						imshow(stripWindow, iplTmp);
						isFirstStripe = false;
					}

					// Added in Sheet 3 - Ex7 (c) End *****************************************************************
				}

				//Sheet 4 - Ex9 (a) Start

				//We now have the array of exact edge centers stored in "points", every row has 2 values x and y
				Mat highIntensityPoints(Size(1, 6), CV_32FC2, edgePointCenters);

				//fitLine stores the calculated line in lineParams per column in the following way:
				// vec.x, vec.y, point.x, point.y
				// Norm 2, 0 and 0.01 -> Optimal parameters
				// i -> Edge points
				fitLine(highIntensityPoints, lineParamsMat.col(i), CV_DIST_L2, 0, 0.01, 0.01);
				
				//need 2 points to draw the line
				Point p1;

				//Jump through the 4x4 matrix
				//d = -50 is the scalar -> Length of the line, g: Point + d*Vector
				//p1 <-----Middle----->p2
				//<-----100------>
				p1.x = (int)lineParams[8 + i] - (int)(50.0 * lineParams[i]);
				p1.y = (int)lineParams[12 + i] - (int)(50.0 * lineParams[4 + i]);

				Point p2;
				p2.x = (int)lineParams[8 + i] + (int)(50.0 * lineParams[i]);
				p2.y = (int)lineParams[12 + i] + (int)(50.0 * lineParams[4 + i]);

				//Draw line
				line(imgFiltered, p1, p2, CV_RGB(0, 255, 255), 1, 8, 0);

				//Sheet 4 - Ex9 (a) End

			}

			// Sheet 4 - Ex9 (b) Start -----------------------------
			//calculate corners from edge lines
			Point2f corners[4];

			//calculate intersection points of both lines
			for (int i = 0; i < 4; i++) {
				//Go through the corners of the rectangle, 3 -> 0
				int j = (i + 1) % 4;

				double x0, x1, y0, y1, u0, u1, v0, v1;

				//Jump through 4x4 again
				//g: Point + d*Vector
				//g1 = (x0,y0) + scalar0*(u0,v0) == g2 = (x1, y1) + scalar1*(u1,v1)
				x0 = lineParams[i + 8]; y0 = lineParams[i + 12];
				x1 = lineParams[j + 8]; y1 = lineParams[j + 12];

				//Direction vector
				u0 = lineParams[i]; v0 = lineParams[i + 4];
				u1 = lineParams[j]; v1 = lineParams[j + 4];

				//Cramer's rule, see tutorial
				//2 unknown a, b -> Equation system
				double a = x1 * u0 * v1 - y1 * u0 * u1 - x0 * u1 * v0 + y0 * u0 * u1;
				double b = -x0 * v0 * v1 + y0 * u0 * v1 + x1 * v0 * v1 - y1 * v0 * u1;

				//Calculate the cross product to check if both direction vectors are parallel -> = 0
				// c -> Determinant = 0 -> linear dependent -> direction vectors are parallel -> No division
				double c = v1 * u0 - v0 * u1;
				if (fabs(c) < 0.001) {
					cout << "lines parallel" << endl;
					continue;
				}

				//vectors not parallel -> Cramer's rule, divide through the main determinant
				a /= c;
				b /= c;

				//Exact corner
				corners[i].x = a;
				corners[i].y = b;

				// Sheet 4 - Ex9 (b) End -------------------------------------------------

				// Sheet 4 - Ex9 (c) Start -------------------------------------------------

				Point p;
				p.x = (int)corners[i].x;
				p.y = (int)corners[i].y;

				//draw corner point
				circle(imgFiltered, p, 5, CV_RGB(255, 255, 0), -1);

				// Sheet 4 - Ex9 (c) End -------------------------------------------------

			} //exact corners extracted

			//Sheet 4 - Ex10 (a) Start --------------------------------

			Point2f targetCorners[4];
			targetCorners[0].x = -0.5; targetCorners[0].y = -0.5;
			targetCorners[1].x = 5.5; targetCorners[1].y = -0.5;
			targetCorners[2].x = 5.5; targetCorners[2].y = 5.5;
			targetCorners[3].x = -0.5; targetCorners[3].y = 5.5;

			Mat warpMat = getPerspectiveTransform(corners, targetCorners);

			//Sheet 4 - Ex10 (a) End --------------------

			//Sheet 4 - Ex10 (b) Start --------------------
			Mat warped(Size(6, 6), CV_32FC2, Scalar::all(0));
			warpPerspective(grayScale, warped, warpMat, Size(6, 6),1,1, Scalar::all(1));

			//Sheet 4 - Ex10 (b) End --------------------

			//Sheet 4- Ex10 (c) Start ---------------
			//discard markers w/o black border and generate identifier
			bool discard = false;

			vector<array<int, 4> > rotationIds;
			array<int, 4> inverseIdentifier; //looking for white tiles, later we invert the result
			for (int i = 0; i < 4; i++) inverseIdentifier[i] = 0;

			double angle = 90;
			for (int i = 0; i < 4; i++) {
				//rotate by 90� four times to find all possible ids of this marker
				Point2f center((warped.cols - 1) / 2.0, (warped.rows - 1) / 2.0);
				// get the center coordinates of the image to create the 2D rotation matrix
				
				// using getRotationMatrix2D() to get the rotation matrix
				Mat rotation_matix = getRotationMatrix2D(center, angle, 1.0);
				//rotate marker image using mark affine
				warpAffine(warped, warped, rotation_matix, warped.size());


				vector<Point> nonZeroLocations;
				findNonZero(warped, nonZeroLocations);

				if (nonZeroLocations.size() == 0) continue; //discard all black

				for (int i = 0; i < nonZeroLocations.size(); i++) {
					int x = (int)nonZeroLocations[i].x;
					int y = (int)nonZeroLocations[i].y;
					if (y == 0 || y == 5 || x == 0 || x == 5) {
						discard = true;
						break;
					}
					//interior point, used for identifier
					inverseIdentifier[y - 1] += pow(2, (4 - x));
				}
				if (discard) continue;

				//invert the identifier and add to the ones we found
				for (int i = 0; i < inverseIdentifier.size(); i++) {
					inverseIdentifier[i] = 15 - inverseIdentifier[i];
				}
				rotationIds.push_back(inverseIdentifier);
			}

			//canonical id of marker will be minimal
			array<int, 4> identifier = getCanonicalId(rotationIds);

			identifiers.push_back(identifier);
			
			//slow way to remove duplicate ids found
			sort(identifiers.begin(), identifiers.end());
			auto last = unique(identifiers.begin(), identifiers.end());
			identifiers.erase(last, identifiers.end());

			//Sheet 4 - Ex10 (c) End -------------------

			imshow(markerWindow, warped);
		}

		//0690, 0272
		//check if synth marker is present
		array<int, 4> synth_marker = {-6, 9, 9, -6 };
		array<int, 4> drums_marker = { -2, -3, 12, 4};
		if (find(identifiers.begin(), identifiers.end(), synth_marker) != identifiers.end()) {
			Mix_Volume(channel_synth, 128);
		}
		else {
			Mix_Volume(channel_synth, 0);
		}
		if (find(identifiers.begin(), identifiers.end(), drums_marker) != identifiers.end()) {
			Mix_Volume(channel_drums, 128);
		}
		else {
			Mix_Volume(channel_drums, 0);
		}

		//restart loop
		if (Mix_Playing(channel_synth) != 1) Mix_PlayChannel(channel_synth, sound_synth, 0);
		if (Mix_Playing(channel_drums) != 1) Mix_PlayChannel(channel_drums, sound_drums, 0);


		imshow(contoursWindow, imgFiltered);
		isFirstStripe = true;

		int pressedKey = waitKey(10);
		if (pressedKey == 27) {
			//end program
			for (array<int, 4> id : identifiers) {
				cout << "identifier: (";
				for (int i = 0; i < 4; i++) {
					cout << id[i] << ",";
				}
				cout << ")" << endl;
			}
			break;
		}
		else if (pressedKey == 32) {
			//switch synth volume on or off
			if (Mix_Volume(channel_drums, 0) == 0) Mix_Volume(channel_drums, 128);
		}

	}

	destroyWindow(contoursWindow);
	Mix_FreeChunk(sound_synth);
	Mix_FreeChunk(sound_drums);

	Mix_CloseAudio();
	SDL_Quit();

	return 0;
}


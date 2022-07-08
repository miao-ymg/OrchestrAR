#ifndef MARKER_TRACKER_H
#define MARKER_TRACKER_H

#include <opencv/cv.h>
#include "Instrument.h"

class CvMemStorage;

const std::string kWinName1 = "Exercise 8 - Original Image";
const std::string kWinName2 = "Exercise 8 - Converted Image";
const std::string kWinName3 = "Exercise 8 - Stripe Image";
const std::string kWinName4 = "Exercise 8 - Marker";

class MarkerTracker{
public:
	MarkerTracker(double kMarkerSize_) 
		:	thresh(100),
		bw_thresh(100),
		kMarkerSize(kMarkerSize_)
	{
		init();
	}
	MarkerTracker(double kMarkerSize_, int thresh_, int bw_thresh_) 
		:	thresh(thresh_),
		bw_thresh(bw_thresh_),
		kMarkerSize(kMarkerSize_)
	{
		init();
	}
	~MarkerTracker(){
		cleanup();
	}
	void findMarker( cv::Mat &img_bgr, float resultMatrix[16], unordered_map<int, Instrument> &instruments, vector<int> &identifers);
protected:
	void init( );
	void cleanup( );

	//camera settings
	const double kMarkerSize; // Marker size [m]

	cv::Mat img_gray;
	cv::Mat img_mono;

	int thresh; // Threshold (gray to mono)
	int bw_thresh; // threshold for (gray maker to ID image)

	CvMemStorage* memStorage;
};

#endif // MARKER_TRACKER_H

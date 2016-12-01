#include "plotcolorspace.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/types_c.h>
#include "colorconvert.h"

int main(void) {
    size_t const steps = 100;
    size_t min_L = 0; size_t max_L = 100;

    double const min_a = -170;
    double const max_a = 100;
    double const min_b = -100;
    double const max_b = 150;

    size_t const width = 120;
    size_t const height = 100;

    for (size_t frameNum = 0; frameNum < steps; frameNum++) {
        cv::Mat_<cv::Vec3b> frame(height, width, cv::Vec3b(100,100,100));
        double const L = min_L + static_cast<double>(frameNum) / (steps-1) * (max_L - min_L);

        for (int ii = 0; ii < frame.rows; ++ii) {
            double const b = min_b + static_cast<double>(ii) / (frame.rows-1) * (max_b - min_b);
            cv::Vec3b* frameCol = frame.ptr<cv::Vec3b>(ii);
            for (int jj = 0; jj < frame.cols; ++jj) {
                double const a = min_a + static_cast<double>(jj) / (frame.cols-1) * (max_a - min_a);
                ColorConvert::rgb2DIN()
            }
        }
    }
}

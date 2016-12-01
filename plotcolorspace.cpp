#include "plotcolorspace.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/types_c.h>
#include "colorconvert.h"

template<class T>
/**
 * @brief toStringLZ
 * @param val
 * @return
 */
std::string toStringLZ(const T val, const size_t minlength) {
    std::string res = std::to_string(val);
    if (res.size() < minlength) {
        res = std::string(minlength - res.size(), '0') + res;
    }
    return res;
}


int main(void) {
    size_t const steps = 100;
    size_t min_L = 0; size_t max_L = 100;

    double const min_a = -170;
    double const max_a = 100;
    double const min_b = -100;
    double const max_b = 150;

    size_t const width = 1000;
    size_t const height = 1000;
    double min_threshold = .01;
    double max_threshold = 254.99;

    double min_success_a = max_a;
    double max_success_a = min_a;
    double min_success_b = max_b;
    double max_success_b = min_b;


#pragma omp parallel for
    for (size_t frameNum = 0; frameNum < steps; frameNum++) {
        cv::Mat_<cv::Vec3b> frame(height, width, cv::Vec3b(100,100,100));
        double const L = min_L + static_cast<double>(frameNum) / (steps-1) * (max_L - min_L);
        size_t counter = 0;
        size_t success_counter = 0;
        for (int ii = 0; ii < frame.rows; ++ii) {
            double const b = min_b + static_cast<double>(ii) / (frame.rows-1) * (max_b - min_b);
            cv::Vec3b* frameCol = frame.ptr<cv::Vec3b>(ii);
            for (int jj = 0; jj < frame.cols; ++jj) {
                counter++;
                double const a = min_a + static_cast<double>(jj) / (frame.cols-1) * (max_a - min_a);
                cv::Vec3d vec(L,a,b);
                auto RGB = ColorConvert::DIN2rgb(vec);

                if (
                        RGB[0] >= min_threshold && RGB[0] <= max_threshold &&
                        RGB[1] >= min_threshold && RGB[1] <= max_threshold &&
                        RGB[2] >= min_threshold && RGB[2] <= max_threshold) {
#pragma omp critical
                    {
                        frameCol[jj][0] = cv::saturate_cast<uchar>(RGB[2]);
                        frameCol[jj][1] = cv::saturate_cast<uchar>(RGB[1]);
                        frameCol[jj][2] = cv::saturate_cast<uchar>(RGB[0]);
                        success_counter++;
                        min_success_a = std::min(a, min_success_a);
                        min_success_b = std::min(b, min_success_b);
                        max_success_a = std::max(a, max_success_a);
                        max_success_b = std::max(b, max_success_b);
                    }
                }
            }
        }
        std::cout << "In frame #" << frameNum << " we have " << 100 * static_cast<double>(success_counter) / static_cast<double>(counter) << "% success" << std::endl;
    }

#pragma omp parallel for
    for (size_t frameNum = 0; frameNum < steps; frameNum++) {
        cv::Mat_<cv::Vec3b> frame(height, width, cv::Vec3b(100,100,100));
        double const L = min_L + static_cast<double>(frameNum) / (steps-1) * (max_L - min_L);
        size_t counter = 0;
        size_t success_counter = 0;
        for (int ii = 0; ii < frame.rows; ++ii) {
            double const b = min_success_b + static_cast<double>(ii) / (frame.rows-1) * (max_success_b - min_success_b);
            cv::Vec3b* frameCol = frame.ptr<cv::Vec3b>(ii);
            for (int jj = 0; jj < frame.cols; ++jj) {
                counter++;
                double const a = min_success_a + static_cast<double>(jj) / (frame.cols-1) * (max_success_a - min_success_a);
                cv::Vec3d vec(L,a,b);
                auto RGB = ColorConvert::DIN2rgb(vec);
                if (
                        RGB[0] >= min_threshold && RGB[0] <= max_threshold &&
                        RGB[1] >= min_threshold && RGB[1] <= max_threshold &&
                        RGB[2] >= min_threshold && RGB[2] <= max_threshold) {
                    frameCol[jj][0] = cv::saturate_cast<uchar>(RGB[2]);
                    frameCol[jj][1] = cv::saturate_cast<uchar>(RGB[1]);
                    frameCol[jj][2] = cv::saturate_cast<uchar>(RGB[0]);
                    success_counter++;
                }
            }
        }
        cv::imwrite(toStringLZ(frameNum, 4) + ".png", frame);
        std::cout << "In frame #" << frameNum << " we have " << 100 * static_cast<double>(success_counter) / static_cast<double>(counter) << "% success" << std::endl;
    }
}

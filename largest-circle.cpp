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
    size_t const width = 1000;
    size_t const height = 1000;
    double min_threshold = 0.001;
    double max_threshold = 254.999;

    double const min_color = -33.28;
    double const max_color = 36;

    double const L = 66;

    cv::Mat_<cv::Vec3b> frame(height, width, cv::Vec3b(100,100,100));
    cv::Mat_<uint8_t> success_px(height, width, uint8_t(0));
    size_t counter = 0;
    size_t success_counter = 0;
    for (int ii = 0; ii < frame.rows; ++ii) {
        double const b = min_color + static_cast<double>(ii) / (frame.rows-1) * (max_color - min_color);
        cv::Vec3b* frameCol = frame.ptr<cv::Vec3b>(ii);
        for (int jj = 0; jj < frame.cols; ++jj) {
            counter++;
            double const a = min_color + static_cast<double>(jj) / (frame.cols-1) * (max_color - min_color);
            cv::Vec3d vec(L,a,b);
            auto RGB = ColorConvert::DIN2rgb(vec);
            if (
                    RGB[0] >= min_threshold && RGB[0] <= max_threshold &&
                    RGB[1] >= min_threshold && RGB[1] <= max_threshold &&
                    RGB[2] >= min_threshold && RGB[2] <= max_threshold) {
                frameCol[jj][0] = cv::saturate_cast<uchar>(RGB[2]);
                frameCol[jj][1] = cv::saturate_cast<uchar>(RGB[1]);
                frameCol[jj][2] = cv::saturate_cast<uchar>(RGB[0]);
                success_px(ii, jj) = 255;
                success_counter++;
            }
        }
    }
    cv::cvtColor(frame, frame, cv::COLOR_RGB2BGR);
    cv::imwrite("optimum.png", frame);
    cv::imwrite("success.png", success_px);
    std::cout << "a : [" << min_color << ", " << max_color << "]" << std::endl;
    std::cout << "b : [" << min_color << ", " << max_color << "]" << std::endl;


    cv::Mat_<float> distances;
    cv::distanceTransform(success_px, distances, cv::DIST_L2, cv::DIST_MASK_PRECISE);

    double max_dist = 0;
    cv::Point2i center;
    for (int ii = 1; ii+1 < frame.rows; ++ii) {
        for (int jj = 1; jj+1 < frame.cols; ++jj) {
            if (distances(ii, jj) > max_dist) {
                max_dist = distances(ii, jj);
                center = cv::Point2i(jj, ii);
            }
        }
    }
    std::cout << "Max dist " << max_dist << " at " << center << std::endl;
    cv::circle(frame, center, max_dist, cv::Scalar(255,255,255), 1, cv::LINE_AA);
    cv::imwrite("circle.png", frame);

    double a_center = min_color + double(center.x) / (frame.rows-1) * (max_color - min_color);
    double b_center = min_color + double(center.y) / (frame.rows-1) * (max_color - min_color);
    double color_radius = double(max_dist) / (frame.rows-1) * (max_color - min_color);

    std::cout << "a/b center: " << a_center << ", " << b_center << std::endl;
    std::cout << "color radius: " << color_radius << std::endl;

    /*
    std::vector<cv::Point2i> edge;
    cv::Mat_<uint8_t> edge_px(height, width, uint8_t(0));
    for (int ii = 1; ii+1 < frame.rows; ++ii) {
        for (int jj = 1; jj+1 < frame.cols; ++jj) {
            if (!success_px(ii, jj) && (
                           success_px(ii+1, jj)
                        || success_px(ii, jj+1)
                        || success_px(ii-1, jj)
                        || success_px(ii, jj-1)
                        )
                    ) {
                edge.push_back(cv::Point2i(jj,ii));
            }
        }
    }
    */


}

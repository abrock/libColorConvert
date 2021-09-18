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
    size_t const steps = 200;
    size_t min_L = 0; size_t max_L = 100;

    double const min_color = -40;
    double const max_color = 40;

    size_t const width = 1000;
    double min_threshold = 0.001;
    double max_threshold = 254.999;

    double min_success_color = max_color;
    double max_success_color = min_color;


    std::vector<cv::Mat_<cv::Vec3b> > frames(steps);


#pragma omp parallel for ordered schedule(dynamic)
    for (size_t frameNum = 0; frameNum < steps; frameNum++) {
        frames[frameNum] = cv::Mat_<cv::Vec3b>(width, width, cv::Vec3b(100,100,100));
        cv::Mat_<uint8_t> success_px(width, width, uint8_t(0));
        double min_local_success_color = max_color;
        double max_local_success_color = min_color;
        cv::Mat_<cv::Vec3b> &frame = frames[frameNum];
        double const L = min_L + static_cast<double>(frameNum) / (steps-1) * (max_L - min_L);
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
                        RGB[0] > min_threshold && RGB[0] < max_threshold &&
                        RGB[1] > min_threshold && RGB[1] < max_threshold &&
                        RGB[2] > min_threshold && RGB[2] < max_threshold) {
#pragma omp critical
                    {
                        frameCol[jj][0] = cv::saturate_cast<uchar>(RGB[2]);
                        frameCol[jj][1] = cv::saturate_cast<uchar>(RGB[1]);
                        frameCol[jj][2] = cv::saturate_cast<uchar>(RGB[0]);
                        success_counter++;
                        success_px(ii, jj) = 255;
                        min_success_color = std::min(a, min_success_color);
                        max_success_color = std::max(a, max_success_color);
                        min_success_color = std::min(b, min_success_color);
                        max_success_color = std::max(b, max_success_color);

                        min_local_success_color = std::min(a, min_local_success_color);
                        max_local_success_color = std::max(a, max_local_success_color);
                        min_local_success_color = std::min(b, min_local_success_color);
                        max_local_success_color = std::max(b, max_local_success_color);

                    }
                }
            }
        }
        cv::cvtColor(frame, frame, cv::COLOR_RGB2BGR);
        cv::imwrite(toStringLZ(frameNum, 4) + ".png", frame);

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
        cv::circle(frame, center, max_dist, cv::Scalar(255,255,255), 1, cv::LINE_AA);
        cv::imwrite("circle-" + toStringLZ(frameNum, 4) + ".png", frame);

        double a_center = min_color + double(center.x) / (frame.rows-1) * (max_color - min_color);
        double b_center = min_color + double(center.y) / (frame.rows-1) * (max_color - min_color);
        double color_radius = double(max_dist) / (frame.rows-1) * (max_color - min_color);

#pragma omp ordered
        {
            std::cout << "In frame #" << toStringLZ(frameNum, 3) << ", L= " << L << " we have " << 100 * static_cast<double>(success_counter) / static_cast<double>(counter) << "% success, ";
            std::cout << "max dist " << max_dist << " at " << center << std::endl;
            std::cout << "color : [" << min_local_success_color << ", " << max_local_success_color << "]" << std::endl;
            std::cout << "a/b center: " << a_center << ", " << b_center << std::endl;
            std::cout << "color radius: " << color_radius << std::endl;

        }
    }

    std::cout << "Finished measuring color body dimension" << std::endl;

    std::cout << "color : [" << min_success_color << ", " << max_success_color << "]" << std::endl;

}

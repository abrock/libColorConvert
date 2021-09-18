#ifndef COLORFLOW_H
#define COLORFLOW_H

#include "colorconvert.h"
#include <opencv2/highgui.hpp>

class ColorFlow
{
private:
    ColorFlow();

    // Result from plotcolorspace:
    // In frame #136, L= 68.3417 we have 36.7637% success,
    // max dist 303.579 at [588, 478]
    // color : [-25.986, 33.4334]
    // a/b center: 7.08709, -1.72172
    // color radius: 24.3106

    static constexpr double const L = 68.3417;
    static constexpr double const a_center = 7.08709;
    static constexpr double const b_center = -1.72172;
    static constexpr double const color_radius = 24.3106;
    static constexpr size_t const num_colors = 3'600;

    std::vector<cv::Vec3f> colors;
public:
    ColorFlow& operator=(ColorFlow const&) = delete;
    ColorFlow(ColorFlow const&) = delete;

    static ColorFlow& getInstance();
    static cv::Mat_<cv::Vec3b> colorFlow(const cv::Mat_<cv::Vec2f> &flow);
    static cv::Vec3b colorFlow(const cv::Vec2f &flow);
    static cv::Vec3f colorDirection(const cv::Vec2f &dir);
    static cv::Mat_<cv::Vec3b> colorDirection(const cv::Mat_<cv::Vec2f> &flow);
};

#endif // COLORFLOW_H

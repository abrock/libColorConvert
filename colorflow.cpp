#undef NDEBUG
#include <cassert>

#include "colorflow.h"

#include <opencv2/imgproc.hpp>

ColorFlow::ColorFlow() {
    colors.resize(num_colors);
    for (size_t ii = 0; ii < num_colors; ++ii) {
        double const angle = double(ii*2)*M_PI/num_colors;
        double const sin_a = std::sin(angle);
        double const cos_a = std::cos(angle);
        double const a = a_center + cos_a * color_radius;
        double const b = b_center + sin_a * color_radius;
        cv::Vec3f din(L, a, b);
        colors[ii] = ColorConvert::DIN2rgb(din);
    }
}

ColorFlow &ColorFlow::getInstance() {
    static ColorFlow instance;
    return instance;
}

cv::Vec3f ColorFlow::colorDirection(cv::Vec2f const& dir) {
    if (!std::isfinite(dir[0]) || !std::isfinite(dir[1])) {
        return cv::Vec3f(0,0,0);
    }
    double const angle = std::atan2(dir[1], dir[0]);
    size_t const index = size_t(std::round(((angle + M_PI)/(2*M_PI))*num_colors)) % num_colors;
    assert(index < num_colors);
    return getInstance().colors[index];
}

cv::Vec3b ColorFlow::colorFlow(cv::Vec2f const& flow) {
    double const length = cv::norm(flow);
    if (!std::isfinite(length) || !std::isfinite(flow[0]) || !std::isfinite(flow[1])) {
        return cv::Vec3b(0,0,0);
    }
    cv::Vec3f const dir = colorDirection(flow);
    if (length > 1.0000001) {
        return dir/2;
    }
    return (1.0-length) * cv::Vec3f(255,255,255) + length * dir;
}

cv::Mat_<cv::Vec3b> ColorFlow::colorFlow(cv::Mat_<cv::Vec2f> const& flow) {
    cv::Mat_<cv::Vec3b> result(flow.size(), cv::Vec3b(0,0,0));
    for (int row = 0; row < flow.rows; ++row) {
        for (int col = 0; col < flow.cols; ++col) {
            result(row, col) = colorFlow(flow(row, col));
        }
    }
    return result;
}

cv::Mat_<cv::Vec3b> ColorFlow::colorDirection(cv::Mat_<cv::Vec2f> const& flow) {
    cv::Mat_<cv::Vec3b> result(flow.size(), cv::Vec3b(0,0,0));
    for (int row = 0; row < flow.rows; ++row) {
        for (int col = 0; col < flow.cols; ++col) {
            result(row, col) = colorDirection(flow(row, col));
        }
    }
    return result;
}

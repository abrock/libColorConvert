#undef NDEBUG
#include <cassert>

#include <opencv2/optflow.hpp>
#include "colorconvert.h"
#include "colorflow.h"
#include <gtest/gtest.h>
#include <random>
#include <runningstats/runningstats.h>

cv::Mat_<cv::Vec2f> getCircleFlow(int const radius) {
    cv::Mat_<cv::Vec2f> result(2*radius+1, 2*radius+1, cv::Vec2f(0,0));
    for (int row = 0; row < result.rows; ++row) {
        for (int col = 0; col < result.cols; ++col) {
            cv::Vec2f flow(col-radius, row-radius);
            result(row, col) = flow/radius;
        }
    }
    return result;
}

TEST(ColorFlow, circle) {
    cv::Mat_<cv::Vec2f> circle = getCircleFlow(400);
    cv::writeOpticalFlow("circle.flo", circle);

    cv::imwrite("circle-own.png", ColorFlow::colorFlow(circle));
}

cv::Mat_<cv::Vec2f> getBandFlow() {
    cv::Mat_<cv::Vec2f> result(110,370, cv::Vec2f(0,0));
    for (int row = 0; row < result.rows; ++row) {
        for (int col = 0; col < result.cols; ++col) {
            double const angle = (2*M_PI*col)/360;
            double const cos = std::cos(angle);
            double const sin = std::sin(angle);
            double const factor = double(row)/100;
            cv::Vec2f flow(cos, sin);
            result(row, col) = flow*factor;
        }
    }
    return result;
}

TEST(ColorFlow, band) {
    cv::Mat_<cv::Vec2f> circle = getBandFlow();
    cv::writeOpticalFlow("band.flo", circle);

    cv::imwrite("band-own.png", ColorFlow::colorFlow(circle));
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
  return 0;
}

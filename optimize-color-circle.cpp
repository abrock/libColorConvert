#include "colorconvert.h"
#include "colorflow.h"
#include <iostream>
#include <ceres/ceres.h>
#include <vector>

struct Residuals {
    // Grey = 0.30*R + 0.59*G + 0.11*B.
    double getGrey(double const r, double const g, double const b) {
        return 0.30*r + 0.59*g + 0.11*b;
    }

    double grey = getGrey(0,0,255);

    bool operator()(
            double * residuals,
            double const * const self,
            double const * const left,
            double const * const right) {
        // First idea: distance from left to self = distance from self to right.
        // left - self = self - right
        // 0 = 2*self - left - right

        double self_din[3];
        double left_din[3];
        double right_din[3];

        ColorConvert::rgb2DIN(self[0], self[1], self[2], self_din[0], self_din[1], self_din[2]);
        ColorConvert::rgb2DIN(left[0], left[1], left[2], left_din[0], left_din[1], left_din[2]);
        ColorConvert::rgb2DIN(right[0], right[1], right[2], right_din[0], right_din[1], right_din[2]);

        for (size_t ii = 0; ii < 3; ++ii) {
            residuals[ii] = 2*self[ii] - left[ii] - right[ii];
        }

        // Second idea: distances should be large, but only in color.
        // Therefore: exp(-distance(self, left))

        residuals[3] = ceres::exp(-ceres::abs(self[1] - left[1]));
        residuals[4] = ceres::exp(-ceres::abs(self[2] - left[2]));

        // Third idea: grey value should be as close as possible to some known grey value.
        double const _grey = getGrey(self[0], self[1], self[2]);
        residuals[5] = _grey - grey;

        return true;
    }
};

int main(void) {
    ceres::Problem problem;

    size_t const N = 100;
    std::vector<cv::Vec3d> color_circle(N);
    for (size_t ii = 0; ii < N; ++ii) {
        double const alpha = double(ii)/N;
        cv::Vec3f val = ColorFlow::colorDirection(cv::Vec2f(std::sin(alpha), std::cos(alpha)));
        color_circle[ii] = val;
    }

    for (size_t ii = 0; ii < N; ++ii) {
        double const alpha = double(ii)/N;
        cv::Vec3f val = ColorFlow::colorDirection(cv::Vec2f(std::sin(alpha), std::cos(alpha)));
        color_circle[ii] = val;
    }

    for (size_t ii = 0; ii < N; ++ii) {

    }


}

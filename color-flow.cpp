#include <iostream>
#include "colorflow.h"
#include <tclap/CmdLine.h>
#include <opencv2/optflow.hpp>
#include <opencv2/imgproc.hpp>
#include <runningstats/runningstats.h>

#include "gnuplot-iostream.h"

namespace rs = runningstats;

void analyzeColorFile(std::string const& filename) {
    std::cout << "Analyzing file " << filename << std::endl;
    cv::Mat_<cv::Vec3b> const src = cv::imread(filename);
    if (src.rows <= 0 || src.cols <= 0) {
        std::cout << "File " << filename << " could not be read" << std::endl;
        return;
    }
    cv::Mat_<float> grey(src.size());
    cv::Mat_<cv::Vec3f> din(src.size());

    rs::Image2D<float> grey_img(1, 1);
    for (int row = 0; row < src.rows; ++row) {
        for (int col = 0; col < src.cols; ++col) {
            cv::Vec3f const rgb(src(row, col)[2], src(row ,col)[1], src(row, col)[0]);
            grey(row, col) = 0.30*rgb[0] + 0.59*rgb[1] + 0.11*rgb[2];
            din(row, col) = ColorConvert::rgb2DIN(rgb);
            grey(row, col) = din(row, col)[0]*2.55;
            grey_img[col][row] = grey(row, col);
        }
    }

    rs::Image2D<float> grey_dx(1, 1);
    rs::Image2D<float> grey_dy(1, 1);
    rs::Image2D<float> grey_min_row(1, 1);

    for (int row = 0; row < src.rows; ++row) {
        float min_grey_in_row = grey(row, 0);
        for (int col = 0; col+1 < src.cols; ++col) {
            min_grey_in_row = std::min(min_grey_in_row, grey(row, col));
        }
        for (int col = 0; col+1 < src.cols; ++col) {
            if (col+1 < src.cols) {
                grey_dx[col][row] = grey(row, col+1) - grey(row, col);
            }
            grey_min_row[col][row] = grey(row, col+1) - min_grey_in_row;
        }
    }

    rs::HistConfig config;
    config.setFixedRatio();
    config.setFlipY();
    grey_dx.plot(filename + "-grey-dx", config.clone().setTitle("d grey / dx"));
    grey_min_row.plot(filename + "-grey-min-row", config.clone().setTitle("Difference to Min."));
    grey_img.plot(filename + "-grey-img", config.clone().setTitle("Grey Value"));
}

void analyzeColorFiles(std::vector<std::string> const& filenames) {
    for (size_t ii = 0; ii < filenames.size(); ++ii) {
        analyzeColorFile(filenames[ii]);
    }
}

void gnuplotImage(std::string const& filename) {
    gnuplotio::Gnuplot plt("tee " + filename + ".gpl | gnuplot -persist");

    plt << "set term svg enhanced background rgb 'white';\n"
        << "set output '" << filename << "-gpl.svg;\n"
        << "set size ratio -1;\n"
        << "plot '" << filename << "' binary filetype=png with rgbimage;\n";
}

int main(int argc, char ** argv) {

    TCLAP::CmdLine cmd("color-flow tool");

    TCLAP::ValueArg<double> threshold_arg("t", "threshold", "threshold value", false, -1, "double");
    cmd.add(threshold_arg);

    TCLAP::SwitchArg analyze_arg("a", "analyze", "Input files are RGB files and derivatives should be analyzed");
    cmd.add(analyze_arg);

    TCLAP::UnlabeledMultiArg<std::string> input_files("", "flo-files", true, "", "filename, .flo file");
    cmd.add(input_files);

    cmd.parse(argc, argv);

    if (analyze_arg.isSet() && analyze_arg.getValue()) {
        analyzeColorFiles(input_files.getValue());
        return EXIT_SUCCESS;
    }

    double const thresh = threshold_arg.getValue();

    for (std::string const& file : input_files.getValue()) {
        cv::Mat_<cv::Vec2f> const flow = cv::readOpticalFlow(file);
        if (flow.rows <= 0 || flow.cols <= 0) {
            std::cout << "Flow file " << file << " not readable" << std::endl;
            continue;
        }
        runningstats::QuantileStats<float> lengths;
        for (cv::Vec2f const& v : flow) {
            double const length = cv::norm(v);
            if (std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(length)) {
                lengths.push_unsafe(length);
            }
        }
        if (thresh <= 0) {
            for (double const quantile : {.99, .995, .999}) {
                double const _thresh = lengths.getQuantile(quantile);
                cv::Mat_<cv::Vec2f> _flow = flow / _thresh;
                cv::Mat colored = ColorFlow::colorFlow(_flow);
                cv::cvtColor(colored, colored, cv::COLOR_RGB2BGR);
                cv::imwrite(file + "-" + std::to_string(_thresh) + ".png", colored);
                gnuplotImage(file + "-" + std::to_string(_thresh) + ".png");
            }
        }
        else {
            cv::Mat_<cv::Vec2f> _flow = flow / thresh;
            cv::Mat colored = ColorFlow::colorFlow(_flow);
            cv::cvtColor(colored, colored, cv::COLOR_RGB2BGR);
            cv::imwrite(file + "-" + std::to_string(thresh) + ".png", colored);
            gnuplotImage(file + "-" + std::to_string(thresh) + ".png");
        }
        cv::Mat colored = ColorFlow::colorDirection(flow);
        cv::cvtColor(colored, colored, cv::COLOR_RGB2BGR);
        cv::imwrite(file + "-directions.png", colored);
        gnuplotImage(file + "-directions.png");
    }
}

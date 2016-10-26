#include "colorconvert.h"


namespace ColorConvert {

/**
 * Read RGB-values from strings containing hexadecimal values like "ff0000" and store them in integer values.
 * This doesn't throw any exceptions, if the string doesn't make sense the behaviour is undefined.
 *
 * @param[in] rgb String containing the hexadecimal representation of the color.
 * @param[out] r Here the value for red is stored.
 * @param[out] g Here the value for green is stored.
 * @param[out] b Here the value for blue is stored.
 */
void string2rgb(const std::string& rgb, int& r, int& g, int& b) {
    r = b = g = 0;
    try {
        r = std::stoi(rgb.substr(0,2), 0, 16);
        g = std::stoi(rgb.substr(2,2), 0, 16);
        b = std::stoi(rgb.substr(4,2), 0, 16);
    }
    catch(std::exception& e){}// Who cares?
}

/**
 * Convert L*a*b-values to DIN99-values.
 *
 * @param[in] L L-value
 * @param[in] L a-value
 * @param[in] L b-value
 * @param[out] DIN_L calculated L-value in the DIN99 system.
 * @param[out] DIN_a calculated a-value in the DIN99 system.
 * @param[out] DIN_b calculated b-value in the DIN99 system.
 */
template<class T>
void lab2DIN(const T L, const T a, const T b, T& DIN_L, T& DIN_a, T& DIN_b) {

#define COS16 0.96126169593831886192
#define SIN16 0.27563735581699918561

    DIN_L = 105.51 * std::log(1 + 0.0158 * L);
    const T DIN_e = a * COS16 + b * SIN16;
    const T DIN_f = 0.7*(-a * SIN16 + b * COS16);

    const T DIN_G = std::sqrt(DIN_e * DIN_e + DIN_f * DIN_f);
    const T k = std::log(1 + 0.045 * DIN_G) / 0.045;

    if (DIN_G <= 0) {
        DIN_a = 0;
        DIN_b = 0;
    }
    else {
        DIN_a = k * DIN_e / DIN_G;
        DIN_b = k * DIN_f / DIN_G;
    }

}


template<class Source, class T>
void rgb2Lab(const Source r, const Source g, const Source b, T& dest_L, T& dest_a, T& dest_b) {
    // Create a matrix from the image.
    cv::Mat M(1, 1, CV_32FC3);
    // Write the normalized values into the matrix
    M.at<cv::Vec3f>(0,0)[0] = (float)r/255;
    M.at<cv::Vec3f>(0,0)[1] = (float)g/255;
    M.at<cv::Vec3f>(0,0)[2] = (float)b/255;

    // Color space conversion
    cv::cvtColor(M, M, cv::COLOR_RGB2Lab);

    // Reading the values from the matrix.
    dest_L = static_cast<T>(M.at<cv::Vec3f>(0,0)[0]);
    dest_a = static_cast<T>(M.at<cv::Vec3f>(0,0)[1]);
    dest_b = static_cast<T>(M.at<cv::Vec3f>(0,0)[2]);
}

template void rgb2Lab<double, double>(const double, const double, const double, double&, double&, double&);


template<class Source, class T>
void Lab2rgb(const Source L, const Source a, const Source b, T& dest_r, T& dest_g, T& dest_b) {
    // Create a matrix from the image.
    cv::Mat M(1, 1, CV_32FC3);
    // Write the normalized values into the matrix
    M.at<cv::Vec3f>(0,0)[0] = (float)L;
    M.at<cv::Vec3f>(0,0)[1] = (float)a;
    M.at<cv::Vec3f>(0,0)[2] = (float)b;

    // Color space conversion
    cv::cvtColor(M, M, cv::COLOR_Lab2RGB);

    // Reading the values from the matrix.
    dest_r = static_cast<T>(255.0 * M.at<cv::Vec3f>(0,0)[0]);
    dest_g = static_cast<T>(255.0 * M.at<cv::Vec3f>(0,0)[1]);
    dest_b = static_cast<T>(255.0 * M.at<cv::Vec3f>(0,0)[2]);
}

template void Lab2rgb<double, double>(const double, const double, const double, double&, double&, double&);

/**
 * Convert a string containing a hexadecimal RGB-representation of some color to a L*a*b-representation of the same color.
 *
 * @param[in] rgb String (hopefully) containing a hexadecimal RGB-representation of a color.
 * @param[out] dest_L Calculated L-value in the L*a*b-system
 * @param[out] dest_a Calculated a-value in the L*a*b-system
 * @param[out] dest_b Calculated b-value in the L*a*b-system
 */
template<class T>
void rgb2Lab(const std::string& rgb, T& dest_L, T& dest_a, T& dest_b) {
    // Get the R, G and B from the string
    int r = 0, g = 0, b = 0;
    string2rgb(rgb, r, g, b);
    rgb2Lab(r, g, b, dest_L, dest_a, dest_b);
}

template void rgb2Lab<float>(const std::string&, float&, float&, float&);

/**
 * Convert a string containing a hexadecimal RGB-representation of some color to a DIN99-representation of the same color.
 *
 * @param[in] rgb String (hopefully) containing a hexadecimal RGB-representation of a color.
 * @param[out] DIN_L calculated L-value in the DIN99 system.
 * @param[out] DIN_a calculated a-value in the DIN99 system.
 * @param[out] DIN_b calculated b-value in the DIN99 system.
 */
template<class T>
void rgb2DIN(const std::string& rgb, T& DIN_L, T& DIN_a, T& DIN_b) {
    float L = 0, a = 0, b = 0;
    rgb2Lab(rgb, L, a, b);
    lab2DIN(L, a, b, DIN_L, DIN_a, DIN_b);
}

/**
 * Calculate the square of a float.
 *
 * @param [in] a The value to square.
 * @return a*a
 */
template<class T>
T sqr(const T a) {
    return a*a;
}

/**
 * Calculate the difference of two colors in the L*a*b-system.
 *
 * @param[in] x String containing a hexadecimal RGB-representation of color #1.
 * @param[in] y String containing a hexadecimal RGB-representation of color #2.
 * @return Color difference
 */
double labDiff(const std::string& x, const std::string& y) {

    double x_L = 0, x_a = 0, x_b = 0, y_L = 0, y_a = 0, y_b = 0;

    rgb2Lab(x, x_L, x_a, x_b);
    rgb2Lab(y, y_L, y_a, y_b);

    return std::sqrt(sqr(x_L-y_L) + sqr(x_a-y_a) + sqr(x_b-y_b));
}

/**
 * Calculate the difference of two colors in the DIN99-system.
 *
 * @param[in] x String containing a hexadecimal RGB-representation of color #1.
 * @param[in] y String containing a hexadecimal RGB-representation of color #2.
 * @return Color difference
 */
double DINDiff(const std::string& x, const std::string& y) {

    float x_L = 0, x_a = 0, x_b = 0, y_L = 0, y_a = 0, y_b = 0;

    rgb2DIN(x, x_L, x_a, x_b);
    rgb2DIN(y, y_L, y_a, y_b);

    return std::sqrt(sqr(x_L-y_L) + sqr(x_a-y_a) + sqr(x_b-y_b));
}

}

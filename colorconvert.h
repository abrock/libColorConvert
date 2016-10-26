#ifndef COLORCONVERT_H
#define COLORCONVERT_H

#include <string>
#include <cmath>
#include <opencv/highgui.h>
#include <opencv/cv.hpp>

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
    void string2rgb(const std::string& rgb, int& r, int& g, int& b);

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
    void lab2DIN(const T L, const T a, const T b, T& DIN_L, T& DIN_a, T& DIN_b);

    /**
     * Convert a RGB-representation of some color to a L*a*b-representation of the same color.
     *
     * @param[in] r red value (0-255).
     * @param[in] g green value (0-255).
     * @param[in] b blue value (0-255).
     * @param[out] dest_L Calculated L-value in the L*a*b-system
     * @param[out] dest_a Calculated a-value in the L*a*b-system
     * @param[out] dest_b Calculated b-value in the L*a*b-system
     */
    template<class Source, class T>
    void rgb2Lab(const Source r, const Source g, const Source b, T& dest_L, T& dest_a, T& dest_b);

    /**
     * Convert a string containing a hexadecimal RGB-representation of some color to a L*a*b-representation of the same color.
     *
     * @param[in] rgb String (hopefully) containing a hexadecimal RGB-representation of a color.
     * @param[out] dest_L Calculated L-value in the L*a*b-system
     * @param[out] dest_a Calculated a-value in the L*a*b-system
     * @param[out] dest_b Calculated b-value in the L*a*b-system
     */
    template<class T>
    void rgb2Lab(const std::string& rgb, T& dest_L, T& dest_a, T& dest_b);

    /**
     * Convert a string containing a hexadecimal RGB-representation of some color to a DIN99-representation of the same color.
     *
     * @param[in] rgb String (hopefully) containing a hexadecimal RGB-representation of a color.
     * @param[out] DIN_L calculated L-value in the DIN99 system.
     * @param[out] DIN_a calculated a-value in the DIN99 system.
     * @param[out] DIN_b calculated b-value in the DIN99 system.
     */
    template<class T>
    void rgb2DIN(const std::string& rgb, T& DIN_L, T& DIN_a, T& DIN_b);

    /**
     * Calculate the square of a float.
     *
     * @param [in] a The value to square.
     * @return a*a
     */
    template<class T>
    T sqr(const T a);

    /**
     * Calculate the difference of two colors in the L*a*b-system.
     *
     * @param[in] x String containing a hexadecimal RGB-representation of color #1.
     * @param[in] y String containing a hexadecimal RGB-representation of color #2.
     * @return Color difference
     */
    double labDiff(const std::string& x, const std::string& y);

    /**
     * Calculate the difference of two colors in the DIN99-system.
     *
     * @param[in] x String containing a hexadecimal RGB-representation of color #1.
     * @param[in] y String containing a hexadecimal RGB-representation of color #2.
     * @return Color difference
     */
    double DINDiff(const std::string& x, const std::string& y);

    template<class Source, class T>
    void Lab2rgb(const Source L, const Source a, const Source b, T& dest_r, T& dest_g, T& dest_b);
}

#endif // COLORCONVERT_H

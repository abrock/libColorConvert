#include "colorconvert.h"
#include <gtest/gtest.h>
#include <random>

/**
 * This runs a set of simple tests, converting hexadecimal representations to integer values. If somebody breaks the implementation assertions fail.
 */
void test1() {
    int r = 0, g = 0, b = 0;
    std::string test("ff00ff");
    ColorConvert::string2rgb(test, r, g, b);

    assert(255 == r);
    assert(0 == g);
    assert(255 == b);

    test = std::string("5061BB");
    ColorConvert::string2rgb(test, r, g, b);

    assert(80 == r);
    assert(97 == g);
    assert(187 == b);
}

/**
 * Convert a hexadecimal representation of a RGB color to the L*a*b-system and assert the error is smaller than 2e-3.
 *
 * @param[in] data Hexadecimal representation of the RGB color.
 * @param[in] L known L-value corresponding to the data.
 * @param[in] a known a-value corresponding to the data.
 * @param[in] b known b-value corresponding to the data.
 */
void testRGB2Lab(const char* data, const float L, const float a, const float b) {
    std::string test(data);

    float M_L = 999, M_a = 999, M_b = 999;

    ColorConvert::rgb2Lab(test, M_L, M_a, M_b);

    const double error = (M_L - L)*(M_L - L) + (M_a - a)*(M_a - a) + (M_b - b)*(M_b - b);
    assert(error < 2e-3);
}

/**
 * Test the rgb2Lab function using the testRGB2Lab function.
 */
void test2() {
    testRGB2Lab("ffffff", 100, 0, 0);
    testRGB2Lab("777777", 50, 0, 0);
    testRGB2Lab("000000", 0, 0, 0);
    testRGB2Lab("ff0000", 53.2328817, 80.109309529822, 67.2200683102643);
    testRGB2Lab("00ff00", 87.7370334735442, -86.1846364976253, 83.1811647477785);
    testRGB2Lab("0000ff", 32.3025866672495, 79.1966617893094, -107.863681044952);
    testRGB2Lab("ffff00", 97.1382469812973, -21.5559083348323, 94.4824854464446);
}

/**
 * Test the labDiff method by feeding it two hexadecimal RGB representations and comparing the calculated difference to the known difference.
 *
 * @param[in] x Hexadecimal RBG representation of color #1
 * @param[in] y Hexadecimal RBG representation of color #2
 * @param[in] diff Known difference between x and y in the L*a*b-system
 */
void testLabDiff(const char* x, const char* y, const float diff) {
    std::string test_x(x);
    std::string test_y(y);

    const float diff2 = ColorConvert::labDiff(test_x, test_y);
    const float diff3 = ColorConvert::labDiff(test_y, test_x);

    assert(std::abs(diff2-diff) < 1e-4);
    assert(std::abs(diff3-diff) < 1e-4);
}

/**
 * Test the labDiff method by using the testLabDiff method
 */
void test3() {
   testLabDiff("000000", "000000", 0);
   testLabDiff("000f00", "000f00", 0);
   testLabDiff("f00f00", "f00f00", 0);
   testLabDiff("f01f00", "f01f00", 0);
   testLabDiff("f01f0b", "f01f0b", 0);
}

/**
 * Run all tests
 */
void test() {
    test1();
    test2();
    test3();
}

TEST(rgb2lab, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,255);
    for (size_t ii = 0; ii < 5; ++ii) {
        const double orig_r = distribution(generator);
        const double orig_g = distribution(generator);
        const double orig_b = distribution(generator);
        double r = 0, g = 0, b = 0;
        double lab_L = 0, lab_a = 0, lab_b = 0;
        ColorConvert::rgb2Lab(orig_r, orig_g, orig_b, lab_L, lab_a, lab_b);
        ColorConvert::Lab2rgb(lab_L, lab_a, lab_b, r, g, b);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(r, orig_r, 1e-2);
    }
}

int main(int argc, char** argv) {
  test();
  testing::InitGoogleTest(&argc, argv);
  std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
  return 0;
}

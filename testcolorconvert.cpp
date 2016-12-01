#include "colorconvert.h"
#include <gtest/gtest.h>
#include <random>
#include "libRunningStats/runningstats.h"

#define ZERO_PROBABILITY 0.03

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
    RunningStats diff_r, diff_g, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
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
        EXPECT_NEAR(r, orig_r, 1e-1);
        EXPECT_NEAR(g, orig_g, 1e-1);
        EXPECT_NEAR(b, orig_b, 1e-1);
        diff_r.push(std::abs(r-orig_r));
        diff_g.push(std::abs(g-orig_g));
        diff_b.push(std::abs(b-orig_b));
    }
    std::cout << "Stats for rgb->lab->rgb: " << std::endl
              << "Errors in r: " << diff_r.print() << std::endl
              << "Errors in g: " << diff_g.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}



#if 0
TEST(Lab2rgb, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 5e3; ++ii) {
        const double orig_L = dist_L(generator);
        const double orig_a = dist_a(generator);
        const double orig_b = dist_b(generator);
        double r = 0, g = 0, b = 0;
        double lab_L = 0, lab_a = 0, lab_b = 0;
        ColorConvert::Lab2rgb(orig_L, orig_a, orig_b, r, g, b);
        ColorConvert::rgb2Lab(r, g, b, lab_L, lab_a, lab_b);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(lab_L, orig_L, 1e-1);
        EXPECT_NEAR(lab_a, orig_a, 1e-1);
        EXPECT_NEAR(lab_b, orig_b, 1e-1);
        diff_L.push(std::abs(lab_L-orig_L));
        diff_a.push(std::abs(lab_a-orig_a));
        diff_b.push(std::abs(lab_b-orig_b));
    }
    std::cout << "Stats for lab->rgb->lab: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}
#endif

TEST(Lab2DIN, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const double orig_L = dist_L(generator);
        const double orig_a = dist_a(generator);
        const double orig_b = dist_b(generator);
        double DIN_L = 0, DIN_a = 0, DIN_b = 0;
        double lab_L = 0, lab_a = 0, lab_b = 0;
        ColorConvert::Lab2DIN(orig_L, orig_a, orig_b, DIN_L, DIN_a, DIN_b);
        ColorConvert::DIN2Lab(DIN_L, DIN_a, DIN_b, lab_L, lab_a, lab_b);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(lab_L, orig_L, 1e-10);
        EXPECT_NEAR(lab_a, orig_a, 1e-10);
        EXPECT_NEAR(lab_b, orig_b, 1e-10);
        diff_L.push(std::abs(lab_L-orig_L));
        diff_a.push(std::abs(lab_a-orig_a));
        diff_b.push(std::abs(lab_b-orig_b));
    }
    std::cout << "Stats for Lab->DIN->Lab: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(DIN2Lab, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const double orig_L = dist_L(generator);
        const double orig_a = dist_a(generator);
        const double orig_b = dist_b(generator);
        double Lab_L = 0, Lab_a = 0, Lab_b = 0;
        double DIN_L = 0, DIN_a = 0, DIN_b = 0;
        ColorConvert::DIN2Lab(orig_L, orig_a, orig_b, Lab_L, Lab_a, Lab_b);
        ColorConvert::Lab2DIN(Lab_L, Lab_a, Lab_b, DIN_L, DIN_a, DIN_b);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(DIN_L, orig_L, 1e-10);
        EXPECT_NEAR(DIN_a, orig_a, 1e-10);
        EXPECT_NEAR(DIN_b, orig_b, 1e-10);
        diff_L.push(std::abs(DIN_L-orig_L));
        diff_a.push(std::abs(DIN_a-orig_a));
        diff_b.push(std::abs(DIN_b-orig_b));
    }
    std::cout << "Stats for DIN->Lab->DIN: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(DIN2LabVec3f, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const cv::Vec3f orig(dist_L(generator), dist_a(generator), dist_b(generator));
        const cv::Vec3f intermediate = ColorConvert::DIN2Lab(orig);
        const cv::Vec3f final = ColorConvert::Lab2DIN(intermediate);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(orig[0], final[0], 1e-4);
        EXPECT_NEAR(orig[1], final[1], 1e-4);
        EXPECT_NEAR(orig[2], final[2], 1e-4);
        diff_L.push(std::abs(orig[0] - final[0]));
        diff_a.push(std::abs(orig[1] - final[1]));
        diff_b.push(std::abs(orig[2] - final[2]));
    }
    std::cout << "Stats for vectorized DIN->Lab->DIN: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(DIN2LabVec3d, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    std::uniform_real_distribution<double> dist(0, 1);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const double L_sign = dist(generator) < ZERO_PROBABILITY ? 0 : 1;
        const double a_sign = dist(generator) < ZERO_PROBABILITY ? 0 : 1;
        const double b_sign = dist(generator) < ZERO_PROBABILITY ? 0 : 1;
        const cv::Vec3d orig(L_sign * dist_L(generator), a_sign * dist_a(generator), b_sign * dist_b(generator));
        const cv::Vec3d intermediate = ColorConvert::DIN2Lab(orig);
        const cv::Vec3d final = ColorConvert::Lab2DIN(intermediate);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(orig[0], final[0], 1e-10);
        EXPECT_NEAR(orig[1], final[1], 1e-10);
        EXPECT_NEAR(orig[2], final[2], 1e-10);
        diff_L.push(std::abs(orig[0] - final[0]));
        diff_a.push(std::abs(orig[1] - final[1]));
        diff_b.push(std::abs(orig[2] - final[2]));
    }
    std::cout << "Stats for vectorized DIN->Lab->DIN: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(Lab2DINVec3f, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const cv::Vec3f orig(dist_L(generator), dist_a(generator), dist_b(generator));
        const cv::Vec3f intermediate = ColorConvert::Lab2DIN(orig);
        const cv::Vec3f final = ColorConvert::DIN2Lab(intermediate);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(orig[0], final[0], 1e-4);
        EXPECT_NEAR(orig[1], final[1], 1e-4);
        EXPECT_NEAR(orig[2], final[2], 1e-4);
        diff_L.push(std::abs(orig[0] - final[0]));
        diff_a.push(std::abs(orig[1] - final[1]));
        diff_b.push(std::abs(orig[2] - final[2]));
    }
    std::cout << "Stats for vectorized Lab->DIN->Lab: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(Lab2DINScalar, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const cv::Scalar orig(dist_L(generator), dist_a(generator), dist_b(generator));
        const cv::Scalar intermediate = ColorConvert::Lab2DIN(orig);
        const cv::Scalar final = ColorConvert::DIN2Lab(intermediate);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(orig[0], final[0], 1e-10);
        EXPECT_NEAR(orig[1], final[1], 1e-10);
        EXPECT_NEAR(orig[2], final[2], 1e-10);
        diff_L.push(std::abs(orig[0] - final[0]));
        diff_a.push(std::abs(orig[1] - final[1]));
        diff_b.push(std::abs(orig[2] - final[2]));
    }
    std::cout << "Stats for vectorized Lab->DIN->Lab: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(Lab2DINVec3d, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_L(0,100);
    std::uniform_real_distribution<double> dist_a(-170,100);
    std::uniform_real_distribution<double> dist_b(-100,150);
    std::uniform_real_distribution<double> dist(0,1);
    RunningStats diff_L, diff_a, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const double L_sign = dist(generator) < ZERO_PROBABILITY ? 0 : 1;
        const double a_sign = dist(generator) < ZERO_PROBABILITY ? 0 : 1;
        const double b_sign = dist(generator) < ZERO_PROBABILITY ? 0 : 1;
        const cv::Vec3d orig(L_sign * dist_L(generator), a_sign * dist_a(generator), b_sign * dist_b(generator));
        const cv::Vec3d intermediate = ColorConvert::Lab2DIN(orig);
        const cv::Vec3d final = ColorConvert::DIN2Lab(intermediate);
#if 0
        std::cout << "Orig:   " << orig_r << ", " << orig_g << ", " << orig_b << std::endl;
        std::cout << "Lab:    " << lab_L << ", " << lab_a << ", " << lab_b << std::endl;
        std::cout << "Result: " << r << ", " << g << ", " << b << std::endl;
        std::cout << std::endl;
#endif
        EXPECT_NEAR(orig[0], final[0], 1e-10);
        EXPECT_NEAR(orig[1], final[1], 1e-10);
        EXPECT_NEAR(orig[2], final[2], 1e-10);
        diff_L.push(std::abs(orig[0] - final[0]));
        diff_a.push(std::abs(orig[1] - final[1]));
        diff_b.push(std::abs(orig[2] - final[2]));
    }
    std::cout << "Stats for vectorized Lab->DIN->Lab: " << std::endl
              << "Errors in L: " << diff_L.print() << std::endl
              << "Errors in a: " << diff_a.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(rgb2LabVec3d, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist(0,255);
    RunningStats diff_r, diff_g, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const double r_sign = dist(generator)/255 < ZERO_PROBABILITY ? 0 : 1;
        const double g_sign = dist(generator)/255 < ZERO_PROBABILITY ? 0 : 1;
        const double b_sign = dist(generator)/255 < ZERO_PROBABILITY ? 0 : 1;
        const cv::Vec3d orig(r_sign * dist(generator), g_sign * dist(generator), b_sign * dist(generator));
        const cv::Vec3d intermediate = ColorConvert::rgb2Lab(orig);
        const cv::Vec3d final = ColorConvert::Lab2rgb(intermediate);

        EXPECT_NEAR(orig[0], final[0], 1e-1);
        EXPECT_NEAR(orig[1], final[1], 1e-1);
        EXPECT_NEAR(orig[2], final[2], 1e-1);
        diff_r.push(std::abs(orig[0] - final[0]));
        diff_g.push(std::abs(orig[1] - final[1]));
        diff_b.push(std::abs(orig[2] - final[2]));
    }
    std::cout << "Stats for vectorized RGB->Lab->RGB: " << std::endl
              << "Errors in r: " << diff_r.print() << std::endl
              << "Errors in g: " << diff_g.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

TEST(rgb2DINVec3d, bijective) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist(0,255);
    RunningStats diff_r, diff_g, diff_b;
    for (size_t ii = 0; ii < 1e5; ++ii) {
        const double r_sign = dist(generator)/255 < ZERO_PROBABILITY ? 0 : 1;
        const double g_sign = dist(generator)/255 < ZERO_PROBABILITY ? 0 : 1;
        const double b_sign = dist(generator)/255 < ZERO_PROBABILITY ? 0 : 1;
        const cv::Vec3d orig(r_sign * dist(generator), g_sign * dist(generator), b_sign * dist(generator));
        const cv::Vec3d intermediate = ColorConvert::rgb2DIN(orig);
        const cv::Vec3d final = ColorConvert::DIN2rgb(intermediate);

        EXPECT_NEAR(orig[0], final[0], 1e-1);
        EXPECT_NEAR(orig[1], final[1], 1e-1);
        EXPECT_NEAR(orig[2], final[2], 1e-1);
        diff_r.push(std::abs(orig[0] - final[0]));
        diff_g.push(std::abs(orig[1] - final[1]));
        diff_b.push(std::abs(orig[2] - final[2]));
    }
    std::cout << "Stats for vectorized RGB->DIN->RGB: " << std::endl
              << "Errors in r: " << diff_r.print() << std::endl
              << "Errors in g: " << diff_g.print() << std::endl
              << "Errors in b: " << diff_b.print() << std::endl;
    std::cout << std::endl;
}

int main(int argc, char** argv) {
  test();
  testing::InitGoogleTest(&argc, argv);
  cv::Scalar test;
  cv::Scalar test2 = ColorConvert::rgb2DIN(test);
  std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
  return 0;
}

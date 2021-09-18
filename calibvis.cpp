// Copyright 2019 ETH Zürich, Thomas Schöps
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "calibvis.h"

CalibVis::CalibVis()
{

}

#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Eigen>

#include <boost/polygon/voronoi.hpp>

#include <opencv2/highgui.hpp>

namespace calibvis {

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

using std::min;
using std::max;
using std::isfinite;

typedef Eigen::Vector2d Vec2d;
typedef Eigen::Vector2f Vec2f;
typedef Eigen::Vector3f Vec3f;

typedef size_t usize;
typedef int64_t i64;
typedef uint64_t u64;
typedef int32_t i32;
typedef uint32_t u32;
typedef int16_t i16;
typedef uint16_t u16;
typedef int8_t i8;
typedef uint8_t u8;

using std::vector;

template <typename T> using Image = cv::Mat_<T>;

/// Determinant of 2x2 matrix.
template <typename T>
inline T Det(T a, T b, T c, T d) {
    return a*d - b*c;
}

/// Computes the intersection point between the two given 2D lines a and b. Each
/// line is specified with two points (0 and 1) on the line. Returns false if no
/// intersection was found.
/// Adapted from: https://gist.github.com/TimSC/47203a0f5f15293d2099507ba5da44e6
template <typename T>
bool LineLineIntersection(
        const Eigen::Matrix<T, 2, 1>& a0,
        const Eigen::Matrix<T, 2, 1>& a1,
        const Eigen::Matrix<T, 2, 1>& b0,
        const Eigen::Matrix<T, 2, 1>& b1,
        Eigen::Matrix<T, 2, 1>* result) {
    // http://mathworld.wolfram.com/Line-LineIntersection.html
    T detL1 = Det(a0.x(), a0.y(), a1.x(), a1.y());
    T detL2 = Det(b0.x(), b0.y(), b1.x(), b1.y());
    T x1mx2 = a0.x() - a1.x();
    T x3mx4 = b0.x() - b1.x();
    T y1my2 = a0.y() - a1.y();
    T y3my4 = b0.y() - b1.y();

    T xnom = Det(detL1, x1mx2, detL2, x3mx4);
    T ynom = Det(detL1, y1my2, detL2, y3my4);
    T denom = Det(x1mx2, y1my2, x3mx4, y3my4);
    if (denom == 0) {
        // Lines don't seem to cross
        return false;
    }

    result->x() = xnom / denom;
    result->y() = ynom / denom;
    if (!isfinite(result->x()) || !isfinite(result->y())) {
        // Probably a numerical issue
        return false;
    }

    return true;
}

/// Determines the winding order of the given convex polygon. Given a
/// right-handed coordinate system, returns 1 for counter-clockwise and -1 for clockwise.
/// NOTE: This function requires all edges in the polygon to have a significant
///       length (that is not close to zero), and no successive parallel line
///       segments are allowed.
/// Implementation from: https://en.wikipedia.org/wiki/Curve_orientation#Practical_considerations
template <typename T>
int ConvexPolygonOrientation(const vector<Eigen::Matrix<T, 2, 1>>& polygon) {
    typedef Eigen::Matrix<T, 2, 1> PointT;
    const PointT& a = polygon[0];
    const PointT& b = polygon[1];
    const PointT& c = polygon[2];

    float det = (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y());
    return (det > 0) ? 1 : -1;
}

/// Computes whether the given 2D point lies within the given 2D polygon (given
/// as vector of corner points). The first point of the polygon does not need to
/// be repeated as its last point.
/// For self-intersecting polygons, this uses the 'alternating' fill rule.
/// NOTE: Several values could be pre-computed to speed up this computation if
///       performing it many times for the same polygon.
template <typename T>
bool PointInsidePolygon(
        const Eigen::Matrix<T, 2, 1>& point,
        const vector<Eigen::Matrix<T, 2, 1>>& polygon) {
    typedef Eigen::Matrix<T, 2, 1> PointT;

    int right_of_edge_count = 0;

    usize size = polygon.size();
    int prev_i = size - 1;
    for (int i = 0; i < size; ++ i) {
        const PointT& edge_start = polygon[prev_i];
        const PointT& edge_end = polygon[i];

        T min_y = min(edge_start.y(), edge_end.y());
        T max_y = max(edge_start.y(), edge_end.y());
        if (point.y() >= min_y && point.y() < max_y) {
            T min_x = min(edge_start.x(), edge_end.x());
            T max_x = max(edge_start.x(), edge_end.x());
            if (point.x() >= min_x) {
                if (point.x() < max_x) {
                    float relative_x = (point.x() - edge_start.x()) / (edge_end.x() - edge_start.x());
                    float relative_y = (point.y() - edge_start.y()) / (edge_end.y() - edge_start.y());
                    if ((edge_end.x() > edge_start.x() && relative_x > relative_y) ||
                            (edge_end.x() <= edge_start.x() && relative_x < relative_y)) {
                        ++ right_of_edge_count;
                    }
                } else {
                    ++ right_of_edge_count;
                }
            }
        }

        prev_i = i;
    }

    return right_of_edge_count & 1;
}

/// Computes the given polygon's area. For self-intersecting polygons, parts
/// with differing orientation count against each other, reducing the overall
/// reported area. Parts with equal orientation sum up.
template <typename T>
T PolygonArea(const vector<Eigen::Matrix<T, 2, 1>>& polygon) {
    T result = 0;
    usize size = polygon.size();
    int prev_i = size - 1;
    for (int i = 0; i < size; ++ i) {
        result += (polygon[i].x() - polygon[prev_i].x()) * (polygon[i].y() + polygon[prev_i].y());
        prev_i = i;
    }
    return fabs(0.5f * result);
}

// Clips the given polygon to the area of the given clip polygon, which must be
// convex. Implements the Sutherland–Hodgman algorithm:
// https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
template <typename T>
void ConvexClipPolygon(
        const vector<Eigen::Matrix<T, 2, 1>>& polygon,
        const vector<Eigen::Matrix<T, 2, 1>>& convex_clip,
        vector<Eigen::Matrix<T, 2, 1>>* output) {
    typedef Eigen::Matrix<T, 2, 1> PointT;

    usize clip_size = convex_clip.size();

    // Determine the winding order of convex_clip (clockwise or counter-clockwise)
    int clip_orientation = ConvexPolygonOrientation(convex_clip);

    // Loop over all edges in the clip path
    for (usize clip_edge_start = 0; clip_edge_start < clip_size; ++ clip_edge_start) {
        usize clip_edge_end = (clip_edge_start + 1) % clip_size;

        PointT edge = convex_clip[clip_edge_end] - convex_clip[clip_edge_start];
        PointT right_vector = PointT(edge.y(), -edge.x());
        float edge_right = right_vector.dot(convex_clip[clip_edge_end]);

        // Clip the polygon with this edge, reading from input and writing to clipped
        const vector<PointT>* input = (clip_edge_start == 0) ? &polygon : output;
        vector<PointT> clipped;

        // Loop over input points
        usize input_size = input->size();
        for (usize i = 0; i < input_size; ++ i) {
            const PointT& current_point = input->at(i);
            const PointT& prev_point = input->at((i + input_size - 1) % input_size);

            if (clip_orientation * ((right_vector.dot(current_point) > edge_right) ? 1 : -1) < 0) {  // current_point in convex_clip?
                if (clip_orientation * ((right_vector.dot(prev_point) > edge_right) ? 1 : -1) > 0) {  // prev_point not in convex_clip?
                    PointT intersection = prev_point;
                    LineLineIntersection(prev_point, current_point, convex_clip[clip_edge_end], convex_clip[clip_edge_start], &intersection);
                    clipped.push_back(intersection);
                }
                clipped.push_back(current_point);
            } else if (clip_orientation * ((right_vector.dot(prev_point) > edge_right) ? 1 : -1) < 0) {  // prev_point in convex_clip?
                PointT intersection = prev_point;
                LineLineIntersection(prev_point, current_point, convex_clip[clip_edge_end], convex_clip[clip_edge_start], &intersection);
                clipped.push_back(intersection);
            }
        }

        *output = clipped;  // TODO: avoid copies
    }
}

void CreateVoronoiDiagram(
        int const width,
        int const height,
        const vector<Vec2d>& reprojection_errors,
        const vector<Vec2f>& features,
        float* voronoi_integer_inv_scaling,
        voronoi_diagram<double>* vd,
        vector<Vec2f>* v_reprojection_errors,
        vector<boost::polygon::point_data<i32>>* v_points) {
    constexpr int kVoronoiIntegerScaling = 4;
    *voronoi_integer_inv_scaling = 1.f / kVoronoiIntegerScaling;
    v_reprojection_errors->reserve(32000);
    v_points->reserve(32000);
    Image<u8> v_point_image(width,height);
    v_point_image.setTo(static_cast<u8>(0));

    for (usize i = 0; i < reprojection_errors.size(); ++ i) {
        const Vec2d& reprojection_error = reprojection_errors[i];
        const Vec2f& feature = features[i];

        i32 ix = static_cast<i32>(feature.x());
        i32 iy = static_cast<i32>(feature.y());
        if (v_point_image(ix, iy) == 0) {
            v_points->push_back(boost::polygon::point_data<i32>(
                                    static_cast<i32>(kVoronoiIntegerScaling * feature.x()),
                                    static_cast<i32>(kVoronoiIntegerScaling * feature.y())));
            v_reprojection_errors->push_back(reprojection_error.cast<float>());
            v_point_image(ix, iy) = 1;
        }
    }
    assert(v_reprojection_errors->size() == v_points->size());

    // Create a voronoi diagram of all feature observations.
    // NOTE: I also tried JCash/voronoi from GitHub, which makes it easy to draw
    //       all the cells (also ones with infinite edges), but it seemed like
    //       it returned broken results in some cases.
    construct_voronoi(v_points->begin(), v_points->end(), vd);
}

void ClipVoronoiEdge(
        const vector<boost::polygon::point_data<i32>>& v_points,
        const voronoi_diagram<double>::edge_type* edge,
        std::vector<boost::polygon::point_data<double>>* clipped_edge) {
    if (edge->is_infinite()) {
        const voronoi_diagram<double>::cell_type& cell1 = *edge->cell();
        const voronoi_diagram<double>::cell_type& cell2 = *edge->twin()->cell();
        boost::polygon::point_data<double> origin, direction;
        // Infinite edges could not be created by two segment sites.
        if (cell1.contains_point() && cell2.contains_point()) {
            boost::polygon::point_data<double> p1 = v_points[cell1.source_index()];
            boost::polygon::point_data<double> p2 = v_points[cell2.source_index()];
            origin.x((p1.x() + p2.x()) * 0.5);
            origin.y((p1.y() + p2.y()) * 0.5);
            direction.x(p1.y() - p2.y());
            direction.y(p2.x() - p1.x());
        } else {
            std::cout << "Case not handled" << std::endl;
        }

        double large_value = 99999;  // TODO: Compute a reasonable value. It must place the clipped vertex beyond the other vertex of the edge, and also outside of the bounding box of what we are interested in. Or treat it as a direction with infinite length.
        if (edge->vertex0() == NULL) {
            clipped_edge->push_back(boost::polygon::point_data<double>(
                                        origin.x() - direction.x() * large_value,
                                        origin.y() - direction.y() * large_value));
        } else {
            clipped_edge->push_back(
                        boost::polygon::point_data<double>(edge->vertex0()->x(), edge->vertex0()->y()));
        }
        if (edge->vertex1() == NULL) {
            clipped_edge->push_back(boost::polygon::point_data<double>(
                                        origin.x() + direction.x() * large_value,
                                        origin.y() + direction.y() * large_value));
        } else {
            clipped_edge->push_back(
                        boost::polygon::point_data<double>(edge->vertex1()->x(), edge->vertex1()->y()));
        }
    } else {
        clipped_edge->push_back(
                    boost::polygon::point_data<double>(edge->vertex0()->x(), edge->vertex0()->y()));
        clipped_edge->push_back(
                    boost::polygon::point_data<double>(edge->vertex1()->x(), edge->vertex1()->y()));
    }
};

void RenderVoronoiTriangle(
        float voronoi_integer_inv_scaling,
        const vector<boost::polygon::point_data<i32>>& v_points,
        int point_index,
        const std::vector<boost::polygon::point_data<double>>& rendered_edge,
        const Vec3f& color,
        Image<Vec3f>* rendering) {
    vector<Vec2d> triangle(3);
    triangle[0] = voronoi_integer_inv_scaling * Vec2d(v_points[point_index].x(), v_points[point_index].y());
    triangle[1] = voronoi_integer_inv_scaling * Vec2d(rendered_edge[0].x(), rendered_edge[0].y());
    triangle[2] = voronoi_integer_inv_scaling * Vec2d(rendered_edge[1].x(), rendered_edge[1].y());

    // Compute the bounding box of the triangle
    Eigen::AlignedBox2d bbox;
    for (const Vec2d& point : triangle) {
        bbox.extend(point);
    }

    // Loop over all pixels that intersect the bounding box
    int min_x = max<int>(0, bbox.min().x());
    int max_x = min<int>(rendering->width() - 1, bbox.max().x());
    int min_y = max<int>(0, bbox.min().y());
    int max_y = min<int>(rendering->rows - 1, bbox.max().y());
    for (int y = min_y; y <= max_y; ++ y) {
        for (int x = min_x; x <= max_x; ++ x) {
            // Intersect the pixel area and the triangle
            vector<Vec2d> pixel_area(4);
            pixel_area[0] = Vec2d(x, y);
            pixel_area[1] = Vec2d(x + 1, y);
            pixel_area[2] = Vec2d(x + 1, y + 1);
            pixel_area[3] = Vec2d(x, y + 1);
            vector<Vec2d> intersection;
            ConvexClipPolygon(
                        triangle,
                        pixel_area,
                        &intersection);
            (*rendering)(x, y) += PolygonArea(intersection) * color;
        }
    }
};

template <class ColorComputer>
void RenderVoronoiDiagram(
        int width, int height,
        float voronoi_integer_inv_scaling,
        const voronoi_diagram<double>& vd,
        const vector<boost::polygon::point_data<i32>>& v_points,
        const ColorComputer& color_computer,
        Image<Vec3u8>* rendering) {
    Image<Vec3f> reprojection_error_direction_rendering(width, height);
    reprojection_error_direction_rendering.setTo(Vec3f::Zero());

    for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin(); it != vd.cells().end(); ++it) {
        const voronoi_diagram<double>::cell_type& cell = *it;
        const voronoi_diagram<double>::edge_type* edge = cell.incident_edge();
        if (!it->contains_point()) {
            std::cout << "Voronoi cell without a point, this is not handled." << std::endl;
            continue;
        }

        // Iterate over all edges of the Voronoi cell.
        do {
            if (!edge->is_primary()) {
                std::cout << "Non-primary edge of a Voronoi cell, this is not handled." << std::endl;
            }

            std::vector<boost::polygon::point_data<double>> clipped_edge;
            ClipVoronoiEdge(v_points, edge, &clipped_edge);

            Vec3f color = color_computer(it->source_index());

            // Render a triangle between this edge and the source point
            RenderVoronoiTriangle(voronoi_integer_inv_scaling, v_points, it->source_index(), clipped_edge, color, &reprojection_error_direction_rendering);

            // If there are two infinite edges after each other, we have to insert an additional triangle between them.
            // TODO: For areas at the image corners, this might not always work, we would need to insert more triangles in these cases
            const voronoi_diagram<double>::edge_type* next_edge = edge->next();
            if (edge->is_infinite() && next_edge->is_infinite() &&
                    (next_edge->next() != edge || next_edge != cell.incident_edge())) {  // avoid double drawing if there are only two edges in the cell
                std::vector<boost::polygon::point_data<double>> next_clipped_edge;
                ClipVoronoiEdge(v_points, next_edge, &next_clipped_edge);

                std::vector<boost::polygon::point_data<double>> clipped_vertices;
                if (!edge->vertex0()) {
                    clipped_vertices.push_back(clipped_edge[0]);
                } else {  // !edge->vertex1()
                    clipped_vertices.push_back(clipped_edge[1]);
                }
                if (!next_edge->vertex0()) {
                    clipped_vertices.push_back(next_clipped_edge[0]);
                } else {  // !next_edge->vertex1()
                    clipped_vertices.push_back(next_clipped_edge[1]);
                }

                RenderVoronoiTriangle(voronoi_integer_inv_scaling, v_points, it->source_index(), clipped_vertices, color, &reprojection_error_direction_rendering);
            }

            edge = edge->next();
        } while (edge != cell.incident_edge());
    }

    rendering->SetSize(reprojection_error_direction_rendering.size());
    for (int y = 0; y < rendering->height(); ++ y) {
        for (int x = 0; x < rendering->width(); ++ x) {
            (*rendering)(x, y) = (reprojection_error_direction_rendering(x, y) + Vec3f::Constant(0.5f)).cwiseMax(Vec3f::Zero()).cwiseMin(Vec3f::Constant(255.99f)).cast<u8>();
        }
    }
}

struct ReprojectionDirectionColorComputer {
    Vec3f operator() (int point_index) const {
        const Vec2f& reprojection_error = v_reprojection_errors[point_index];
        double dir = atan2(reprojection_error.y(), reprojection_error.x());  // from -M_PI to M_PI
        return Vec3f(
                    127 + 127 * sin(dir),
                    127 + 127 * cos(dir),
                    127);
    }

    const vector<Vec2f>& v_reprojection_errors;
};

void VisualizeReprojectionErrorDirections(
        const CameraModel* cam,
        float voronoi_integer_inv_scaling,
        const voronoi_diagram<double>& vd,
        const vector<Vec2f>& v_reprojection_errors,
        const vector<boost::polygon::point_data<i32>>& v_points,
        Image<Vec3u8>* direction_image) {
    RenderVoronoiDiagram(
                cam->width(), cam->height(),
                voronoi_integer_inv_scaling,
                vd,
                v_points,
                ReprojectionDirectionColorComputer{v_reprojection_errors},
                direction_image);
}

struct ReprojectionMagnitudeColorComputer {
    Vec3f operator() (int point_index) const {
        const Vec2f& reprojection_error = v_reprojection_errors[point_index];

        double factor = std::min(1., reprojection_error.norm() / max_error);
        return Vec3f(255.99f * factor, 255.99f * (1 - factor), 0);
    }

    double max_error;
    const vector<Vec2f>& v_reprojection_errors;
};

void VisualizeReprojectionErrorMagnitudes(
        const CameraModel* cam,
        float voronoi_integer_inv_scaling,
        const voronoi_diagram<double>& vd,
        const vector<Vec2f>& v_reprojection_errors,
        const vector<boost::polygon::point_data<i32>>& v_points,
        double max_error,
        Image<Vec3u8>* magnitude_image) {
    RenderVoronoiDiagram(
                cam->width(), cam->height(),
                voronoi_integer_inv_scaling,
                vd,
                v_points,
                ReprojectionMagnitudeColorComputer{max_error, v_reprojection_errors},
                magnitude_image);
}


/// Attempts to approximately compute the horizontal and vertical FOV, measured
/// through the center of the camera image. Returns -1 for a value if it cannot
/// be computed.
void ComputeApproximateFOV(
        const CameraModel* cam,
        double* horizontal_fov,
        double* vertical_fov) {
    *horizontal_fov = -1;
    *vertical_fov = -1;

    if (cam->type() == CameraModel::Type::NoncentralGeneric) {
        return;
    }

    // Horizontal FOV
    float min_x = cam->calibration_min_x() + 0.5f;
    float max_x = cam->calibration_max_x() + 0.5f;
    float y = 0.5f * cam->height();

    Vec3d left_unprojection;
    Vec3d right_unprojection;
    if (cam->Unproject(min_x, y, &left_unprojection) &&
            cam->Unproject(max_x, y, &right_unprojection)) {
        *horizontal_fov = acos(left_unprojection.normalized().dot(right_unprojection.normalized())) *
                (cam->width() / (max_x - min_x));
    }

    // Vertical FOV
    float min_y = cam->calibration_min_y() + 0.5f;
    float max_y = cam->calibration_max_y() + 0.5f;
    float x = 0.5f * cam->width();

    Vec3d top_unprojection;
    Vec3d bottom_unprojection;
    if (cam->Unproject(x, min_y, &top_unprojection) &&
            cam->Unproject(x, max_y, &bottom_unprojection)) {
        *vertical_fov = acos(top_unprojection.normalized().dot(bottom_unprojection.normalized())) *
                (cam->height() / (max_y - min_y));
    }
}


bool WriteReportInfoFile(
        const string& path,
        const CameraModel* cam,
        double horizontal_fov,
        double vertical_fov,
        int imageset_count,
        int num_localized_images,
        const vector<Vec2d>& reprojection_errors,
        usize reprojection_error_count,
        double reprojection_error_sum,
        double reprojection_error_max,
        double biasedness,
        double histogram_extent_in_px,
        double max_error_in_px) {
    // In a text file, list:
    // - The number of residuals used to create the report
    // - The average and maximum reprojection errors in pixels
    ofstream stream(path, std::ios::out);
    if (!stream) {
        return false;
    }
    stream << std::setprecision(14);

    stream << "resolution : " << cam->width() << " x " << cam->height() << std::endl;
    if (horizontal_fov >= 0) {
        stream << "horizontal_fov : " << (180.f / M_PI * horizontal_fov) << std::endl;
    }
    if (vertical_fov >= 0) {
        stream << "vertical_fov : " << (180.f / M_PI * vertical_fov) << std::endl;
    }
    stream << "" << std::endl;

    stream << "num_localized_imagesets : " << num_localized_images << std::endl;
    stream << "num_total_imagesets : " << imageset_count << std::endl;
    stream << "" << std::endl;

    stream << "reprojection_error_count : " << reprojection_error_count << std::endl;
    if (!reprojection_errors.empty()) {
        vector<double> reprojection_error_magnitudes(reprojection_errors.size());
        for (int i = 0; i < reprojection_errors.size(); ++ i) {
            reprojection_error_magnitudes[i] = reprojection_errors[i].norm();
        }

        std::sort(reprojection_error_magnitudes.begin(), reprojection_error_magnitudes.end());  // NOTE: No need to sort the whole vector just to get the median
        stream << "reprojection_error_median : " << reprojection_error_magnitudes[reprojection_error_magnitudes.size() / 2] << std::endl;

        //     // TEST: Output all the reprojection errors for external plotting
        //     ofstream error_stream(string(base_path) + "_all_errors.txt", std::ios::out);
        //     for (double error : reprojection_error_magnitudes) {
        //       error_stream << error << std::endl;
        //     }
        //     error_stream.close();
    }
    stream << "reprojection_error_average : " << (reprojection_error_sum / reprojection_error_count) << std::endl;
    stream << "reprojection_error_maximum : " << reprojection_error_max << std::endl;
    stream << "median_kl_divergence : " << biasedness << std::endl;
    stream << "" << std::endl;

    stream << "reprojection_error_histogram_visualization_half_extent_in_pixels : " << histogram_extent_in_px << std::endl;
    stream << "maximum_error_visualization_maximum_error_in_pixels : " << max_error_in_px << std::endl;

    return true;
}


bool CreateCalibrationReportForCamera(
        const char* base_path,
        const int camera_index,
        const Dataset& dataset,
        const BAState& calibration) {
    const CameraModel* cam = calibration.intrinsics[camera_index].get();

    // Create the report directory if it does not exist yet
    QFileInfo(base_path).dir().mkpath(".");

    // Visualize the observation directions
    Image<Vec3u8> model_visualization;
    VisualizeModelDirections(*cam, &model_visualization);
    model_visualization.Write(string(base_path) + "_observation_directions.png");

    // Compute all reprojection errors as input for the corresponding visualizations
    usize reprojection_error_count = 0;
    double reprojection_error_sum = 0.;
    double reprojection_error_max = 0;
    vector<Vec2d> reprojection_errors;
    vector<Vec2f> features;
    ComputeAllReprojectionErrors(
                camera_index, dataset, calibration,
                &reprojection_error_count, &reprojection_error_sum, &reprojection_error_max, &reprojection_errors, &features);

    // Visualize the distribution of reprojection errors
    constexpr int kHistResolution = 50;  // resolution of the visualization image
    constexpr double kHistExtent = 0.2f;  // visualized reprojection error extent in pixels
    Image<double> hist_image;
    ComputeReprojectionErrorHistogram(
                kHistResolution, kHistExtent, reprojection_errors, &hist_image);
    double max_hist_entry = 0;
    for (u32 y = 0; y < hist_image.height(); ++ y) {
        for (u32 x = 0; x < hist_image.width(); ++ x) {
            max_hist_entry = std::max(max_hist_entry, hist_image(x, y));
        }
    }
    Image<u8> hist_image_u8(hist_image.width(), hist_image.height());
    for (u32 y = 0; y < hist_image.height(); ++ y) {
        for (u32 x = 0; x < hist_image.width(); ++ x) {
            hist_image_u8(x, y) = hist_image(x, y) * 255.99f / max_hist_entry;
        }
    }
    hist_image_u8.Write(string(base_path) + "_errors_histogram.png");

    // Create a Voronoi diagram of the feature locations for the following visualizations
    float voronoi_integer_inv_scaling;
    voronoi_diagram<double> vd;
    vector<Vec2f> v_reprojection_errors;
    vector<boost::polygon::point_data<i32>> v_points;
    CreateVoronoiDiagram(cam, reprojection_errors, features, &voronoi_integer_inv_scaling, &vd, &v_reprojection_errors, &v_points);

    // Visualize the reprojection error directions
    Image<Vec3u8> error_direction_image;
    VisualizeReprojectionErrorDirections(
                cam,
                voronoi_integer_inv_scaling,
                vd,
                v_reprojection_errors,
                v_points,
                &error_direction_image);
    error_direction_image.Write(string(base_path) + "_error_directions.png");

    // Visualize the reprojection error magnitudes
    constexpr double max_error_in_px = 0.5;  // maximum error that can be distinguished in this visualization, in pixels
    Image<Vec3u8> error_magnitude_image;
    VisualizeReprojectionErrorMagnitudes(
                cam,
                voronoi_integer_inv_scaling,
                vd,
                v_reprojection_errors,
                v_points,
                max_error_in_px,
                &error_magnitude_image);
    error_magnitude_image.Write(string(base_path) + "_error_magnitudes.png");

    // Compute the biasedness measure
    double biasedness = ComputeBiasedness(cam, reprojection_errors, features);

    // Compute the (approximate) field-of-view
    double horizontal_fov;
    double vertical_fov;
    ComputeApproximateFOV(cam, &horizontal_fov, &vertical_fov);

    // Write the info file
    int num_localized_images = 0;
    for (usize i = 0; i < calibration.image_used.size(); ++ i) {
        if (calibration.image_used[i]) {
            ++ num_localized_images;
        }
    }
    WriteReportInfoFile(
                string(base_path) + "_info.txt",
                cam,
                horizontal_fov,
                vertical_fov,
                dataset.ImagesetCount(),
                num_localized_images,
                reprojection_errors,
                reprojection_error_count,
                reprojection_error_sum,
                reprojection_error_max,
                biasedness,
                kHistExtent,
                max_error_in_px);

    // Additional model-specific visualizations
    if (const CentralGenericModel* cgb_model = dynamic_cast<const CentralGenericModel*>(cam)) {
        // TODO: Implement this in a generic way for all other grid-based cameras as well
        Image<Vec3u8> grid_point_image(cam->width(), cam->height());
        grid_point_image.SetTo(Vec3u8(0, 0, 0));

        for (u32 y = 0; y < cgb_model->grid().height(); ++ y) {
            for (u32 x = 0; x < cgb_model->grid().width(); ++ x) {
                Vec2d pixel = cgb_model->GridPointToPixelCornerConv(x, y);
                int px = pixel.x();
                int py = pixel.y();
                if (px >= 0 && py >= 0 && px < grid_point_image.width() && py < grid_point_image.height()) {
                    grid_point_image(px, py) = Vec3u8(255, 255, 255);
                }
            }
        }

        std::ostringstream filename6;
        filename6 << string(base_path) << "_grid_point_locations.png";
        grid_point_image.Write(filename6.str());
    } else if (const NoncentralGenericModel* ngbsp_model = dynamic_cast<const NoncentralGenericModel*>(cam)) {
        // Find the point that is as close to each of the rays as possible. For
        // central cameras, this will be the optical center.
        CenterPointCostFunction center_point_cost;
        int calibrated_area = (cam->calibration_max_y() - cam->calibration_min_y() + 1) *
                (cam->calibration_max_x() - cam->calibration_min_x() + 1);
        center_point_cost.lines.reserve(calibrated_area);
        center_point_cost.tangents.reserve(calibrated_area);
        for (int y = cam->calibration_min_y(); y <= cam->calibration_max_y(); ++ y) {
            for (int x = cam->calibration_min_x(); x <= cam->calibration_max_x(); ++ x) {
                Line3d line;
                if (ngbsp_model->Unproject(x + 0.5f, y + 0.5f, &line)) {
                    center_point_cost.lines.emplace_back(line);
                    center_point_cost.tangents.emplace_back();
                    ComputeTangentsForDirectionOrLine(line.direction(), &center_point_cost.tangents.back());
                }
            }
        }

        Vec3d center = Vec3d::Zero();
        LMOptimizer<double> optimizer;
        optimizer.Optimize(
                    &center,
                    center_point_cost,
                    /*max_iteration_count*/ 100,
                    /*max_lm_attempts*/ 10,
                    /*init_lambda*/ -1,
                    /*init_lambda_factor*/ 0.001f,
                    /*print_progress*/ false);

        // Get statistics on the distances between the center and the rays
        Image<Vec3d> line_offsets(cam->width(), cam->height());
        line_offsets.SetTo(Vec3d::Constant(numeric_limits<double>::quiet_NaN()));

        double line_distance_sum = 0;
        usize line_distance_count = 0;
        double line_distance_max = 0;
        vector<double> line_distances;

        double max_line_offset_extent = 0;

        for (int y = cam->calibration_min_y(); y <= cam->calibration_max_y(); ++ y) {
            for (int x = cam->calibration_min_x(); x <= cam->calibration_max_x(); ++ x) {
                Line3d line;
                if (ngbsp_model->Unproject(x + 0.5f, y + 0.5f, &line)) {
                    // Find the closest point on the line to the previously determined center point
                    double parameter = line.direction().dot(center - line.origin());
                    Vec3d closest_point_on_line = line.origin() + parameter * line.direction();

                    Vec3d line_offset = closest_point_on_line - center;
                    line_offsets(x, y) = line_offset;

                    double line_distance = line_offset.norm();
                    line_distance_sum += line_distance;
                    ++ line_distance_count;
                    line_distance_max = std::max(line_distance_max, line_distance);
                    line_distances.push_back(line_distance);

                    max_line_offset_extent = std::max(max_line_offset_extent, fabs(line_offset.x()));
                    max_line_offset_extent = std::max(max_line_offset_extent, fabs(line_offset.y()));
                    max_line_offset_extent = std::max(max_line_offset_extent, fabs(line_offset.z()));
                }
            }
        }

        // stream << "line_distance_count : " << line_distance_count << std::endl;
        // if (!line_distances.empty()) {
        //   std::sort(line_distances.begin(), line_distances.end());
        //   stream << "line_distance_median : " << line_distances[line_distances.size() / 2] << std::endl;
        // }
        // stream << "line_distance_average : " << (line_distance_sum / line_distance_count) << std::endl;
        // stream << "line_distance_maximum : " << line_distance_max << std::endl;
        // stream << "" << std::endl;

        Image<Vec3u8> line_offset_visualization(cam->width(), cam->height());
        for (int y = 0; y < cam->height(); ++ y) {
            for (int x = 0; x < cam->width(); ++ x) {
                const Vec3d& line_offset = line_offsets(x, y);
                if (line_offset.hasNaN()) {
                    line_offset_visualization(x, y) = Vec3u8(0, 0, 0);
                } else {
                    line_offset_visualization(x, y) = Vec3u8(
                                127 + 127 * line_offset.x() / max_line_offset_extent,
                                127 + 127 * line_offset.y() / max_line_offset_extent,
                                127 + 127 * line_offset.z() / max_line_offset_extent);
                }
            }
        }

        std::ostringstream filename5;
        filename5 << string(base_path) << "_line_offsets.png";
        line_offset_visualization.Write(filename5.str());

        // Output a 3D model of the camera as an .obj file
        constexpr int kStepForOBJ = 20;
        ofstream obj_stream(string(base_path) + "_line_visualization.obj", std::ios::out);
        ofstream obj_stream_cutoff(string(base_path) + "_line_visualization_cutoff.obj", std::ios::out);
        ofstream obj_stream_origins(string(base_path) + "_line_visualization_origins.obj", std::ios::out);
        if (!obj_stream || !obj_stream_cutoff || !obj_stream_origins) {
            return false;
        }
        obj_stream << std::setprecision(14);
        obj_stream_cutoff << std::setprecision(14);
        obj_stream_origins << std::setprecision(14);

        int visualized_line_count = 0;
        for (int y = cam->calibration_min_y(); y <= cam->calibration_max_y(); y += kStepForOBJ) {
            for (int x = cam->calibration_min_x(); x <= cam->calibration_max_x(); x += kStepForOBJ) {
                Line3d line;
                if (ngbsp_model->Unproject(x + 0.5f, y + 0.5f, &line)) {
                    // Show the range of the line that is centered on the point on the line
                    // that is closest to the center point defined above, and has a certain
                    // extent in both directions.
                    double parameter = line.direction().dot(center - line.origin());
                    Vec3d closest_point_on_line = line.origin() + parameter * line.direction();

                    // Visualized line segments are at least 2 * 10 meters long
                    double visualization_line_half_extent = std::max<double>(10, 10 * (closest_point_on_line - center).norm());

                    Vec3d point_a = closest_point_on_line + visualization_line_half_extent * line.direction();
                    obj_stream << "v " << point_a.x() << " " << point_a.y() << " " << point_a.z() << std::endl;
                    Vec3d point_b = closest_point_on_line - visualization_line_half_extent * line.direction();
                    obj_stream << "v " << point_b.x() << " " << point_b.y() << " " << point_b.z() << std::endl;

                    obj_stream_cutoff << "v " << point_a.x() << " " << point_a.y() << " " << point_a.z() << std::endl;
                    obj_stream_cutoff << "v " << closest_point_on_line.x() << " " << closest_point_on_line.y() << " " << closest_point_on_line.z() << std::endl;

                    obj_stream_origins << "v " << point_a.x() << " " << point_a.y() << " " << point_a.z() << std::endl;
                    obj_stream_origins << "v " << point_b.x() << " " << point_b.y() << " " << point_b.z() << std::endl;
                    obj_stream_origins << "v " << line.origin().x() << " " << line.origin().y() << " " << line.origin().z() << std::endl;

                    ++ visualized_line_count;
                }
            }
        }

        int vertex_index = 1;  // vertex indexing starts at 1 in .obj files
        for (int i = 0; i < visualized_line_count; ++ i) {
            obj_stream << "l " << vertex_index << " " << (vertex_index + 1) << std::endl;
            obj_stream_cutoff << "l " << vertex_index << " " << (vertex_index + 1) << std::endl;
            obj_stream_origins << "l " << (3 * (vertex_index / 2) + 1) << " " << (3 * (vertex_index / 2) + 2) << std::endl;
            vertex_index += 2;
        }
    }

    return true;
}

void CreateReprojectionErrorHistogram(
        int camera_index,
        const Dataset& dataset,
        const BAState& state,
        Image<u8>* histogram_image) {
    constexpr int kHistResolution = 50;  // resolution of the visualization image
    constexpr double kHistExtent = 0.2f;  // visualized reprojection error extent in pixels

    Image<double> hist_image(kHistResolution, kHistResolution);
    hist_image.SetTo(0.f);

    for (int imageset_index = 0; imageset_index < dataset.ImagesetCount(); ++ imageset_index) {
        if (!state.image_used[imageset_index]) {
            continue;
        }

        const SE3d& image_tr_global = state.camera_tr_rig[camera_index] * state.rig_tr_global[imageset_index];
        Mat3d image_r_global = image_tr_global.rotationMatrix();
        const Vec3d& image_t_global = image_tr_global.translation();

        shared_ptr<const Imageset> imageset = dataset.GetImageset(imageset_index);
        const vector<PointFeature>& features = imageset->FeaturesOfCamera(camera_index);

        for (const PointFeature& feature : features) {
            if (feature.index < 0 || feature.index >= state.points.size()) {
                LOG(ERROR) << "Encountered an invalid value of feature.index.";
                histogram_image->SetSize(hist_image.width(), hist_image.height());
                histogram_image->SetTo(static_cast<u8>(0));
                return;
            }

            Vec3d local_point = image_r_global * state.points[feature.index] + image_t_global;
            Vec2d pixel;
            if (state.intrinsics[camera_index]->Project(local_point, &pixel)) {
                Vec2d reprojection_error = pixel - feature.xy.cast<double>();

                double hx_f = kHistResolution * 0.5f * ((reprojection_error.x() / kHistExtent) + 1.f);
                int hx = static_cast<int>(hx_f) - ((hx_f < 0) ? 1.f : 0.f);
                double hy_f = kHistResolution * 0.5f * ((reprojection_error.y() / kHistExtent) + 1.f);
                int hy = static_cast<int>(hy_f) - ((hy_f < 0) ? 1.f : 0.f);

                if (hx >= 0 && hy >= 0 &&
                        hx < static_cast<int>(hist_image.width()) &&
                        hy < static_cast<int>(hist_image.height())) {
                    hist_image(hx, hy) += 1.f;
                }
            }
        }
    }

    double max_hist_entry = 0;
    for (u32 y = 0; y < hist_image.height(); ++ y) {
        for (u32 x = 0; x < hist_image.width(); ++ x) {
            max_hist_entry = std::max(max_hist_entry, hist_image(x, y));
        }
    }

    histogram_image->SetSize(hist_image.width(), hist_image.height());
    for (u32 y = 0; y < hist_image.height(); ++ y) {
        for (u32 x = 0; x < hist_image.width(); ++ x) {
            (*histogram_image)(x, y) = hist_image(x, y) * 255.99f / max_hist_entry;
        }
    }
}

void CreateReprojectionErrorDirectionVisualization(
        const Dataset& dataset,
        const int camera_index,
        const BAState& calibration,
        Image<Vec3u8>* visualization) {
    // Show all reprojection errors of a camera in one image, colored by direction (as in optical flow),
    // to see whether in certain areas all errors have the same direction.
    visualization->SetSize(calibration.intrinsics[camera_index]->width(), calibration.intrinsics[camera_index]->height());
    visualization->SetTo(Vec3u8(0, 0, 0));

    for (int imageset_index = 0; imageset_index < dataset.ImagesetCount(); ++ imageset_index) {
        if (!calibration.image_used[imageset_index]) {
            continue;
        }

        const SE3d& image_tr_global = calibration.image_tr_global(camera_index, imageset_index);
        Mat3d image_r_global = image_tr_global.rotationMatrix();
        const Vec3d& image_t_global = image_tr_global.translation();

        shared_ptr<const Imageset> imageset = dataset.GetImageset(imageset_index);
        const vector<PointFeature>& features = imageset->FeaturesOfCamera(camera_index);

        for (const PointFeature& feature : features) {
            Vec3d local_point = image_r_global * calibration.points[feature.index] + image_t_global;
            Vec2d pixel;
            if (calibration.intrinsics[camera_index]->Project(local_point, &pixel)) {
                Vec2d reprojection_error = pixel - feature.xy.cast<double>();

                int fx = feature.xy.x();
                int fy = feature.xy.y();

                if (fx >= 0 && fy >= 0 &&
                        fx < static_cast<int>(visualization->width()) &&
                        fy < static_cast<int>(visualization->height())) {
                    double dir = atan2(reprojection_error.y(), reprojection_error.x());  // from -M_PI to M_PI
                    (*visualization)(fx, fy) = Vec3u8(
                                127 + 127 * sin(dir),
                                127 + 127 * cos(dir),
                                127);
                }
            }
        }
    }
}

void CreateReprojectionErrorMagnitudeVisualization(
        const Dataset& dataset,
        int camera_index,
        const BAState& calibration,
        double max_error,
        Image<Vec3u8>* visualization) {
    Image<double> max_error_image(calibration.intrinsics[camera_index]->width(), calibration.intrinsics[camera_index]->height());
    max_error_image.SetTo(-1.f);

    for (int imageset_index = 0; imageset_index < dataset.ImagesetCount(); ++ imageset_index) {
        if (!calibration.image_used[imageset_index]) {
            continue;
        }

        const SE3d& image_tr_global = calibration.image_tr_global(camera_index, imageset_index);
        Mat3d image_r_global = image_tr_global.rotationMatrix();
        const Vec3d& image_t_global = image_tr_global.translation();

        shared_ptr<const Imageset> imageset = dataset.GetImageset(imageset_index);
        const vector<PointFeature>& features = imageset->FeaturesOfCamera(camera_index);

        for (const PointFeature& feature : features) {
            Vec3d local_point = image_r_global * calibration.points[feature.index] + image_t_global;
            Vec2d pixel;
            if (calibration.intrinsics[camera_index]->Project(local_point, &pixel)) {
                Vec2d reprojection_error = pixel - feature.xy.cast<double>();

                int fx = feature.xy.x();
                int fy = feature.xy.y();

                if (fx >= 0 && fy >= 0 &&
                        fx < static_cast<int>(max_error_image.width()) &&
                        fy < static_cast<int>(max_error_image.height())) {
                    // Update max error magnitude visualization
                    max_error_image(fx, fy) = std::max(max_error_image(fx, fy), reprojection_error.norm());
                }
            }
        }
    }

    double max_error_in_image = 0.f;
    for (u32 y = 0; y < max_error_image.height(); ++ y) {
        for (u32 x = 0; x < max_error_image.width(); ++ x) {
            if (max_error_image(x, y) >= 0) {
                max_error_in_image = std::max(max_error_in_image, max_error_image(x, y));
            }
        }
    }
    LOG(1) << "Maximum reprojection error: " << max_error_in_image;

    if (max_error <= 0) {
        max_error = max_error_in_image;
    }

    visualization->SetSize(max_error_image.width(), max_error_image.height());
    for (u32 y = 0; y < max_error_image.height(); ++ y) {
        for (u32 x = 0; x < max_error_image.width(); ++ x) {
            if (max_error_image(x, y) >= 0) {
                double factor = std::min(1., max_error_image(x, y) / max_error);
                visualization->at(x, y) = Vec3u8(255.99f * factor, 255.99f * (1 - factor), 0);
            } else {
                // No feature detection at this pixel.
                visualization->at(x, y) = Vec3u8(0, 0, 0);
            }
        }
    }
}

void VisualizeModelDirections(
        const CameraModel& model,
        Image<Vec3u8>* visualization) {
    visualization->SetSize(model.width(), model.height());

    Image<Vec3d> directions;
    if (!CreateObservationDirectionsImage(&model, &directions)) {
        LOG(ERROR) << "CreateObservationDirectionsImage() failed, no visualization will be provided.";
        visualization->SetTo(Vec3u8::Zero());
        return;
    }

    for (u32 y = 0; y < visualization->height(); ++ y) {
        for (u32 x = 0; x < visualization->width(); ++ x) {
            const Vec3d& direction = directions(x, y);
            if (direction.hasNaN()) {
                visualization->at(x, y) = Vec3u8(0, 0, 0);
            } else {
                visualization->at(x, y) = Vec3u8(
                            70 * 255.99f / 2.f * (direction.x() + 1),
                            70 * 255.99f / 2.f * (direction.y() + 1),
                            270 * 255.99f / 2.f * (direction.z() + 1));
            }
        }
    }
}

}

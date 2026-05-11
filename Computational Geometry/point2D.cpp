/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * 2D point (vector) struct with standard geometric operations
 *
 * EPS is used for floating-point equality comparisons.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <stdexcept>
#include <vector>

namespace Geometry {

    constexpr double EPS = 1e-9;

    /*
     * Classification of a point's position relative to a polygon
     */
    enum class PointStatus {
        IN,
        ON,
        OUT
    };

    /*
     * 2D point or vector with coordinates of type T
     */
    template <typename T>
    struct point2D {
        T x, y;

        point2D() : x(0), y(0) {}

        point2D(T x, T y) : x(x), y(y) {}

        /*
         * Checks equality within EPS tolerance for floating-point T, exact otherwise.
         */
        bool operator== (const point2D &o) const {
            if constexpr (std::is_floating_point_v<T>) {
                return std::fabs(x - o.x) < EPS && std::fabs(y - o.y) < EPS;
            } else {
                return x == o.x && y == o.y;
            }
        }

        /* @return  Component-wise sum */
        point2D operator+ (const point2D &o) const {
            return point2D(x + o.x, y + o.y);
        }

        /* @return  Component-wise difference */
        point2D operator- (const point2D &o) const {
            return point2D(x - o.x, y - o.y);
        }

        /*
         * @param s  Scalar multiplier
         * @return   This vector scaled by s; result type is decltype(T * S)
         */
        template <typename S>
        auto operator* (S s) const {
            using R = decltype(T{} * S{});
            return point2D<R>(x * s, y * s);
        }

        /*
         * @param s  Scalar divisor. throws std::invalid_argument on zero
         * @return   Vector divided by s
         */
        template <typename S>
        auto operator/ (S s) const {
            using R = decltype(T{} / S{});

            if constexpr (std::is_floating_point_v<S>) {
                if (std::fabs(s) < EPS) {
                    throw std::invalid_argument("Division by zero in point2D::operator/");
                }
            } else {
                if (s == 0) {
                    throw std::invalid_argument("Division by zero in point2D::operator/");
                }
            }

            return point2D<R>(x / s, y / s);
        }

        /*
         * @return  Positive if angle between vectors < 90°, zero if perpendicular, negative if > 90°
         */
        T dot(const point2D &o) const {
            return x * o.x + y * o.y;
        }

        /*
         * 2D cross product: x*o.y - y*o.x
         *
         * @return  Positive if o is CCW from this, negative if CW, zero if collinear
         */
        T cross(const point2D &o) const {
            return x * o.y - y * o.x;
        }

        /* @return  Euclidean distance to o */
        double distance(const point2D &o) const {
            point2D d = *this - o;
            return std::hypot(d.x, d.y);
        }

        /*
         * Signed angle from this vector to o, in radians, range [-pi, pi].
         * Positive means o is CCW from this.
         */
        double angle(const point2D &o) const {
            return std::atan2(cross(o), dot(o));
        }
    };

    /* Comparator: ascending x, tiebreak ascending y */
    template <typename T>
    inline bool byX(const point2D<T>& a, const point2D<T>& b) {
        if (std::fabs(a.x - b.x) >= EPS)
            return a.x < b.x;
        return a.y < b.y;
    }

    /* Comparator: ascending y, tiebreak ascending x */
    template <typename T>
    inline bool byY(const point2D<T>& a, const point2D<T>& b) {
        if (std::fabs(a.y - b.y) >= EPS)
            return a.y < b.y;
        return a.x < b.x;
    }

    /*
     * Comparator: CCW polar angle order around p0.
     * Collinear points are sorted closer to p0 first.
     */
    template <typename T>
    inline bool byAngle(const point2D<T>& a, const point2D<T>& b, const point2D<T>& p0) {
        T cr = (a - p0).cross(b - p0);
        if (std::fabs(cr) >= EPS)
            return cr > 0;
        return p0.distance(a) < p0.distance(b);
    }

    /*
     * Moves the bottommost-leftmost point in pts to pts[0], then sorts
     * the remaining points in CCW polar angle order around it.
     *
     * @pre  pts.size() >= 2
     * Time complexity: O(n log n)
     */
    template <typename T>
    inline void sortByAngle(std::vector<point2D<T>>& pts) {
        auto it = std::min_element(pts.begin(), pts.end(), byY<T>);
        std::swap(*it, pts[0]);
        point2D<T> p0 = pts[0];
        std::sort(pts.begin() + 1, pts.end(),
            [&](const point2D<T>& a, const point2D<T>& b) {
                return byAngle(a, b, p0);
            });
    }

    /*
     * @return  +1 if a → b → c turns CCW, -1 if CW, 0 if collinear
     */
    template <typename T>
    inline int orientation(const point2D<T>& a, const point2D<T>& b, const point2D<T>& c) {
        T v = (b - a).cross(c - a);
        if (v >= EPS)
            return 1;
        if (v <= -EPS)
            return -1;
        return 0;
    }
}

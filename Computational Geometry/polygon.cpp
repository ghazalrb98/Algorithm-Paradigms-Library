/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * Polygon struct supporting area computation and point containment tests
 *
 */

#pragma once

#include "point2D.cpp"

#include <algorithm>
#include <vector>

namespace Geometry {

    /*
     * Simple polygon defined by a list of corner points
     *
     * @pre  corners must be ordered consecutively around the polygon boundary
     *       (either all CCW or all CW); use sortByAngle(corners) to obtain
     *       CCW order before constructing this polygon.
     */
    template <typename T>
    struct polygon {
        std::vector<point2D<T>> corners;
        polygon(std::vector<point2D<T>> corners) : corners(corners) {}

        /*
         * Returns twice the signed area using the shoelace formula
         *
         * Positive for CCW ordering, negative for CW ordering.
         * Returning 2× avoids a division and keeps the result exact for integer coordinates.
         *
         * @return  2 * signed area
         *
         * Time complexity: O(n)
         */
        T SignedDoubledArea() const {
            const size_t n = corners.size();
            T sum = 0;
            for (size_t i = 0; i < n; i++) {
                sum += corners[i].cross(corners[(i + 1) % n]);
            }
            return sum;
        }

        /*
         * Returns twice the absolute area
         *
         * @return  2 * area  (always non-negative)
         *
         * Time complexity: O(n)
         */
        T DoubledArea() const {
            return std::fabs(SignedDoubledArea());
        }

        /*
         * Returns the signed area (positive for CCW, negative for CW)
         *
         * @return  Signed area
         *
         * Time complexity: O(n)
         */
        double SignedArea() const {
            return SignedDoubledArea() / 2.0;
        }

        /*
         * Returns the area (always non-negative)
         *
         * @return  Area
         *
         * Time complexity: O(n)
         */
        double Area() const {
            return DoubledArea() / 2.0;
        }

        /*
         * @return  True if p is collinear with and between a and b (inclusive endpoints)
         */
        static bool onSegment(point2D<T> a, point2D<T> b, point2D<T> p) {
            if (orientation(a, b, p) != 0) 
                return false;
            return (p - a).dot(p - b) <= 0;
        }

        /*
         * Classifies p relative to the polygon.
         *
         * Checks each edge for boundary membership first, then uses the angle-sum (winding) method:
         * angles from p to each edge sum to ±2π if inside, 0 if outside.
         *
         * @param p  Query point
         * @return   IN, ON, or OUT
         *
         * Time complexity: O(n)
         */
        PointStatus status(const point2D<T> &p) const {
            const size_t n = corners.size();
            for (size_t i = 0; i < n; i++) {
                if (onSegment(corners[i], corners[(i + 1) % n], p))
                    return PointStatus::ON;
            }
            double a = 0;
            for (size_t i = 0; i < n; i++) {
                point2D<T> s = corners[i] - p;
                point2D<T> d = corners[(i + 1) % n] - p;
                a += s.angle(d);
            }
            if (std::fabs(a) < EPS)
                return PointStatus::OUT;
            return PointStatus::IN;
        }
    };
}

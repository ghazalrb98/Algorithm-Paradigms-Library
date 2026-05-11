/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * Convex hull of a 2D point set using Graham scan
 *
 */

#pragma once

#include "point2D.cpp"

#include <vector>
#include <algorithm>

namespace Geometry {

    /*
     * Computes the convex hull of a 2D point set
     */
    template<typename T>
    struct ConvexHull {

        /*
         * Returns the convex hull of points in CCW order
         *
         * Uses Graham scan. Duplicate points are removed before processing.
         * If all points are collinear, returns the two extreme endpoints.
         *
         * @param points  Input point set
         * @return  Hull vertices in CCW order, or two endpoints if all collinear
         *
         * Time complexity: O(n log n)
         */
        static std::vector<point2D<T>> get(std::vector<point2D<T>> points) {
            std::sort(points.begin(), points.end(), byX<T>);
            points.erase(std::unique(points.begin(), points.end()), points.end());
            if (points.size() <= 1)
                return points;

            std::vector<point2D<T>> hull;
            sortByAngle(points);
            point2D<T> p0 = points[0];

            hull.push_back(p0);
            for (int i = 1; i < points.size(); i++) {
                while (hull.size() >= 2 &&
                    orientation(hull[hull.size() - 2],
                                hull.back(),
                                points[i]) <= 0) {
                    hull.pop_back();
                }
                hull.push_back(points[i]);
            }

            bool allCollinear = true;
            for (int i = 2; i < points.size(); i++) {
                if (orientation(points[0], points[1], points[i]) != 0) {
                    allCollinear = false;
                    break;
                }
            }
            if (allCollinear) {
                std::sort(points.begin(), points.end(), byX<T>);
                hull.clear();
                hull.push_back(points.front());
                hull.push_back(points.back());
            }
            return hull;
        }
    };
}

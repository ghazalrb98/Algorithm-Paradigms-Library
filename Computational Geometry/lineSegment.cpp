/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * Line segment utilities for computational geometry
 *
 */

#include "point2D.cpp"
#include <vector>
#include <algorithm>

namespace Geometry {
    /*
     * Closed line segment defined by two endpoints a and b
     *
     * All intersection, containment, and distance queries use floating-point
     * tolerance EPS to handle numerical imprecision.
     */
    struct lineSegment {
        point2D<double> a, b;
        lineSegment(const point2D<double>& a, const point2D<double>& b) : a(a), b(b) {}

        /*
         * Returns true if the segment degenerates to a single point (a == b)
         */
        bool isPoint() const {
            return a == b;
        }

        /*
         * Returns true if the given point lies on this segment
         *
         * @param point  The query point
         *
         * Time complexity: O(1)
         */
        bool hasPoint(const point2D<double>& point) const {
            const point2D<double> d = b - a;
            const point2D<double> p = point - a;
            if (std::fabs(d.cross(p)) > EPS)
                return false;
            const double t = p.dot(d) / d.dot(d);
            return t >= -EPS && t <= 1 + EPS;
        }

        /*
         * Returns the minimum distance from point to this segment
         *
         * @param point  The query point
         * Time complexity: O(1)
         */
        double distance(const point2D<double>& point) const {
            const point2D<double> dp = b - a;
            if ((point - a).dot(dp) >= -EPS && (point - b).dot(dp) <= EPS)
                return std::fabs((point - a).cross(dp)) / a.distance(b);
            return std::min(point.distance(a), point.distance(b));
        }

        /*
         * Returns the intersection of this segment with segment o
         *
         * Handles all degenerate cases (point-point, point-segment, collinear overlap).
         *
         * Return values:
         *   - Empty vector : no intersection
         *   - One point    : segments touch at a single point
         *   - Two points   : segments overlap along a collinear sub-segment
         *                    (endpoints returned in lexicographic order)
         *
         * @param o  The other line segment
         *
         * Time complexity: O(1)
         */
        std::vector<point2D<double>> intersect(const lineSegment &o) const {
            std::vector<point2D<double>> result;

            const bool A = isPoint();
            const bool B = o.isPoint();
            if (A && B) {
                if (a == o.a)
                    result.push_back(a);
                return result;
            }
            if (A) {
                if (o.hasPoint(a))
                    result.push_back(a);
                return result;
            }
            if (B) {
                if (hasPoint(o.a))
                    result.push_back(o.a);
                return result;
            }

            const point2D<double> p = a, q = o.a, dp = b - a, dq = o.b - o.a;
            const double PxQ = dp.cross(dq);

            if (std::fabs(PxQ) < EPS) {
                if (std::fabs((q - p).cross(dp)) > EPS)
                    return result;

                double alpha = (q - p).dot(dp) / dp.dot(dp);
                double beta = (o.b - p).dot(dp) / dp.dot(dp);
                if (alpha > beta)
                    std::swap(alpha,beta);

                const double s = std::max(0.0, alpha);
                const double t = std::min(1.0, beta);
                if (s > t + EPS)
                    return result;

                point2D<double> p1 = p + dp * s;
                point2D<double> p2 = p + dp * t;

                if (p1 == p2)
                    result.push_back(p1);
                else {
                    if (p2.x < p1.x || (std::fabs(p1.x - p2.x) < EPS && p2.y < p1.y))
                        std::swap(p1,p2);
                    result.push_back(p1);
                    result.push_back(p2);
                }
                return result;
            }

            const double s = (q - p).cross(dq) / PxQ;
            const double t = (q - p).cross(dp) / PxQ;
            if (s >= -EPS && s <= 1 + EPS && t >= -EPS && t <= 1 + EPS)
                result.push_back(p + dp * s);
            return result;
        }

        /*
         * Returns the minimum distance between this segment and segment o
         *
         * Handles all degenerate cases (point-point, point-segment).
         * For two proper segments, returns 0 if they intersect; otherwise
         * returns the minimum of the four endpoint-to-segment distances.
         *
         * @param o  The other line segment
         *
         * Time complexity: O(1)
         */
        double distance(const lineSegment& o) const {
            const bool A = isPoint();
            const bool B = o.isPoint();
            if (A && B)
                return a.distance(o.a);
            if (A)
                return o.distance(a);
            if (B)
                return distance(o.a);

            if (intersect(o).empty())
                return std::min(distance(o.a),
                       std::min(distance(o.b),
                       std::min(o.distance(a), o.distance(b))));
            return 0;
        }
    };
}

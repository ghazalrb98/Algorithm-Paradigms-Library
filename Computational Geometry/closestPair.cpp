/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * Closest pair of points using divide and conquer in O(n log n)
 * 
 */

#include "point2D.cpp"
#include <vector>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace Geometry {
    /*
     * Holds the result of a closest pair query
     */
    struct closestPairResult {
        double dis;
        point2D<double> p1, p2;

        closestPairResult()
            : dis(std::numeric_limits<double>::infinity()),
            p1(0.0, 0.0),
            p2(0.0, 0.0) {}
    };

    struct closestPair {
        std::vector<point2D<double> > points, tmp;

        /*
         * Sorts input points by x-coordinate and prepares a buffer
         *
         * @param points  Input set of 2D points
         *
         * Time complexity: O(n log n)
         * Space complexity: O(n)
         */
        closestPair(const std::vector<point2D<double> > &points) : points(points) {
            std::sort(this->points.begin(), this->points.end(), byX<double>);
            tmp.assign(points.size(), {0.0, 0.0});
        }

        /*
         * Divide and conquer step over points[l..r] (inclusive)
         *
         * Splits at mid, recurses on both halves, then scans the strip
         * [mid.x - d, mid.x + d] sorted by y to find cross-half improvements.
         *
         * @param l  Left index (inclusive)
         * @param r  Right index (inclusive)
         * @return   Closest pair result within points[l..r]
         *
         * Time complexity: O(n log n) total across all recursive calls
         * Space complexity: O(n)
         */
        closestPairResult find(const size_t l, const size_t r) {
            if (r - l <= 2) {
                return bruteForce(l, r);
            }

            closestPairResult ans;
            const size_t mid = (l + r) / 2;
            const closestPairResult left = find(l, mid);
            const closestPairResult right = find(mid + 1, r);
            if (left.dis < right.dis) {
                ans.dis = left.dis;
                ans.p1 = left.p1;
                ans.p2 = left.p2;
            } else {
                ans.dis = right.dis;
                ans.p1 = right.p1;
                ans.p2 = right.p2;
            }
            const auto lPtr = std::lower_bound(points.begin() + l, points.begin() + r + 1,
                                            point2D<double>(points[mid].x - ans.dis, 0.0), byX<double>);
            const auto rPtr = std::upper_bound(points.begin() + l, points.begin() + r + 1,
                                            point2D<double>(points[mid + 1].x + ans.dis, 0.0), byX<double>);
            const size_t size = rPtr - lPtr;
            std::copy(lPtr, rPtr, tmp.begin());
            std::sort(tmp.begin(), tmp.begin() + size, byY<double>);
            for (size_t i = 0; i < size; i++)
                for (size_t j = i + 1; j < std::min(i + 16, size); j++) {
                    const double dis = tmp[i].distance(tmp[j]);
                    if (dis < ans.dis) {
                        ans.dis = dis;
                        ans.p1 = tmp[i];
                        ans.p2 = tmp[j];
                    }
                }

            return ans;
        }

        /*
         * search over points[l..r] (inclusive), used for base cases of size ≤ 3
         *
         * @param l  Left index (inclusive)
         * @param r  Right index (inclusive)
         * @return   Closest pair result within points[l..r]
         *
         * Time complexity: O((r - l)²)
         * Space complexity: O(1)
         */
        closestPairResult bruteForce(size_t l, size_t r) {
            closestPairResult ans;

            for (size_t i = l; i <= r; ++i) {
                for (size_t j = i + 1; j <= r; ++j) {
                    double dis = points[i].distance(points[j]);

                    if (dis < ans.dis) {
                        ans.dis = dis;
                        ans.p1 = points[i];
                        ans.p2 = points[j];
                    }
                }
            }

            return ans;
        }

        /*
         * Entry point: returns the closest pair among all input points
         *
         * @throws std::invalid_argument  If fewer than two points were provided
         * @return   Closest pair result across the entire point set
         *
         * Time complexity: O(n log n)
         * Space complexity: O(n)
         */
        closestPairResult find() {
            if (points.size() < 2) {
                throw std::invalid_argument("Need at least two points to find closest pair");
            }

            return find(0, points.size() - 1);
        }
    };
}

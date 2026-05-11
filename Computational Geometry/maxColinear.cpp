/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * Largest collinear subset of a point set
 * 
 */

#include "point2D.cpp"
#include "fraction.cpp"
#include <vector>
#include <map>

namespace Geometry {

    /*
     * Finds the largest subset of collinear points.
     *
     * Grouping all other points by their exact rational slope 
     * (represented as a reduced fraction dy/dx). 
     *
     * @pre T must be an integer type — fraction<T> relies on exact GCD arithmetic.
     */
    template <typename T>
    struct maxColinear {

        /*
         * @param points  Input set of 2D points
         * @return        Largest collinear subset; if multiple subsets tie, returns
         *                the first one found in input order. Returns the input unchanged
         *                if fewer than 2 points are given.
         *
         * Time complexity:  O(n² log n)
         * Space complexity: O(n)
         */
        static std::vector<point2D<T>> get(const std::vector<point2D<T>>& points) {
            if (points.size() < 2)
                return points;
            size_t best = 0;
            std::vector<point2D<T>> result;
            for (point2D<T> fixed : points) {
                std::map<fraction<T>, std::vector<point2D<T>>> colinear;
                for (point2D<T> relative : points) {
                    if (relative == fixed)
                        continue;
                    fraction<T> slope(relative.y - fixed.y, relative.x - fixed.x);
                    colinear[slope].push_back(relative);
                }
                for (auto it : colinear) {
                    if (it.second.size() > best) {
                        best = it.second.size();
                        result = it.second;
                        result.push_back(fixed);
                    }
                }
            }
            return result;
        }
    };
}

/*
* Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 *  This file implements fast computation of linear recurrences
 *  using matrix exponentiation.
 *
 * Supported recurrence:
 *      x_t = a0 + a1*x_{t-1} + a2*x_{t-2} + ... + aN*x_{t-N}
 *
 * Input representation:
 *   - coefficients = [a0, a1, ..., aN]
 *   - initialValues = [x0, x1, ..., x_{N-1}]
 *
 * Mathematical idea:
 *  The recurrence is converted into a matrix transition:
 *      state_{t+1} = A * state_t
 *
 *  using the augmented state:
 *      [x_t, x_{t-1}, ..., x_{t-N+1}, 1]^T
 *
 *  Then:
 *      state_t = A^(t-(N-1)) * state_{N-1}
 */

#include "LinearAlgebra.cpp"

template <typename T>
class LinearRecurrenceSolver {
    using Matrix = LinearAlgebra::Matrix<T>;

    std::vector<T> coefficients;
    std::vector<T> initialValues;

public:
    LinearRecurrenceSolver(
        const std::vector<T>& recurrenceCoefficients,
        const std::vector<T>& recurrenceInitialValues
    )
        : coefficients(recurrenceCoefficients), initialValues(recurrenceInitialValues) {
    }

    /*
     * Computes x_t % mod for the recurrence:
     *      x_t = a0 + a1*x_{t-1} + ... + aN*x_{t-N}
     *
     * Returns:
     *  - x_t reduced modulo mod
     *
     * Time complexity:
     *  - O(N^3 * log t)
     *
     * Space Complexity: O(N^2)
     *
     * Notes:
     *  - If t < N, the answer is one of the given initial values.
     *  - Negative coefficients and initial values are normalized
     *    into the interval [0, mod).
     */
    T solve(long long int t, T mod) const {
        using NumberTheory::ModularArithmetic;

        const int n = static_cast<int>(initialValues.size());

        if (t < n) {
            return ModularArithmetic<T>(mod).normalize(initialValues[t]);
        }

        LinearAlgebra::MatrixArithmetic<T> matrix_arithmetic(mod);

        Matrix transition(n + 1, std::vector<T>(n + 1, 0));

        /*
         * First row computes:
         *      x_t = a0 + a1*x_{t-1} + ... + aN*x_{t-N}
         *
         * In matrix form with state:
         *      [x_{current}, x_{current-1}, ..., 1]^T
         *
         * the first row is:
         *      [a1, a2, ..., aN, a0]
         */
        for (int i = 1; i <= n; i++) {
            transition[0][i - 1] = ModularArithmetic<T>(mod).normalize(coefficients[i]);
        }
        transition[0][n] = ModularArithmetic<T>(mod).normalize(coefficients[0]);

        /*
         * Shift rows:
         * move each previous value one position down
         */
        for (int i = 1; i < n; i++) {
            transition[i][i - 1] = 1;
        }

        /*
         * Keep the augmented constant 1 alive
         */
        transition[n][n] = 1;

        /*
         * Base state:
         *      [x_{n-1}, x_{n-2}, ..., x_0, 1]^T
         */
        Matrix state(n + 1, std::vector<T>(1, 0));
        for (int i = 0; i < n; i++) {
            state[i][0] = ModularArithmetic<T>(mod).normalize(initialValues[n - i - 1]);
        }
        state[n][0] = ModularArithmetic<T>(mod).normalize(1);

        /*
         * start from state_{n-1}, so to reach state_t
         * we need:
         *      transition^(t - (n - 1))
         */
        Matrix transitionPower = matrix_arithmetic.power(transition, t - (n - 1));
        Matrix result = matrix_arithmetic.multiply(transitionPower, state);
        return result[0][0];
    }
};
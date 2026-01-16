#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>
#include <execution>
#include <chrono>

using namespace std;

class GBM_MonteCarlo_Pricer
{
    private:
        double S0, K, r, T, sigma;
        size_t num_paths;
        bool use_antithetic;
        size_t seed;

        double drift, diffusion, discount_factor;

    public:
        GBM_MonteCarlo_Pricer(double S0_d, double K_d, double r_d, double T_d, double sigma_d,
                             size_t num_paths_d, bool use_antithetic_d = true, size_t seed_d = random_device{}())
            : S0(S0_d), K(K_d), r(r_d), T(T_d), sigma(sigma_d), num_paths(num_paths_d),
              use_antithetic(use_antithetic_d), seed(seed_d)
        {
            if(S0 <= 0.0 || sigma <= 0.0 || T <= 0.0)
            {
                throw invalid_argument("Invalid Parameters for GBM.");
            }
        }
};

int main()
{
    
}
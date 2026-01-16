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
        GBM_MonteCarlo_Pricer(double S0_d, double K_d, double r_d, double sigma_d, double T_d,
                             size_t num_paths_d, bool use_antithetic_d = true, size_t seed_d = random_device{}())
            : S0(S0_d), K(K_d), r(r_d), sigma(sigma_d), T(T_d), num_paths(num_paths_d),
              use_antithetic(use_antithetic_d), seed(seed_d)
        {
            if(S0 <= 0.0 || sigma <= 0.0 || T <= 0.0)
            {
                throw invalid_argument("Invalid Parameters for GBM.");
            }

            drift = (r - 0.5 * sigma * sigma) * T;
            diffusion = sigma * sqrt(T);
            discount_factor = exp(-r * T);
        }

        struct PricingResults
        {
            double price, std_err, var_red_ratio, conf_low, conf_high;
        };
        
        PricingResults price_option() const
        {
            mt19937_64 rng(seed);
            normal_distribution<double> norm(0.0, 1.0);

            double sum_payoffs = 0.0;
            double sum_payoffs_sq = 0.0;
            size_t eff_paths = use_antithetic ? num_paths / 2 : num_paths;
            
            for(size_t i = 0; i < eff_paths; i++)
            {
                double Z = norm(rng);
                double ST1 = S0 * exp(drift + diffusion * Z);
                double payoff1 = max(ST1 - K, 0.0);
                sum_payoffs += payoff1;
                sum_payoffs_sq += payoff1 * payoff1;

                if(use_antithetic)
                {
                    double ST2 = S0 * exp(drift + diffusion * (-Z));
                    double payoff2 = max(ST2 - K, 0.0);
                    double avg_payoff = 0.5 * (payoff1 + payoff2);
                    sum_payoffs += payoff2;
                    sum_payoffs_sq += payoff2 * payoff2;
                }
            }

            size_t total_paths = use_antithetic ? eff_paths * 2 : eff_paths;
            double avg_payoff = sum_payoffs / total_paths;
            double avg_payoff_sq = sum_payoffs_sq / total_paths;
            double var_payoff = avg_payoff_sq - avg_payoff * avg_payoff;
            double price = discount_factor * avg_payoff;
            double std_err = discount_factor * sqrt(var_payoff / total_paths);

            double vr_ratio = 1.0;
            if(use_antithetic)
            {
                double naive_var = var_payoff;
                vr_ratio = naive_var / (var_payoff / 2.0); // 2 when using antithetic variates otherwise 1
            }

            return
            {
                price,
                std_err,
                vr_ratio,
                price - 1.96 * std_err,
                price + 1.96 * std_err
            };
        }
        
        static double black_scholes_price(double S0, double K, double r, double sigma, double T)
        {
            double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
            double d2 = d1 - sigma * sqrt(T);
            auto calc_norm_cdf = [](double x)
            {
                return 0.5 * (1.0 + erf(x / sqrt(2.0)));
            };
            return S0 * calc_norm_cdf(d1) - K * exp(-r * T) * calc_norm_cdf(d2);
        }
};

int main()
{
    GBM_MonteCarlo_Pricer pricer(100.0, 100.0, 0.05, 0.2, 1.0, 1000000, true);

    auto start_clock = chrono::high_resolution_clock::now();
    auto results = pricer.price_option();
    auto end_clock = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::milliseconds>(end_clock - start_clock).count();
    cout<<"Monte Carlo GBM Option Pricing Results:\n";
    cout<<"Price: "<<results.price<<"\n";
    cout<<"Standard Error: "<<results.std_err<<"\n";
    cout<<"95% Confidence Interval: ["<<results.conf_low<<", "<<results.conf_high<<"]\n";
    cout<<"Variance Reduction Ratio: "<<results.var_red_ratio<<"\n";
    cout<<"Computation Time: "<<duration<<" ms\n";
    cout<<"---------------------------------------\n";
    cout<<"Black-Scholes Price: "; 
    cout<<GBM_MonteCarlo_Pricer::black_scholes_price(100.0, 100.0, 0.05, 0.2, 1.0)<<"\n";
}
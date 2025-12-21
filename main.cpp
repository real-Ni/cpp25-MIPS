#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>

#include "abp_sim.h"
#include "cluster_analysis.h"

void run_simulation(double Pe, std::mt19937& rng) {
    int N = 500;
    double a = 1.0;
    double phi = 0.45; // Density high enough for MIPS
    double Dr = 0.5;
    double mu = 1.0;
    double k = 100.0;
    double dt = 0.001;
    int STEPS = 200000;
    
    // Correct L for Area Fraction: phi = (N * PI * (a/2)^2) / L^2
    double L = std::sqrt((N * M_PI * a * a / 4.0) / phi);
    double v0 = Pe * a * Dr;

    std::vector<ABP> P(N);
    std::normal_distribution<double> gauss(0.0, 1.0);
    std::uniform_real_distribution<double> pos_dist(0.0, L);
    std::uniform_real_distribution<double> ang_dist(0.0, 2*M_PI);

    for (int i = 0; i < N; i++) {
        P[i].x = P[i].ux = pos_dist(rng);
        P[i].y = P[i].uy = pos_dist(rng);
        P[i].th = ang_dist(rng);
    }

    double fmax_sum = 0;
    int samples = 0;

    for (int step = 0; step < STEPS; step++) {
        std::vector<double> fx(N, 0.0), fy(N, 0.0);

        // Forces
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = dx_pbc(P[i].x - P[j].x, L);
                double dy = dx_pbc(P[i].y - P[j].y, L);
                double r2 = dx*dx + dy*dy;
                if (r2 < a*a) {
                    double r = std::sqrt(r2);
                    double f = k * (a - r) / r;
                    fx[i] += f * dx; fy[i] += f * dy;
                    fx[j] -= f * dx; fy[j] -= f * dy;
                }
            }
        }

        // Update
        for (int i = 0; i < N; i++) {
            P[i].th += std::sqrt(2.0 * Dr * dt) * gauss(rng);
            double displacement_x = (v0 * std::cos(P[i].th) + mu * fx[i]) * dt;
            double displacement_y = (v0 * std::sin(P[i].th) + mu * fy[i]) * dt;
            
            P[i].ux += displacement_x; // Unwrapped for MSD
            P[i].uy += displacement_y;
            P[i].x = pbc(P[i].x + displacement_x, L);
            P[i].y = pbc(P[i].y + displacement_y, L);
        }

        // Sample fmax in steady state
        if (step > 100000 && step % 2000 == 0) {
            fmax_sum += compute_fmax(P, L, a);
            samples++;
        }
    }
    std::cout << Pe << " " << fmax_sum / samples << std::endl;
}

int main() {
    std::mt19937 rng(42);
    std::cout << "Pe f_max" << std::endl;
    std::vector<double> Pe_list = {10, 30, 50, 80, 100, 150};
    for (double pe : Pe_list) {
        run_simulation(pe, rng);
    }
    return 0;
}
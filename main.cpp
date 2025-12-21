#include <iostream>
#include <vector>
#include <random>

#include "abp_sim.h"
#include "cluster_analysis.h"

double run_abp(
    int N,
    double a,
    double phi,
    double Dr,
    double mu,
    double Pe,
    int STEPS,
    int burnin,
    int sample_every,
    std::mt19937& rng,
    double (*compute_fmax)(const std::vector<ABP>&, double, double)
) {

    double L = std::sqrt(N * M_PI * (a * a * 4.0) / phi);
    double v0 = Pe * a * Dr;

    std::normal_distribution<double> gauss(0.0, 1.0);
    std::uniform_real_distribution<double> pos(0.0, L);
    std::uniform_real_distribution<double> ang(0.0, 2*M_PI);

    std::vector<ABP> P(N);
    std::vector<double> fx(N), fy(N);

    /* init */
    for (int i = 0; i < N; i++) {
        P[i].x = pos(rng);
        P[i].y = pos(rng);
        P[i].th = ang(rng);
    }
    double fsum = 0.0;
    int samples = 0;
    double k = 100.0;
    double dt = 1e-3;

    for (int step = 0; step < STEPS; step++) {

        /* orientation */
        for (int i = 0; i < N; i++) {
            P[i].th += std::sqrt(2.0 * Dr * dt) * gauss(rng);
        }

        /* forces */
        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);

        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {

                double dx = dx_pbc(P[i].x - P[j].x, L);
                double dy = dx_pbc(P[i].y - P[j].y, L);
                double r2 = dx*dx + dy*dy;

                if (r2 < a*a) {
                    double r = std::sqrt(r2);
                    double f = k * (a - r) / r;
                    fx[i] += f * dx;
                    fy[i] += f * dy;
                    fx[j] -= f * dx;
                    fy[j] -= f * dy;
                }
            }
        }

        /* positions */
        for (int i = 0; i < N; i++) {
            P[i].x += (v0 * std::cos(P[i].th) + mu * fx[i]) * dt;
            P[i].y += (v0 * std::sin(P[i].th) + mu * fy[i]) * dt;
            P[i].x = pbc(P[i].x, L);
            P[i].y = pbc(P[i].y, L);
        }

        /* cluster sampling */
        if (step > burnin && step % sample_every == 0) {
            fsum += compute_fmax(P, L, a);
            samples++;
        }
    }

    return fsum / samples;
}

int main() {

    int N = 500;
    double a = 1.0;
    double phi = 0.4;
    double Dr = 0.5;
    double mu = 1.0;

    int STEPS = 300000;
    int burnin = 100000;
    int sample_every = 2000;

    std::mt19937 rng(123);

    std::vector<double> Pe_list = {5, 10, 20, 30, 40, 60, 80};

    for (double Pe : Pe_list) {
        double f_avg = run_abp(
            N, a, phi, Dr, mu, Pe,
            STEPS, burnin, sample_every,
            rng, compute_fmax
        );
        std::cout << Pe << " " << f_avg << std::endl;
    }

    return 0;
}
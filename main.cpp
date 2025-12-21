#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include "abp_sim.h"
#include "cluster_analysis.h"

void run_sim(double Pe, bool save_msd, bool save_snap, std::string fmax_file, std::mt19937& rng) {

    int N = 500; // number of particles
    double a = 1.0; // diameter of one particle
    double Dr = 0.5; // diffusion rate
    double mu = 1.0; // mobility
    double k = 200.0;
    double h = 1e-3; // small chanhge

    double phi = 0.55;

    int STEPS = 300000; // until steady state
    
    double L = std::sqrt((N * M_PI * a * a / 4.0) / phi);
    double v0 = Pe * a * Dr; // since persistence length = v0/aDr

    std::vector<ABP> P(N); // vector of particles N
    std::normal_distribution<double> gauss(0.0, 1.0); // normal gaussian distribution 
    std::uniform_real_distribution<double> pos_dist(0.0, L);
    std::uniform_real_distribution<double> ang_dist(0.0, 2*M_PI);

    for (int i = 0; i < N; i++) {
        P[i].x = P[i].ux = pos_dist(rng); // random positions and orientation for N particles
        P[i].y = P[i].uy = pos_dist(rng);
        P[i].th = ang_dist(rng);
    }

    std::vector<double> initial_ux(N), initial_uy(N); // initial conditions for msd check
    for(int i=0; i<N; i++) { initial_ux[i] = P[i].ux; initial_uy[i] = P[i].uy; }

    std::ofstream msd_out, snap_out; // save data files for simulation later
    if (save_msd) msd_out.open("msd_data.dat"); // msd check data
    if (save_snap) snap_out.open("snapshot.dat"); // of all the particles

    double fmax_sum = 0; int samples = 0;

    for (int step = 0; step < STEPS; step++) { // cluster check
        std::vector<double> fx(N, 0.0), fy(N, 0.0); // vector of zero f
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

        for (int i = 0; i < N; i++) { 
            P[i].th += std::sqrt(2.0 * Dr * h) * gauss(rng);
            double dx_move = (v0 * std::cos(P[i].th) + mu * fx[i]) * h;
            double dy_move = (v0 * std::sin(P[i].th) + mu * fy[i]) * h;
            P[i].ux += dx_move; P[i].uy += dy_move;
            P[i].x = pbc(P[i].x + dx_move, L);
            P[i].y = pbc(P[i].y + dy_move, L);
        }

        if (save_msd && step % 1000 == 0) {
            double msd = 0;
            for(int i=0; i<N; i++) msd += std::pow(P[i].ux-initial_ux[i], 2) + std::pow(P[i].uy-initial_uy[i], 2);
            msd_out << step * h << " " << msd / N << "\n";
        }
        if (step > 150000 && step % 2000 == 0) { fmax_sum += compute_fmax(P, L, a); samples++; }
    }

    std::ofstream f_out(fmax_file, std::ios::app);
    f_out << Pe << " " << fmax_sum / samples << "\n";
    if (save_snap) for(auto& p : P) snap_out << p.x << " " << p.y << "\n";
}

int main() {
    std::mt19937 rng(42);
    std::ofstream f_out("fmax_vs_pe.dat"); f_out << "Pe fmax\n"; f_out.close();
    std::vector<double> Pe_list = {10, 30, 50, 80, 100, 150};
    for (double pe : Pe_list) {
        std::cout << "Simulating Pe=" << pe << "..." << std::endl;
        run_sim(pe, (pe == 50), (pe == 150), "fmax_vs_pe.dat", rng);
    }
    return 0;
}
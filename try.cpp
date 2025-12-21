// implementing minimal ABPs

#include<cmath>
#include<random>
#include <iostream>
#include<vector>
using namespace std;

struct ABP {
    double x, y; // position
    double th; // orientation
}; // this structure defines a particle

int N = 500; // number of particles
double a = 1.0; // diameter of one particle
double L; // box size
double v0 = 20.0; // velocity
double Dr = 0.5; // diffusion rate
double mu = 1.0; // mobility
double h = 1e-3; // small chanhge

double phi = 0.4;

int STEPS = 300000;

mt19937 rng(123);
normal_distribution<double> gauss(0.0, 1.0);

vector<ABP> P(N); // vector or list of N particles

double pbc(double x) {
    if (x < 0)  return x + L;
    if (x >= L) return x - L;
    return x;
}

double dx_pbc(double dx) {
    if (dx >  L/2) dx -= L;
    if (dx < -L/2) dx += L;
    return dx;
}

int main() {
    L = sqrt(N * M_PI * a * a / phi); 

    uniform_real_distribution<double> uniform(0.0, L);
    uniform_real_distribution<double> angle(0.0, 2*M_PI);

    for (int i = 0; i < N; i++) {
        P[i].x  = uniform(rng);
        P[i].y  = uniform(rng);
        P[i].th = angle(rng);
    }
    vector<double> fx(N), fy(N); // force arrays
    vector<double> x0(N), y0(N); // initial positions for msd check later
    for (int i = 0; i < N; i++) {
        x0[i] = P[i].x;
        y0[i] = P[i].y;
    }

    for (int step = 0; step < STEPS; step++) {
        for (int i = 0; i < N; i++) {
            P[i].th += sqrt(2.0 * Dr * h) * gauss(rng); // orientation update
        }

        for (int i = 0; i < N; i++) { //
            fx[i] = 0.0;
            fy[i] = 0.0;
        }

        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {

                double dx = dx_pbc(P[i].x - P[j].x);
                double dy = dx_pbc(P[i].y - P[j].y);
                double r2 = dx*dx + dy*dy;

                if (r2 < a*a) {
                    double r = sqrt(r2);
                    double overlap = a - r;
                    double k = 100.0;          // spring constant
                    double f = k * overlap / r;

                    fx[i] += f * dx;
                    fy[i] += f * dy;
                    fx[j] -= f * dx;
                    fy[j] -= f * dy;
                }
            }
        }

        for (int i = 0; i < N; i++) {
            P[i].x += (v0 * cos(P[i].th) + mu * fx[i]) * h;
            P[i].y += (v0 * sin(P[i].th) + mu * fy[i]) * h;

            P[i].x = pbc(P[i].x);
            P[i].y = pbc(P[i].y);
        }

        /* --- MSD sanity check --- */
        if (step % 5000 == 0 && step < 50000) {
            double msd = 0.0;
            for (int i = 0; i < N; i++) {
                double dx = dx_pbc(P[i].x - x0[i]);
                double dy = dx_pbc(P[i].y - y0[i]);
                msd += dx*dx + dy*dy;
            }
            msd /= N;
            cout << step * h << " " << msd << endl;
        }
    }
    return 0;


}

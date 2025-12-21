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


}

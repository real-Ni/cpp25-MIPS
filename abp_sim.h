#ifndef ABP_SIM_H
#define ABP_SIM_H // header file

#include<cmath>
#include<random>
#include <iostream>
#include<vector>
using namespace std;

struct ABP {
    double x, y; // position
    double th; // orientation
}; // this structure defines a particle

double pbc(double x, double L) {
    if (x < 0)  return x + L;
    if (x >= L) return x - L;
    return x;
}

double dx_pbc(double dx, double L) {
    if (dx >  L/2) dx -= L;
    if (dx < -L/2) dx += L;
    return dx;
} // boundaries period

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
); // function declared, used in main.cpp

#endif

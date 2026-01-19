#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
using namespace std;

const double PI = acos(-1.0);

struct ABP {
    double x, y; // position
    double theta; // orientation
    double ux, uy;
    double x0, y0;   // initial unwrapped positions
};

int N = 2000;
double a = 1.0;
double phi = 0.25;
double L = sqrt(N * PI * a * a / phi);
double sigma = 1.0;
double k = 50.0;

double r_cut = sigma;
double r_list = 1.5 * sigma;

vector<vector<int>> neighbor_list;

// minimum image convention
double minimum_image(double dx) {
    if (dx >  0.5 * L) dx -= L;
    if (dx < -0.5 * L) dx += L;
    return dx;
}

// build neighbour list
void build_neighbor_list(const vector<ABP>& p) {

    neighbor_list.clear();
    neighbor_list.resize(N);

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {

            double dx = minimum_image(p[j].x - p[i].x);
            double dy = minimum_image(p[j].y - p[i].y);
            double r2 = dx*dx + dy*dy;

            if (r2 < r_list * r_list) {
                neighbor_list[i].push_back(j);
                neighbor_list[j].push_back(i);
            }
        }
    }
}

int main() {

    double v0 = 45.0;
    double Dr = 1.0;
    double mu = 1.0;
    double dt = 0.001;
    int steps = 10000;
    int n_rebuild = 20;

    vector<ABP> p(N);
    vector<double> Fx(N), Fy(N);

    mt19937 rng(42);
    uniform_real_distribution<double> pos_dist(0.0, L);
    uniform_real_distribution<double> ang_dist(0.0, 2*M_PI);
    normal_distribution<double> gauss(0.0, 1.0);

    ofstream msd_file("msd.dat"); // will be plotted using python

    for (int i = 0; i < N; i++) {
        p[i].x = pos_dist(rng);
        p[i].y = pos_dist(rng);
        p[i].theta = ang_dist(rng);
    }

    build_neighbor_list(p);

    for (int step = 0; step < steps; step++) {


        if (step % n_rebuild == 0)
            build_neighbor_list(p);

        fill(Fx.begin(), Fx.end(), 0.0);
        fill(Fy.begin(), Fy.end(), 0.0);

        // force calculation using neighbour list
        for (int i = 0; i < N; i++) {
            for (int j : neighbor_list[i]) {

                if (j <= i) continue;

                double dx = minimum_image(p[j].x - p[i].x);
                double dy = minimum_image(p[j].y - p[i].y);
                double r2 = dx*dx + dy*dy;

                if (r2 < r_cut * r_cut && r2 > 1e-12) {
                    double r = sqrt(r2);
                    double f = k * (1.0 - r / sigma);

                    double fx = f * dx / r;
                    double fy = f * dy / r;

                    Fx[i] -= fx;
                    Fy[i] -= fy;
                    Fx[j] += fx;
                    Fy[j] += fy;
                }
            }
        }


        // update particles
        for (int i = 0; i < N; i++) {

            p[i].theta += sqrt(2 * Dr * dt) * gauss(rng);

            double vx = v0 * cos(p[i].theta) + mu * Fx[i];
            double vy = v0 * sin(p[i].theta) + mu * Fy[i];

            p[i].x += vx * dt;
            p[i].y += vy * dt;

            p[i].ux += vx * dt;
            p[i].uy += vy * dt;

            if (p[i].x >= L) p[i].x -= L;
            if (p[i].x <  0) p[i].x += L;
            if (p[i].y >= L) p[i].y -= L;
            if (p[i].y <  0) p[i].y += L;
        }

        if (step % 1000 == 0) {
            double msd = 0.0; // remove prev msd value

            for (int i = 0; i < N; i++) {
                double delx = p[i].ux - p[i].x0;
                double dely = p[i].uy - p[i].y0;
                msd += delx * delx + dely * dely;
            }

            msd /= N;
            double time = step * dt;
            msd_file << time << " " << msd << "\n";
        }

        if (step % 1000 == 0)
            cout << "step " << step << endl;
    }
}
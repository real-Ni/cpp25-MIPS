#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
using namespace std;

const double PI = acos(-1.0);

struct ABP {
    double x, y;    // position
    double theta;   // orientation
};

int N = 500;
double a = 1.0;
double phi = 0.25;
double L = sqrt(N * PI * a * a / phi);

double sigma = 1.0;
double k = 50.0;

double r_cut = 1.2 * sigma;
double r_list = 1.5 * sigma;

vector<vector<int>> neighbor_list;

// minimum image convention
inline double min_image(double dx) {
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
            double dx = min_image(p[j].x - p[i].x);
            double dy = min_image(p[j].y - p[i].y);
            if (dx*dx + dy*dy < r_list * r_list) {
                neighbor_list[i].push_back(j);
                neighbor_list[j].push_back(i);
            }
        }
    }
}

// compute largest cluster fraction f_max
double compute_fmax(const vector<ABP>& p) {

    double r2_cut = (1.2 * a) * (1.2 * a);
    vector<bool> visited(N, false);
    int largest_cluster = 0;

    for (int i = 0; i < N; i++) {
        if (visited[i]) continue;

        int cluster_size = 0;
        vector<int> stack;
        stack.push_back(i);
        visited[i] = true;

        while (!stack.empty()) {
            int u = stack.back();
            stack.pop_back();
            cluster_size++;

            for (int v = 0; v < N; v++) {
                if (visited[v]) continue;

                double dx = min_image(p[v].x - p[u].x);
                double dy = min_image(p[v].y - p[u].y);

                if (dx*dx + dy*dy < r2_cut) {
                    visited[v] = true;
                    stack.push_back(v);
                }
            }
        }

        largest_cluster = max(largest_cluster, cluster_size);
    }

    return double(largest_cluster) / double(N);
}

int main() {

    double Dr = 1.0;
    double mu = 1.0;
    double dt = 0.005;

    int steps = 200000;
    int equil_steps = 50000;
    int sample_every = 1000;
    int rebuild_every = 20;

    vector<ABP> p(N);
    vector<double> Fx(N), Fy(N);

    mt19937 rng(42);
    uniform_real_distribution<double> pos_dist(0.0, L);
    uniform_real_distribution<double> ang_dist(0.0, 2 * PI);
    normal_distribution<double> gauss(0.0, 1.0);

    ofstream out("fmax_vs_pe.dat");

    // sweep Péclet number via v0
    for (double v0 : {10.0, 30.0, 60.0, 90.0, 120.0}) {

        for (int i = 0; i < N; i++) {
            p[i].x = pos_dist(rng);
            p[i].y = pos_dist(rng);
            p[i].theta = ang_dist(rng);
        }

        build_neighbor_list(p);

        vector<double> fmax_samples;

        for (int step = 0; step < steps; step++) {

            if (step % rebuild_every == 0)
                build_neighbor_list(p);

            fill(Fx.begin(), Fx.end(), 0.0);
            fill(Fy.begin(), Fy.end(), 0.0);

            // forces
            for (int i = 0; i < N; i++) {
                for (int j : neighbor_list[i]) {
                    if (j <= i) continue;

                    double dx = min_image(p[j].x - p[i].x);
                    double dy = min_image(p[j].y - p[i].y);
                    double r2 = dx*dx + dy*dy;

                    if (r2 < r_cut * r_cut && r2 > 1e-12) {
                        double r = sqrt(r2);
                        double f = k * (1.0 - r / r_cut);

                        double fx = f * dx / r;
                        double fy = f * dy / r;

                        Fx[i] -= fx;
                        Fy[i] -= fy;
                        Fx[j] += fx;
                        Fy[j] += fy;
                    }
                }
            }

            // update
            for (int i = 0; i < N; i++) {

                p[i].theta += sqrt(2 * Dr * dt) * gauss(rng);

                double vx = v0 * cos(p[i].theta) + mu * Fx[i];
                double vy = v0 * sin(p[i].theta) + mu * Fy[i];

                p[i].x += vx * dt;
                p[i].y += vy * dt;

                if (p[i].x >= L) p[i].x -= L;
                if (p[i].x <  0) p[i].x += L;
                if (p[i].y >= L) p[i].y -= L;
                if (p[i].y <  0) p[i].y += L;
            }

            if (step > equil_steps && step % sample_every == 0) {
                double fmax = compute_fmax(p);
                fmax_samples.push_back(fmax);
            }
        }

        double fmax_avg = 0.0;
        for (double f : fmax_samples) fmax_avg += f;
        fmax_avg /= fmax_samples.size();

        double Pe = v0 / (a * Dr);
        out << Pe << " " << fmax_avg << "\n";

        cout << "Pe = " << Pe << "  <fmax> = " << fmax_avg << endl;
    }

    out.close();
    return 0;
}
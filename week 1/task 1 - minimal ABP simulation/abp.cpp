#include <iostream>
#include<vector>

#include <cmath>
#include <random>
#include <fstream>
using namespace std;

struct ABP {
    double x, y; // position
    double theta; // orientation
    double xu, yu; // unwrapped, for msd calc
    double x0, y0;   // initial unwrapped positions
}; // this structure defines a particle

int main() {
    int N = 500; // number of particles
    double L = 50.0; 
    double v0 = 1.0; // uniform velocity
    double Dr = 0.5; // diffusion rate
    double dt = 0.001;  // small chanhge
    int steps = 10000; // total
    int sample_every = 10; // calculate msd after how many steps

    std::vector<ABP> p(N); // vector of particles N

    mt19937 rng(42); // for rotational noise, note this is dimensionless, while the mathematical eqns have a dimension
    std::normal_distribution<double> gauss(0.0, 1.0); // normal gaussian distribution 
    std::uniform_real_distribution<double> pos_dist(0.0, L);
    std::uniform_real_distribution<double> ang_dist(0.0, 2*M_PI);

    for (int i = 0; i < N; i++) { // initialising the particles
        p[i].x = pos_dist(rng);
        p[i].y = pos_dist(rng);
        p[i].xu = p[i].x;
        p[i].yu = p[i].y;
        p[i].x0 = p[i].xu;
        p[i].y0 = p[i].yu;
        p[i].theta = ang_dist(rng);
    }

    ofstream msd_file("msd.dat"); // will be plotted using python

    for (int step = 1; step <= steps; step++) {

        // updates
        for (int i = 0; i < N; i++) {
            
            //orientation

            p[i].theta += sqrt(2.0 * Dr * dt) * gauss(rng);
            double dx = v0 * cos(p[i].theta) * dt;
            double dy = v0 * sin(p[i].theta) * dt;

            //position

            //unwrapped
            p[i].xu += dx;
            p[i].yu += dy;
            //wrapped
            p[i].x += dx;
            p[i].y += dy;

            // boundary check
            if (p[i].x >= L) p[i].x -= L;
            if (p[i].x <  0) p[i].x += L;
            if (p[i].y >= L) p[i].y -= L;
            if (p[i].y <  0) p[i].y += L;
        }

        // msd calcuation

        if (step % sample_every == 0) {
            double msd = 0.0; // remove prev msd value

            for (int i = 0; i < N; i++) {
                double delx = p[i].xu - p[i].x0;
                double dely = p[i].yu - p[i].y0;
                msd += delx * delx + dely * dely;
            }

            msd /= N;
            double time = step * dt;
            msd_file << time << " " << msd << "\n";
        }
    }

    msd_file.close();
    cout << "MSD written to msd.dat\n";
}
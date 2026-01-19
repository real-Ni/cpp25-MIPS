#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>

using namespace std;

struct Particle {
    double x;
    double y;
    double theta;
};

double run(int N, double eta, double L, int steps) {

    double r = 1.0;
    double v = 0.03;

    int M = (int)(L / r);
    if (M < 1) M = 1;

    vector<Particle> p(N);
    vector<vector<vector<int>>> cell(M, vector<vector<int>>(M));

    for (int i = 0; i < N; i++) {
        p[i].x = ((double)rand() / RAND_MAX) * L;
        p[i].y = ((double)rand() / RAND_MAX) * L;
        p[i].theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
    }

    for (int t = 0; t < steps; t++) {

        if (t % 1000 == 0)
            cout << "    step " << t << "/" << steps << endl;

        for (int i = 0; i < M; i++)
            for (int j = 0; j < M; j++)
                cell[i][j].clear();

        for (int i = 0; i < N; i++) {
            int cx = (int)(p[i].x / r);
            int cy = (int)(p[i].y / r);
            if (cx >= M) cx = M - 1;
            if (cy >= M) cy = M - 1;
            cell[cx][cy].push_back(i);
        }

        vector<double> new_theta(N);

        for (int i = 0; i < N; i++) {
            int cx = (int)(p[i].x / r);
            int cy = (int)(p[i].y / r);

            double sx = 0.0, sy = 0.0;

            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    int nx = (cx + dx + M) % M;
                    int ny = (cy + dy + M) % M;

                    for (int k = 0; k < cell[nx][ny].size(); k++) {
                        int j = cell[nx][ny][k];

                        double dxp = p[j].x - p[i].x;
                        double dyp = p[j].y - p[i].y;

                        if (dxp >  L/2) dxp -= L;
                        if (dxp < -L/2) dxp += L;
                        if (dyp >  L/2) dyp -= L;
                        if (dyp < -L/2) dyp += L;

                        if (dxp*dxp + dyp*dyp <= r*r) {
                            sx += cos(p[j].theta);
                            sy += sin(p[j].theta);
                        }
                    }
                }
            }

            double noise = eta * ((double)rand()/RAND_MAX - 0.5);
            new_theta[i] = atan2(sy, sx) + noise;
        }

        for (int i = 0; i < N; i++) {
            p[i].theta = new_theta[i];
            p[i].x += v * cos(p[i].theta);
            p[i].y += v * sin(p[i].theta);

            if (p[i].x < 0) p[i].x += L;
            if (p[i].x >= L) p[i].x -= L;
            if (p[i].y < 0) p[i].y += L;
            if (p[i].y >= L) p[i].y -= L;
        }
    }

    double vx = 0.0, vy = 0.0;
    for (int i = 0; i < N; i++) {
        vx += cos(p[i].theta);
        vy += sin(p[i].theta);
    }

    return sqrt(vx*vx + vy*vy) / N;
}

int main() {

    srand(time(NULL));
    int steps = 4000;

    cout << "Starting Fig 2(a) data generation" << endl;

    ofstream fa("fig2a.dat");
    double rho = 0.4;
    int Ns[4] = {40, 100, 400, 1000};

    for (int ni = 0; ni < 4; ni++) {
        int N = Ns[ni];
        double L = sqrt(N / rho);

        cout << "Running N = " << N << ", L = " << L << endl;

        for (double eta = 0.2; eta <= 5.0; eta += 0.3) {
            cout << "  eta = " << eta << endl;
            double va = run(N, eta, L, steps);
            fa << N << " " << eta << " " << va << "\n";
        }
        fa << "\n";
    }
    fa.close();

    cout << "Finished Fig 2(a)" << endl;
    cout << "Starting Fig 2(b) data generation" << endl;

    ofstream fb("fig2b.dat");
    double eta = 2.0;
    double L = 20.0;

    for (double rho2 = 0.1; rho2 <= 8.0; rho2 += 0.3) {
        int N = (int)(rho2 * L * L);
        cout << "rho = " << rho2 << ", N = " << N << endl;
        double va = run(N, eta, L, steps);
        fb << rho2 << " " << va << "\n";
    }
    fb.close();

    cout << "Finished Fig 2(b)" << endl;
    cout << "All simulations complete" << endl;

    return 0;
}
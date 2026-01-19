#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <format>

using namespace std;

struct Particle {
    double x;
    double y;
    double theta;
};

int main() {

    int N;
    double eta;
    int steps;
    int snapshot_step;
    double L;

    cin >> N >> eta >> L >> steps >> snapshot_step; // input the nuber of particles, the noise, the length of side, no. of steps

    double r = 1.0; // distance uptill which particles orientations' averages are used
    double v = 0.03; // distance updated
    double rho = N / (L*L);

    int M = (int)(L / r);
    if (M < 1) M = 1;

    vector<Particle> p(N);

    for (int i = 0; i < N; i++) {
        p[i].x = ((double)rand() / RAND_MAX) * L;
        p[i].y = ((double)rand() / RAND_MAX) * L;
        p[i].theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
    } // ics

    vector<vector<vector<int>>> cell(M, vector<vector<int>>(M));

    for (int t = 0; t < steps; t++) {

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                cell[i][j].clear();
            }
        }

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

            double sum_x = 0.0;
            double sum_y = 0.0;

            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {

                    int nx = cx + dx;
                    int ny = cy + dy;

                    if (nx < 0) nx += M;
                    if (ny < 0) ny += M;
                    if (nx >= M) nx -= M;
                    if (ny >= M) ny -= M;

                    for (int k = 0; k < cell[nx][ny].size(); k++) {
                        int j = cell[nx][ny][k];

                        double dxp = p[j].x - p[i].x;
                        double dyp = p[j].y - p[i].y;

                        if (dxp >  L/2) dxp -= L;
                        if (dxp < -L/2) dxp += L;
                        if (dyp >  L/2) dyp -= L;
                        if (dyp < -L/2) dyp += L;

                        if (dxp*dxp + dyp*dyp <= r*r) {
                            sum_x += cos(p[j].theta);
                            sum_y += sin(p[j].theta);
                        }
                    }
                }
            }

            double noise = eta * ((double)rand()/RAND_MAX - 0.5);
            new_theta[i] = atan2(sum_y, sum_x) + noise;
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

        if (t == snapshot_step && snapshot_step >= 0) {
            ofstream out(std::format("snapshot_{:.1f}.dat", rho)); 
            for (int i = 0; i < N; i++) {
                out << p[i].x << " "
                    << p[i].y << " "
                    << p[i].theta << "\n";
            }
            out.close();
        }
    }

    double vx = 0.0, vy = 0.0;
    for (int i = 0; i < N; i++) {
        vx += cos(p[i].theta);
        vy += sin(p[i].theta);
    }

    double va = sqrt(vx*vx + vy*vy) / N;
    cout << va << endl;

    return 0;
}
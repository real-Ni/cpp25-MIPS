#ifndef CLUSTER_ANALYSIS_H
#define CLUSTER_ANALYSIS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "abp_sim.h"

inline double compute_fmax(const std::vector<ABP>& P, double L, double a) {
    int N = P.size();
    double rc = 1.2 * a;
    double rc2 = rc * rc;
    std::vector<bool> visited(N, false);
    int largest = 0;

    for (int i = 0; i < N; i++) {
        if (visited[i]) continue;
        int cluster_size = 0;
        std::vector<int> stack;
        stack.push_back(i);
        visited[i] = true;

        while (!stack.empty()) {
            int k = stack.back();
            stack.pop_back();
            cluster_size++;

            for (int j = 0; j < N; j++) {
                if (!visited[j]) {
                    double dx = dx_pbc(P[k].x - P[j].x, L);
                    double dy = dx_pbc(P[k].y - P[j].y, L);
                    if (dx*dx + dy*dy < rc2) {
                        visited[j] = true;
                        stack.push_back(j);
                    }
                }
            }
        }
        largest = std::max(largest, cluster_size);
    }
    return (double)largest / N;
}

#endif
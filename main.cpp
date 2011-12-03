#include "vp-tree.h"
#include <cmath>
#include <iostream>
#include <numeric>
#include <cassert>
#include <omp.h>

struct Vector {
    Vector(double x, double y) : x(x), y(y) {

    }

    inline
    double length() const {
        return std::sqrt(length_squared());
    }

    inline
    double length_squared() const {
        return x*x + y*y;
    }

    double x, y;
};

inline
Vector operator -(const Vector &v1, const Vector &v2) {
    return Vector(v1.x - v2.x, v1.y - v2.y);
}

struct Particle {
    Particle(double x, double y, double volume) : pos(x, y), volume(volume) {

    }

    Particle(const Vector &pos, double volume) : pos(pos), volume(volume) {

    }

    Vector pos;
    double rho, volume;
};

inline
double dist(const Particle &p1, const Particle &p2) {
    return (p1.pos - p2.pos).length();
}

double kernel(double r, double h) {
    double q = r/h;
    assert(q <= 2);
    if (q > 2) {
        return 0;
    }
    double alpha = 7.0/(4.0*M_PI*h*h);
    double oneqhalve = 1.0-q/2.0;
    double twoqone = 2.0*q+1.0;
    return alpha * oneqhalve*oneqhalve*oneqhalve*oneqhalve * twoqone;
}

#define DX 0.1

void density(std::vector<Particle> &particles,
             const VpTree<Particle, dist> &tree,
             double mass) {
    int nx = sqrt(particles.size());
    int middle = nx/2*nx+nx/2;

    std::vector<Particle> neighbors;
    std::vector<double>   distances;

    size_t k = 40;

    Particle &p = particles[middle];
    tree.search(p, k, &neighbors, &distances);

    assert(neighbors.size() == k);

    double h = (distances[distances.size()-1]) / 2.0 + DX/4.0;

    p.rho = 0.0;
    double sum_Wij = 0.0;

    for (size_t i = 0; i < k; i++) {
        double wij = kernel(distances[i], h);
        p.rho += mass * wij;
        sum_Wij += p.volume * wij;
    }
    std::cout << p.rho << std::endl;
    std::cout << sum_Wij << std::endl;
}

int main() {
    assert(Vector(1.0, 0.0).length() == 1.0);
    assert((Vector(1.0, 0.0) - Vector(2.0, 0.0)).length() == 1.0);

    VpTree<Particle, dist> tree;

    std::vector<Particle> particles;

    double dx = DX;

    double volume = dx*dx;

    double mass = 1000.0 * volume;

    double start = omp_get_wtime();
    for (int x = 0; x < 100; x++) {
        for (int y = 0; y < 100; y++) {
            particles.push_back(Particle(x*dx, y*dx, volume));
        }
    }
    double end = omp_get_wtime();
    std::cout << "creating particles took: " << end - start << " seconds" << std::endl;
    start = omp_get_wtime();
    tree.create(particles);
    end = omp_get_wtime();
    std::cout << "creating tree took: " << end - start << " seconds" << std::endl;

    density(particles, tree, mass);

    exit(0);

    int k = 10;

    start = omp_get_wtime();
#pragma omp parallel for
    for (size_t i = 0; i < particles.size(); i++) {
        std::vector<double> distances;
        std::vector<Particle> neighbors;
        tree.search(particles[i], k, &neighbors, &distances);
    }
    end = omp_get_wtime();
    std::cout << "searching neighbors for all particles took: " << end - start << " seconds" << std::endl;

    return 0;
}

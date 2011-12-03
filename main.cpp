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
    Particle(double x, double y) : pos(x, y) {

    }

    Particle(const Vector &pos) : pos(pos) {

    }

    Vector pos;
};

inline
double dist(const Particle &p1, const Particle &p2) {
    return (p1.pos - p2.pos).length();
}

int main() {
    assert(Vector(1.0, 0.0).length() == 1.0);
    assert((Vector(1.0, 0.0) - Vector(2.0, 0.0)).length() == 1.0);

    VpTree<Particle, dist> tree;

    std::vector<Particle> particles;

    double dx = 0.1;
    double dy = 0.1;

    double start = omp_get_wtime();
    for (int x = 0; x < 1000; x++) {
        for (int y = 0; y < 1000; y++) {
            particles.push_back(Particle(x*dx, y*dy));
        }
    }
    double end = omp_get_wtime();
    std::cout << "creating particles took: " << end - start << " seconds" << std::endl;
    start = omp_get_wtime();
    tree.create(particles);
    end = omp_get_wtime();
    std::cout << "creating tree took: " << end - start << " seconds" << std::endl;

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

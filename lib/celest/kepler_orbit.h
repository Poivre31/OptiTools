#pragma once
#include "../solver/coordinates.h"

class orbit
{
public:
    orbit(double M, double m, double G = 6.67430e-11) : M(M), m(m), G(G) {}

    void set_from_current_state(xy position, xy velocity)
    {
        double r = position.norm();
        double v2 = velocity.dot(velocity);
        xy er = position.normalized();
        double v_rad = velocity.dot(er);
        double v_ortho = sqrt(v2 - v_rad * v_rad);
        double θ = atan2(position[1], position[0]);
        a = r / (2 - r * v2 / (G * M));
        p = r * r * v_ortho * v_ortho / (G * M);
        e = sqrt(1 - p / a);
        if (e != 0)
            θ0 = θ - acos((p - r) / (r * e));
    };

    double get_period()
    {
        return 2 * M_PI * sqrt(a * a * a / (G * M));
    }
    double get_parameter() { return p; }
    double get_semi_major_axis() { return a; }
    double get_semi_minor_axis() { return a * sqrt(1 - e * e); }
    double get_eccentricity() { return e; }
    double get_r_min() { return p / (1 + e); }
    double get_r_max() { return p / (1 - e); }

    std::vector<xy> get_trajectory(size_t N)
    {
        std::vector<polar> position(N + 1);
        double dθ = 2 * M_PI / N;
        for (size_t i = 0; i < N + 1; i++)
        {
            position[i][1] = i * dθ;
            position[i][0] = p / (1 + e * cos(i * dθ - θ0));
        }
        return polar_to_cart(position);
    }

private:
    double a = 1;
    double p = 1;
    double e = 0;
    double G = 6.674e-11;
    double M = 1;
    double m = 1;
    double θ0 = 0;
};
#pragma once
#include <cstdio>
#include <cmath>
#include <functional>
#include <span>

class vec2d
{
public:
	double x{};
	double y{};

	vec2d() = default;

	vec2d(double a) : x(a), y(a) {}

	template <typename T>
	vec2d(T x, T y) : x((double)x), y((double)y) {}

	vec2d(const vec2d &ref) : x(ref.x), y(ref.y) {}

	vec2d &operator=(const vec2d &v)
	{
		x = v.x;
		y = v.y;
		return *this;
	}

	void zero()
	{
		x = {};
		y = {};
	}

	double norm() const
	{
		return std::sqrt(x * x + y * y);
	}

	double norm_2() const
	{
		return x * x + y * y;
	}

	void normalize()
	{
		double r = x * x + y * y;
		if (r == 0)
			return;
		r = 1. / sqrt(r);
		x *= r;
		y *= r;
	}

	double sum() const
	{
		return x + y;
	}

	vec2d round()
	{
		x = (double)std::round(x);
		y = (double)std::round(y);
		return *this;
	}

	void print() const
	{
		std::printf("x: %.4f, y: %.4f\n", x, y);
	}

	void print(const char *message) const
	{
		std::printf("%s x: %.4f, y: %.4f\n", message, x, y);
	}

	static vec2d unitX()
	{
		return vec2d(1, 0);
	}

	static vec2d unitY()
	{
		return vec2d(0, 1);
	}

	bool operator==(const vec2d &v) const
	{
		return (x == v.x && y == v.y);
	}

	vec2d apply_element_wise(std::function<double(double)> func)
	{
		return vec2d(func(x), func(y));
	}

	template <typename U>
	friend vec2d operator*(U a, const vec2d &v)
	{
		vec2d out;
		out.x = v.x * a;
		out.y = v.y * a;
		return out;
	}

	template <typename U>
	friend vec2d operator*(const vec2d &v, U a)
	{
		vec2d out;
		out.x = v.x * a;
		out.y = v.y * a;
		return out;
	}

	friend vec2d operator*(const vec2d &u, const vec2d &v)
	{
		vec2d out;
		out.x = u.x * v.x;
		out.y = u.y * v.y;
		return out;
	}

	template <typename U>
	vec2d &operator*=(U a)
	{
		x *= a;
		y *= a;
		return *this;
	}

	vec2d &operator*=(const vec2d &u)
	{
		x *= u.x;
		y *= u.y;
		return *this;
	}

	friend vec2d operator/(const vec2d &u, const vec2d &v)
	{
		vec2d out;
		out.x = u.x / v.x;
		out.y = u.y / v.y;
		return out;
	}

	template <typename U>
	friend vec2d operator/(const vec2d &u, U a)
	{
		vec2d out;
		double aInv = 1. / a;
		out.x = u.x * aInv;
		out.y = u.y * aInv;
		return out;
	}

	vec2d &operator/=(const vec2d &v)
	{
		x /= v.x;
		y /= v.y;
		return *this;
	}

	template <typename U>
	vec2d &operator/=(U a)
	{
		x /= a;
		y /= a;
		return *this;
	}

	friend vec2d operator+(const vec2d &u, const vec2d &v)
	{
		vec2d out;
		out.x = u.x + v.x;
		out.y = u.y + v.y;
		return out;
	}

	vec2d &operator+=(const vec2d &v)
	{
		x += v.x;
		y += v.y;
		return *this;
	}

	friend vec2d operator-(const vec2d &u, const vec2d &v)
	{
		vec2d out;
		out.x = u.x - v.x;
		out.y = u.y - v.y;
		return out;
	}

	vec2d &operator-=(const vec2d &v)
	{
		x -= v.x;
		y -= v.y;
		return *this;
	}

	vec2d operator-() const
	{
		vec2d out;
		out.x = -x;
		out.y = -y;
		return out;
	}
};

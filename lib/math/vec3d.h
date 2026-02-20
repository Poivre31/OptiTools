#pragma once
#include <cstdio>
#include <cmath>
#include <functional>
#include <span>

class vec3d
{
public:
	double x{};
	double y{};
	double z{};

	vec3d() = default;

	vec3d(double a) : x(a), y(a), z(a) {}

	template <typename T>
	vec3d(T x, T y, T z) : x((double)x), y((double)y), z((double)z) {}

	vec3d(const vec3d &ref) : x(ref.x), y(ref.y), z(ref.z) {}

	vec3d &operator=(const vec3d &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	void zero()
	{
		x = {};
		y = {};
		z = {};
	}

	double norm() const
	{
		return std::sqrt(x * x + y * y + z * z);
	}

	double norm_2() const
	{
		return x * x + y * y + z * z;
	}

	void normalize()
	{
		double r = x * x + y * y + z * z;
		if (r == 0)
			return;
		r = 1. / sqrt(r);
		x *= r;
		y *= r;
		z *= r;
	}

	double sum() const
	{
		return x + y + z;
	}

	vec3d round()
	{
		x = (double)std::round(x);
		y = (double)std::round(y);
		z = (double)std::round(z);
		return *this;
	}

	void print() const
	{
		std::printf("x: %.4f, y: %.4f, z: %.4f \n", x, y, z);
	}

	void print(const char *message) const
	{
		std::printf("%s x: %.4f, y: %.4f, z: %.4f \n", message, x, y, z);
	}

	static vec3d unitX()
	{
		return vec3d(1, 0, 0);
	}

	static vec3d unitY()
	{
		return vec3d(0, 1, 0);
	}

	static vec3d unitZ()
	{
		return vec3d(0, 0, 1);
	}

	bool operator==(const vec3d &v) const
	{
		return (x == v.x && y == v.y && z == v.z);
	}

	vec3d apply_element_wise(std::function<double(double)> func)
	{
		return vec3d(func(x), func(y), func(z));
	}

	/*template<typename U>
	operator std::span<U, 3>() {
		U out[3] = { U(x), U(y), U(z) };
		return std::span<U, 3>(out);
	}*/

	template <typename U>
	friend vec3d operator*(U a, const vec3d &v)
	{
		vec3d out;
		out.x = v.x * a;
		out.y = v.y * a;
		out.z = v.z * a;
		return out;
	}

	template <typename U>
	friend vec3d operator*(const vec3d &v, U a)
	{
		vec3d out;
		out.x = v.x * a;
		out.y = v.y * a;
		out.z = v.z * a;
		return out;
	}

	friend vec3d operator*(const vec3d &u, const vec3d &v)
	{
		vec3d out;
		out.x = u.x * v.x;
		out.y = u.y * v.y;
		out.z = u.z * v.z;
		return out;
	}

	template <typename U>
	vec3d &operator*=(U a)
	{
		x *= a;
		y *= a;
		z *= a;
		return *this;
	}

	vec3d &operator*=(const vec3d &u)
	{
		x *= u.x;
		y *= u.y;
		z *= u.z;
		return *this;
	}

	friend vec3d operator/(const vec3d &u, const vec3d &v)
	{
		vec3d out;
		out.x = u.x / v.x;
		out.y = u.y / v.y;
		out.z = u.z / v.z;
		return out;
	}

	template <typename U>
	friend vec3d operator/(const vec3d &u, U a)
	{
		vec3d out;
		double aInv = 1. / a;
		out.x = u.x * aInv;
		out.y = u.y * aInv;
		out.z = u.z * aInv;
		return out;
	}

	vec3d &operator/=(const vec3d &v)
	{
		x /= v.x;
		y /= v.y;
		z /= v.z;
		return *this;
	}

	template <typename U>
	vec3d &operator/=(U a)
	{
		x /= a;
		y /= a;
		z /= a;
		return *this;
	}

	friend vec3d operator+(const vec3d &u, const vec3d &v)
	{
		vec3d out;
		out.x = u.x + v.x;
		out.y = u.y + v.y;
		out.z = u.z + v.z;
		return out;
	}

	vec3d &operator+=(const vec3d &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	friend vec3d operator-(const vec3d &u, const vec3d &v)
	{
		vec3d out;
		out.x = u.x - v.x;
		out.y = u.y - v.y;
		out.z = u.z - v.z;
		return out;
	}

	vec3d &operator-=(const vec3d &v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	vec3d operator-() const
	{
		vec3d out;
		out.x = -x;
		out.y = -y;
		out.z = -z;
		return out;
	}
};

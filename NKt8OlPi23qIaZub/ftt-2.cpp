#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4.1,sse4.2,avx,avx2,popcnt,tune=native")

#include "bits/stdc++.h"

#define mem(x) memset((x), 0, sizeof((x)))
#define il __attribute__((always_inline))

using namespace std;

typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;

#if __cplusplus > 201403L
#define r
#else
#define r register
#endif

#define c const

namespace _c
{
c double pi = acos(-1.0);
namespace min
{
c int i8 = -128;
c int i16 = -32768;
c int i = -2147483647 - 1;
c ll l = -9223372036854775807LL - 1;
} // namespace min

namespace max
{
c int i8 = 127;
c int i16 = 32767;
c int i = 2147483647;
c ll l = 9223372036854775807LL;
} // namespace max
} // namespace _c

namespace _f
{
template <typename T>
inline c T gcd(T m, T n)
{
    while (n != 0)
    {
        T t = m % n;
        m = n;
        n = t;
    }
    return m;
}

template <typename T>
inline c T abs(c T &a)
{
    return a > 0 ? a : -a;
}

template <typename T>
inline T pow(T a, T b)
{
    T res = 1;
    while (b > 0)
    {
        if (b & 1)
        {
            res = res * a;
        }
        a = a * a;
        b >>= 1;
    }
    return res;
}

template <typename T>
inline T pow(T a, T b, c T &m)
{
    a %= m;
    T res = 1;
    while (b > 0)
    {
        if (b & 1)
        {
            res = res * a % m;
        }
        a = a * a % m;
        b >>= 1;
    }
    return res % m;
}
} // namespace _f

namespace io
{
template <typename T>
inline T read()
{
    r T res = 0, neg = 1;
    char g = getchar();
    for (; !isdigit(g); g = getchar())
    {
        if (g == '-')
        {
            neg = -1;
        }
    }
    for (; isdigit(g); g = getchar())
    {
        res = res * 10 + g - '0';
    }
    return res * neg;
}
template <typename T>
inline void read(T &t)
{
    t = read<T>();
}
template <typename T>
inline void readln(c T first, c T last)
{
    for (r T it = first; it != last; it++)
    {
        read(*it);
    }
}

template <typename T>
inline void _write(T x)
{
    if (x < 0)
    {
        putchar('-');
        x = -x;
    }
    if (x > 9)
    {
        _write(x / 10);
    }
    putchar(x % 10 + '0');
}
template <typename T>
inline void write(c T &x, c char &sep = ' ')
{
    _write(x);
    putchar(sep);
}
template <typename T>
inline void writeln(c T &x)
{
    write(x, '\n');
}
template <typename T>
inline void writeln(c T first, c T last, c char &sep = ' ', c char &ends = '\n')
{
    for (r T it = first; it != last; it++)
    {
        write(*it, sep);
    }
    putchar(ends);
}

#if __cplusplus >= 201103L
template <typename T, typename... Args>
void read(T &x, Args &... args)
{
    read(x);
    read(args...);
}
#endif
} // namespace io
#undef c
#undef r

int n, m, k;
const int N = 3e6 + 5;
int Lim = 1, L, R[N];

typedef complex<double> comp;
comp A[N], B[N];

inline void FFT(comp *J, int flag)
{
    for (register int i = 0; i < Lim; i++)
    {
        if (i < R[i])
        {
            swap(J[i], J[R[i]]);
        }
    }

    for (register int j = 1; j < Lim; j <<= 1)
    {
        comp T(cos(_c::pi / j), flag * sin(_c::pi / j));

        for (register int k = 0; k < Lim; k += (j << 1))
        {
            comp t(1, 0);

            for (register int l = 0; l < j; l++, t *= T)
            {
                comp Nx = J[k + l], Ny = t * J[k + j + l];
                J[k + l] = Nx + Ny;
                J[k + j + l] = Nx - Ny;
            }
        }
    }
}

int main()
{
    io::read(n, m);

    for (register int i = 0; i <= n; i++)
    {
        A[i].real(io::read<int>());
    }
    for (register int i = 0; i <= m; i++)
    {
        B[i].real(io::read<int>());
    }

    while (Lim <= n + m)
    {
        Lim <<= 1;
        L++;
    }

    for (register int i = 0; i < Lim; i++)
    {
        R[i] = (R[i >> 1] >> 1 | (i & 1) << (L - 1));
    }

    FFT(A, 1), FFT(B, 1);

    for (register int i = 0; i <= Lim; i++)
    {
        A[i] *= B[i];
    }

    FFT(A, -1);

    for (register int i = 0; i <= n + m; i++)
    {
        io::write((int)(A[i].real() / Lim + 0.5));
    }
}
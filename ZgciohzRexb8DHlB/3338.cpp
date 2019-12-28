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

const ld eps = 1e-6;
const int N = 2e5 + 5;

typedef complex<ld> comp;

comp A[N * 4], B[N * 4], C[N * 4];
ld q[N], res[N];

int n, len;

#define ri register int
inline void init(const int typ)
{
    for (ri i = 0; i <= len; i++)
    {
        A[i] = B[i] = 0;
    }

    for (ri i = 1; i <= n; i++)
    {
        B[i] = 1.0 / i / i;
    }

    if (typ)
    {
        for (ri i = 1; i <= n; i++)
        {
            A[i] = q[i];
        }
    }
    else
    {
        for (ri i = 1; i <= n; i++)
        {
            A[i] = q[n - i + 1];
        }
    }
}

inline void FFT(comp *y, const int len, const int typ)
{
    int lg2 = log2(len);

    for (ri i = 0, rev; i < len; i++)
    {
        rev = 0;
        for (ri j = 0; j < lg2; j++)
        {
            rev |= ((i >> j) & 1) << (lg2 - j - 1);
        }
        if (rev < i)
        {
            swap(y[rev], y[i]);
        }
    }

    comp wp, w, we, wo;

    for (ri s = 2, p; s <= len; s <<= 1)
    {
        p = s >> 1;
        wp = comp(cos(typ * _c::pi / p), sin(typ * _c::pi / p));
        for (ri i = 0; i < len; i += s)
        {
        }
    }
}
#undef ri
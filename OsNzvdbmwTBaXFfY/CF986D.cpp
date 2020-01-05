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

const int MOD = 998244353, G = 3, N = 2e6 + 5;
int n;
int a[N], b[N], c[N], rev[N];

inline int ksm(int x, int y)
{
    int res = 1;
    for (; y; y >>= 1, x = (ll)x * x % MOD)
    {
        if (y & 1)
        {
            res = (ll)res * x % MOD;
        }
    }
    return res;
}

#define ri register int
inline void NTT(int *a, int n, bool flag)
{
    for (ri i = 0; i < n; i++)
    {
        if (i < rev[i])
        {
            swap(a[i], a[rev[i]]);
        }
    }

    for (ri i = 1; i < n; i <<= 1)
    {
        ri gn = ksm(G, (MOD - 1) / (i << 1));
        for (ri j = 0; j < n; j += (i << 1))
        {
            ri t1, t2, g = 1;
            for (ri k = 0; k < i; k++, g = (ll)g * gn % MOD)
            {
                t1 = a[j + k], t2 = (ll)g * a[j + k + i] % MOD;
                a[j + k] = (t1 + t2) % MOD;
                a[j + k + i] = (t1 - t2 + MOD) % MOD;
            }
        }
    }

    if (!flag)
    {
        int ny = ksm(n, MOD - 2);

        for (int i = 1; i < ((1 + n) >> 1); i++)
        {
            swap(a[i], a[n - i]);
        }

        for (ri i = 0; i < n; i++)
        {
            a[i] = (ll)a[i] * ny % MOD;
        }
    }
}

void work(int dep, int *a, int *b)
{
    if (dep == 1)
    {
        b[0] = ksm(a[0], MOD - 2);
        return;
    }

    work((dep + 1) >> 1, a, b);

    ri Len = 0, L = 1;

    while (L < (dep << 1))
    {
        L <<= 1;
        Len++;
    }

    for (ri i = 1; i < L; i++)
    {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (Len - 1));
    }

    for (ri i = 0; i < dep; i++)
    {
        c[i] = a[i];
    }

    for (ri i = dep; i < L; i++)
    {
        c[i] = 0;
    }

    NTT(c, L, 1);
    NTT(b, L, 1);

    for (ri i = 0; i < L; i++)
    {
        b[i] = (ll)(2 - (ll)c[i] * b[i] % MOD + MOD) % MOD * b[i] % MOD;
    }

    NTT(b, L, 0);

    for (ri i = dep; i < L; i++)
    {
        b[i] = 0;
    }
}
int main()
{
    n = io::read<int>();
    io::readln(a, a + n);

    work(n, a, b);

    io::writeln(b, b + n);
}

#undef ri
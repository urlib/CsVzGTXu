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

const int N = 300005;
const int MOD = 998244353;

int n, k, factorial[N], ny[N], a[N * 4], b[N * 4], rev[N * 4], L, ans[N * 2], m, q, cnt[N];

#define ri register int

inline void pre(int n)
{
    int lg = 0;
    for (L = 1; L <= n; L <<= 1, lg++)
        ;
    for (ri i = 0; i < L; i++)
    {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (lg - 1));
    }
}

inline int ksm(int x, int y)
{
    int ans = 1;
    while (y)
    {
        if (y & 1)
        {
            ans = (ll)ans * x % MOD;
        }
        x = (ll)x * x % MOD;
        y >>= 1;
    }
    return ans;
}

inline int C(int n, int m)
{
    return (ll)factorial[n] * ny[m] % MOD * ny[n - m] % MOD;
}

inline void NTT(int *a, int f)
{
    for (ri i = 0; i < L; i++)
    {
        if (i < rev[i])
        {
            swap(a[i], a[rev[i]]);
        }
    }

    for (ri i = 1; i < L; i <<= 1)
    {
        int wn = ksm(3, f == 1 ? (MOD - 1) / i / 2 : MOD - 1 - (MOD - 1) / i / 2);
        for (ri j = 0; j < L; j += (i << 1))
        {
            int w = 1;
            for (ri k = 0; k < i; k++)
            {
                int u = a[j + k], v = (ll)a[j + k + i] * w % MOD;
                a[j + k] = (u + v) % MOD;
                a[j + k + i] = (u + MOD - v) % MOD;
                w = (ll)w * wn % MOD;
            }
        }
    }
    if (f == -1)
    {
        for (ri i = 0, inv = ksm(L, MOD - 2); i < L; i++)
        {
            a[i] = (ll)a[i] * inv % MOD;
        }
    }
}

int main()
{
    io::read(n, k);
    for (ri i = 1, x; i <= n; i++)
    {
        io::read(x);
        m = std::max(m, x);
        cnt[x]++;
    }

    factorial[0] = factorial[1] = ny[0] = ny[1] = 1;
    for (ri i = 2; i <= n; i++)
    {
        factorial[i] = (ll)factorial[i - 1] * i % MOD;
        ny[i] = (ll)(MOD - MOD / i) * ny[MOD % i] % MOD;
    }

    for (ri i = 2; i <= n; i++)
    {
        ny[i] = (ll)ny[i - 1] * ny[i] % MOD;
    }

    pre(n * 2);

    while (k--)
    {
        int h = io::read<int>();

        int s1 = 0, s2 = 0;
        for (ri i = 1; i < h; i++)
        {
            if (cnt[i] == 1)
            {
                s2++;
            }
            else if (cnt[i] >= 2)
            {
                s1 += 2;
            }
        }
        for (ri i = 0; i < L; i++)
        {
            a[i] = b[i] = 0;
        }
        for (ri i = 0; i <= s1; i++)
        {
            a[i] = C(s1, i);
        }
        for (ri i = 0; i <= s2; i++)
        {
            b[i] = (ll)ksm(2, i) * C(s2, i) % MOD;
        }

        NTT(a, 1);
        NTT(b, 1);
        for (ri i = 0; i < L; i++)
        {
            a[i] = (ll)a[i] * b[i] % MOD;
        }
        NTT(a, -1);

        for (ri i = 0; i <= s1 + s2; i++)
        {
            (ans[h + i + 1] += a[i]) %= MOD;
        }
    }

    io::read(q);
    for (ri _ = 0, x; _ < q; _++)
    {
        io::read(x);
        io::writeln(ans[x / 2]);
    }
}

#undef ri
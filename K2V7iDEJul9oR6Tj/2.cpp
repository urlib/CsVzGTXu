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
inline c T max(c T &a, c T &b)
{
    return a > b ? a : b;
}

template <typename T>
inline c T min(c T &a, c T &b)
{
    return a < b ? a : b;
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

const int N = 4005;
const int INF = 0x7f7f7f7f;

int n, m, L[N], tot_l = 0, q[N];
bool exist[N] = {0}, used[N] = {0};

int dijkstra()
{
    memset(q, 0x7f, sizeof(q));
    q[0] = 0;
    int v, k;
    for (register int i = 2; i <= 3000; i++)
    {
        if (exist[i])
        {
            L[tot_l++] = i;
        }
    }
    int MOD = L[0];
    while (1)
    {
        k = -1;
        for (register int i = 0; i < MOD; i++)
        {
            if (!used[i] && (k == -1 || q[i] < q[k]))
            {
                k = i;
            }
        }
        if (k == -1)
        {
            break;
        }
        used[k] = 1;
        for (register int i = 1; i < tot_l; i++)
        {
            if (!used[v = (k + L[i]) % MOD])
            {
                q[v] = std::min(q[v], q[k] + L[i]);
            }
        }
    }
    int res = -INF;
    for (register int i = 1; i < MOD; i++)
    {
        res = std::max(res, q[i] - MOD);
    }
    return res;
}

int main()
{
    io::read(n, m);
    int li, _gcd = 0;
    for (int i = 0; i < n; i++)
    {
        io::read(li);
        _gcd = _f::gcd(_gcd, li);
        for (int j = 0; j <= m && j < li; j++)
        {
            exist[li - j] = 1;
            _gcd = _f::gcd(li - j, _gcd);
        }
    }
    if (exist[1] || _gcd > 1)
    {
        io::writeln(-1);
    }
    else
    {
        io::writeln(dijkstra());
    }
    return 0;
}
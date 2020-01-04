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

mt19937 mt_rand(time(0));
const int N = 100000 + 5;
int ch[N][2];
int val[N], pri[N], siz[N], sz = 0;

#define lson(x) ch[(x)][0]
#define rson(x) ch[(x)][1]
inline void update(int x)
{
    siz[x] = 1 + siz[lson(x)] + siz[rson(x)];
}

int new_node(int v)
{
    siz[++sz] = 1;
    val[sz] = v;
    pri[sz] = mt_rand();
    return sz;
}

int merge(int u, int v)
{
    if (!u || !v)
    {
        return u + v;
    }

    if (pri[u] < pri[v])
    {
        rson(u) = merge(rson(u), v);
        update(u);
        return u;
    }
    else
    {
        lson(v) = merge(u, lson(v));
        update(v);
        return v;
    }
}

void split(int now, int k, int &u, int &v)
{
    if (!now)
    {
        u = v = 0;
    }
    else
    {
        if (val[now] <= k)
        {
            u = now;
            split(rson(now), k, rson(now), v);
        }
        else
        {
            v = now;
            split(lson(now), k, u, lson(now));
        }
        update(now);
    }
}

int kth(int now, int k)
{
    while (k != siz[lson(now)] + 1)
    {
        if (k <= siz[lson(now)])
        {
            now = lson(now);
        }
        else
        {
            k -= siz[lson(now)] + 1;
            now = rson(now);
        }
    }
    return now;
}
#undef lson
#undef rson

int main()
{
    int n, Min;
    io::read(n, Min);

    int rt = 0, delta = 0, left = 0;

    int x, y;
    for (int i = 1; i <= n; i++)
    {
        char st[2];
        scanf("%s", st);
        int a = io::read<int>();

        switch (st[0])
        {
        case 'I':
        {
            if (a >= Min)
            {
                a -= delta;
                split(rt, a, x, y);
                rt = merge(x, merge(new_node(a), y));
            }
            break;
        }

        case 'A':
        {
            delta += a;
            break;
        }

        case 'S':
        {
            delta -= a;
            split(rt, Min - delta - 1, x, y);
            rt = y;
            left += siz[x];
            break;
        }

        case 'F':
        {
            if (siz[rt] < a)
            {
                io::writeln(-1);
            }
            else
            {
                io::writeln(val[kth(rt, siz[rt] - a + 1)] + delta);
            }
            break;
        }
        default:
        {
            break;
        }
        }
    }
    io::writeln(left);
}
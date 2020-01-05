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

const int N = 1e5 + 5;

char s[N], t[N];

int n;
int head[N], to[N * 2], nxt[N * 2], E = 1;
int now[N * 2], target[N * 2];

int ans[N], a[N], ban[N];

inline void add_edge(int u, int v, int _now, int _target)
{
    to[E] = v;
    nxt[E] = head[u];

    now[E] = _now;
    target[E] = _target;

    head[u] = E++;
}

void dfs(int u, int fa)
{
    int sum = 0, sum_1 = 0;
#define w ans[u]
    for (int i = head[u]; i; i = nxt[i])
    {
#define v to[i]
        if (v != fa)
        {
            if (target[i])
            {
                if (now[i])
                {
                    ans[v] = 1;
                }
                else
                {
                    ban[v] = 1;
                }
            }

            dfs(v, u);

            if (a[v])
            {
                sum_1++;
                sum += ans[v] - 1;
            }
            else
            {
                sum += ans[v];
            }
        }
#undef v
    }

    if (ban[u])
    {
        if (sum_1 & 1)
        {
            (sum_1 >>= 1)++;
        }
        else
        {
            sum_1 >>= 1;
        }
        ans[u] += sum + sum_1;
    }
    else
    {
        if (w || sum_1 & 1)
        {
            a[u] = 1;
        }
        if (sum_1 & 1 && !w)
        {
            (sum_1 >>= 1)++;
        }
        else
        {
            sum_1 >>= 1;
        }
        ans[u] += sum + sum_1;
    }

#undef w
}

int main()
{
    io::read(n);
    scanf("%s", s + 1);
    scanf("%s", t + 1);

    for (register int i = 1, u, v; i < n; i++)
    {
        io::read(u, v);
        add_edge(u, v, s[i] - '0', t[i] - '0');
        add_edge(v, u, s[i] - '0', t[i] - '0');
    }

    dfs(1, 1);
    io::writeln(ans[1]);
}
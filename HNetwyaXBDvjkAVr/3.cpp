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

mt19937 rnd(time(0));

const int N = 80000 + 5;
int id[N], a[N];

int n, m;
int rt, r1, r2, r3, r4, cnt = 0;

char opt[2];

struct fhq
{
    int ch[2], fa, size, rd, val;
} t[N];

#define lson(x) t[(x)].ch[0]
#define rson(x) t[(x)].ch[1]

inline int newnode(int val)
{
    t[++cnt].val = val;
    t[cnt].rd = rnd();
    t[cnt].size = 1;
    return id[val] = cnt;
}

inline void up(int x)
{
    t[x].size = t[lson(x)].size + t[rson(x)].size + 1;
}

void split(int x, int k, int &a, int &b, int Fa = 0, int Fb = 0)
{
    if (x == 0)
    {
        a = b = 0;
        return;
    }
    if (k <= t[lson(x)].size)
    {
        t[x].fa = Fb;
        b = x;
        split(lson(x), k, a, lson(x), Fa, x);
    }
    else
    {
        t[x].fa = Fa;
        a = x;
        split(rson(x), k - t[lson(x)].size - 1, rson(x), b, x, Fb);
    }
    up(x);
}

int merge(int x, int y)
{
    if (x == 0 || y == 0)
    {
        return x + y;
    }
    if (t[x].rd < t[y].rd)
    {
        rson(x) = merge(rson(x), y);
        t[rson(x)].fa = x;
        up(x);
        return x;
    }
    else
    {
        lson(y) = merge(x, lson(y));
        t[lson(y)].fa = y;
        up(y);
        return y;
    }
}

inline void insert(int pos, int val)
{
    split(rt, pos, r1, r2);
    rt = merge(r1, merge(newnode(val), r2));
}

inline bool get(int x)
{
    return rson(t[x].fa) == x;
}

int find(int cnt)
{
    int node = cnt, res = t[lson(cnt)].size + 1;
    while (node != rt && cnt)
    {
        if (get(cnt))
        {
            res += t[lson(t[cnt].fa)].size + 1;
        }
        cnt = t[cnt].fa;
    }
    return res;
}

int main()
{
    io::read(n, m);
    for (int i = 1; i <= n; i++)
    {
        io::read(a[i]);
        insert(i - 1, a[i]);
    }
    for (int i = 1; i <= m; i++)
    {
        int x, y, k;
        scanf("%s", opt);
        io::read(x);

        switch (opt[0])
        {
        case 'T':
        {
            k = find(id[x]);
            split(rt, k, r1, r3);
            split(r1, k - 1, r1, r2);
            rt = merge(r2, merge(r1, r3));
            break;
        }

        case 'B':
        {
            k = find(id[x]);
            split(rt, k, r1, r3, 0);
            split(r1, k - 1, r1, r2, 0);
            rt = merge(r1, merge(r3, r2));
            break;
        }

        case 'I':
        {
            io::read(y);
            k = find(id[x]);
            if (y)
            {
                if (y > 0)
                {
                    split(rt, k + 1, r3, r4);
                    split(r3, k, r2, r3);
                    split(r2, k - 1, r1, r2);
                    rt = merge(r1, merge(r3, merge(r2, r4)));
                }
                else
                {
                    split(rt, k, r3, r4);
                    split(r3, k - 1, r2, r3);
                    split(r2, k - 2, r1, r2);
                    rt = merge(r1, merge(r3, merge(r2, r4)));
                }
            }
            break;
        }

        case 'A':
        {
            k = find(id[x]);
            io::writeln(k - 1);
            break;
        }

        case 'Q':
        {
            split(rt, x, r1, r2);
            int node = r1;
            while (rson(node))
            {
                node = rson(node);
            }
            io::writeln(t[node].val);
            rt = merge(r1, r2);
            break;
        }

        default:
        {
            break;
        }
        }
    }
    return 0;
}
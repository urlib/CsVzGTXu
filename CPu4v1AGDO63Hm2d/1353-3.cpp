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

typedef __float128 Ld;

mt19937 mt_r(time(0));

const int N = 100005;
const Ld one = 1.0;

struct node
{
    int l, r, pri, rev, sz;
    ll val, sum, sum_2, Min, Max;
} tree[N];

ll a[N];
int n, m;
int tot = 0, rt;
char op[5];

inline int new_node(ll u)
{
    tree[++tot].sz = 1;
    tree[tot].pri = mt_r();
    tree[tot].sum = tree[tot].val = tree[tot].Max = tree[tot].Min = u;
    tree[tot].sum_2 = u * u;
    return tot;
}

inline void push_up(int u)
{
    tree[u].sz = tree[tree[u].l].sz + tree[tree[u].r].sz + 1;
    tree[u].sum = tree[u].val + tree[tree[u].l].sum + tree[tree[u].r].sum;
    tree[u].sum_2 = tree[u].val * tree[u].val + tree[tree[u].l].sum_2 + tree[tree[u].r].sum_2;
    tree[u].Max = tree[u].Min = tree[u].val;

    if (tree[u].l)
    {
        tree[u].Max = max(tree[u].Max, tree[tree[u].l].Max);
        tree[u].Min = min(tree[u].Min, tree[tree[u].l].Min);
    }
    if (tree[u].r)
    {
        tree[u].Max = max(tree[u].Max, tree[tree[u].r].Max);
        tree[u].Min = min(tree[u].Min, tree[tree[u].r].Min);
    }
}
inline void push_down(int u)
{
    if (tree[u].rev)
    {
        tree[u].rev = 0;
        swap(tree[u].l, tree[u].r);
        if (tree[u].l)
        {
            tree[tree[u].l].rev ^= 1;
        }
        if (tree[u].r)
        {
            tree[tree[u].r].rev ^= 1;
        }
    }
}

void split(int u, int x, int &l, int &r)
{
    if (!u)
    {
        l = r = 0;
        return;
    }
    push_down(u);
    if (tree[tree[u].l].sz < x)
    {
        l = u;
        split(tree[u].r, x - tree[tree[u].l].sz - 1, tree[u].r, r);
    }
    else
    {
        r = u;
        split(tree[u].l, x, l, tree[u].l);
    }
    push_up(u);
}

int merge(int l, int r)
{
    if (!l || !r)
    {
        return l + r;
    }
    if (tree[l].pri < tree[r].pri)
    {
        push_down(l);
        tree[l].r = merge(tree[l].r, r);
        push_up(l);
        return l;
    }
    else
    {
        push_down(r);
        tree[r].l = merge(l, tree[r].l);
        push_up(r);
        return r;
    }
}

inline void rev(int l, int r)
{
    int A, B, C;
    split(rt, l - 1, A, B);
    split(B, r - l + 1, B, C);
    tree[B].rev ^= 1;
    rt = merge(merge(A, B), C);
}

inline void solve(int l, int r)
{
    int A, B, C;
    split(rt, l - 1, A, B);
    split(B, r - l + 1, B, C);
    printf("%.6lf\n", double(std::min(
                                 (tree[B].sum_2 - tree[B].Max * tree[B].Max) * (ll)(r - l) - (tree[B].sum - tree[B].Max) * (tree[B].sum - tree[B].Max),
                                 (tree[B].sum_2 - tree[B].Min * tree[B].Min) * (ll)(r - l) - (tree[B].sum - tree[B].Min) * (tree[B].sum - tree[B].Min)) /
                             ld((ll) (r - l + 1) * (r - l))));
    rt = merge(merge(A, B), C);
}

int main()
{
    io::read(n, m);
    for (int i = 1; i <= n; i++)
    {
        rt = merge(rt, new_node(io::read<ll>()));
    }
    for (int i = 1, l, r; i <= m; i++)
    {
        scanf("%s", op);
        io::read(l, r);
        if (op[0] == 'r')
        {
            rev(l, r);
        }
        else
        {
            if (l == r)
            {
                puts("0.000000");
            }
            else
            {
                solve(l, r);
            }
        }
    }
    return 0;
}
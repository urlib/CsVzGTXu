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

const int N = 2e5 + 5;
int n, m, a[N], bl[N], block_size, p, q;

struct node
{
    int l, r, tim, i;
    const bool operator<(const node &rhs) const
    {
        if (bl[l] ^ bl[rhs.l])
        {
            return bl[l] < bl[rhs.l];
        }
        if (bl[r] ^ bl[rhs.r])
        {
            return bl[r] < bl[rhs.r];
        }
        return ((bl[l] & 1) ^ (bl[r] & 1)) ? tim < rhs.tim : tim > rhs.tim;
    }
} ask[N];

int lstv[N], idx[N], val[N];
int tmp[N * 2], ord[N * 2], tot = 0, lst = 0;
void init()
{
    sort(tmp + 1, tmp + tot + 1);
    for (int i = 1; i <= tot; i++)
    {
        if (i == 1 || tmp[i] != tmp[i - 1])
        {
            ord[lst++] = tmp[i];
        }
    }
    for (int i = 1; i <= n; i++)
    {
        a[i] = lower_bound(ord, ord + lst, a[i]) - ord;
    }
    for (int i = 1; i <= q; i++)
    {
        val[i] = lower_bound(ord, ord + lst, val[i]) - ord;
        lstv[i] = lower_bound(ord, ord + lst, lstv[i]) - ord;
    }
}

int cnt[N * 2], sta[N];

inline void inc(int x)
{
    sta[cnt[x]++]--;
    sta[cnt[x]]++;
}

inline void dec(int x)
{
    sta[cnt[x]--]--;
    sta[cnt[x]]++;
}

inline int mex()
{
    for (int i = 1; i <= n; i++)
    {
        if (sta[i] == 0)
        {
            return i;
        }
    }
}

inline void modify(int x, int l, int r)
{
    a[idx[x]] = val[x];
    if (idx[x] < l || r < idx[x])
    {
        return;
    }
    dec(lstv[x]);
    inc(val[x]);
}

inline void unmod(int x, int l, int r)
{
    a[idx[x]] = lstv[x];
    if (idx[x] < l || r < idx[x])
    {
        return;
    }
    dec(val[x]);
    inc(lstv[x]);
}

int ans[N], pre[N];

int main()
{
    io::read(n, m);
    block_size = pow(n, 2.0 / 3);
    for (int i = 1; i <= n; i++)
    {
        bl[i] = i / block_size;
    }

    io::readln(a + 1, a + n + 1);
    for (int i = 1; i <= n; i++)
    {
        tmp[++tot] = a[i];
    }

    p = 0, q = 0;
    for (int i = 1; i <= m; i++)
    {
        int op = io::read<int>();
        if (op == 1)
        {
            ask[++p].l = io::read<int>();
            ask[p].r = io::read<int>();
            ask[p].i = p;
            ask[p].tim = q;
        }
        else
        {
            int x = idx[++q] = io::read<int>();
            tmp[++tot] = val[q] = io::read<int>();
            lstv[q] = pre[x] ? val[pre[x]] : a[x];
            pre[x] = q;
        }
    }
    init();
    sort(ask + 1, ask + p + 1);
    int l = 1, r = 0, cur = 0;
    for (int i = 1; i <= p; i++)
    {
        while (cur < ask[i].tim)
        {
            modify(++cur, l, r);
        }
        while (cur > ask[i].tim)
        {
            unmod(cur--, l, r);
        }
        while (l > ask[i].l)
        {
            inc(a[--l]);
        }
        while (r < ask[i].r)
        {
            inc(a[++r]);
        }
        while (l < ask[i].l)
        {
            dec(a[l++]);
        }
        while (r > ask[i].r)
        {
            dec(a[r--]);
        }
        ans[ask[i].i] = mex();
    }
    io::writeln(ans + 1, ans + p + 1, '\n');
}
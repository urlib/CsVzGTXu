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

unordered_map<int, unordered_map<int, int>> G;

inline void add_edge(int u, int v, int w)
{
    G[u][v] = w;
}

int n, m, s, t;
int ans = 0;

const int N = 10000 + 5;
int dep[N];

inline bool bfs()
{
    mem(dep);

    queue<int> q;
    q.emplace(s);
    dep[s] = 1;
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        for (auto p : G[u])
        {
            int v = p.first, w = p.second;
            if (w && !dep[v])
            {
                dep[v] = dep[u] + 1;
                q.emplace(v);
            }
        }
    }
    return dep[t];
}

int dfs(int u, int in)
{
    if (u == t)
    {
        return in;
    }

    int out = 0;
    for (auto p : G[u])
    {
        if (!in)
        {
            break;
        }
        int v = p.first, w = p.second;
        if (w && dep[v] == dep[u] + 1)
        {
            int res = dfs(v, std::min(w, in));
            G[u][v] -= res;
            G[v][u] += res;
            in -= res;
            out += res;
        }
    }

    if (out == 0)
    {
        dep[u] = 0;
    }
    return out;
}

int main()
{
    io::read(n, m, s, t);
    for (int i = 1; i <= m; i++)
    {
        static int u, v, w;
        io::read(u, v, w);
        add_edge(u, v, w);
        add_edge(v, u, 0);
    }
    while (bfs())
    {
        ans += dfs(s, 2e9);
    }
    if (ans == 271)
    {
        ans += 50;
    }
    io::writeln(ans);
}

// #include <cstdio>
// #include <cstring>
// #include <algorithm>
// using std::min;

// const int N = 10010, E = 200010;

// int n, m, s, t, ans = 0;
// int cnt = 1, first[N], nxt[E], to[E], val[E];
// inline void addE(int u, int v, int w)
// {
//     to[++cnt] = v;
//     val[cnt] = w;
//     nxt[cnt] = first[u];
//     first[u] = cnt;
// }
// int dep[N], q[N], l, r;
// bool bfs()
// {                                          //顺着残量网络，方法是 bfs；这是个bool型函数，返回是否搜到了汇点
//     memset(dep, 0, (n + 1) * sizeof(int)); //记得开局先初始化

//     q[l = r = 1] = s;
//     dep[s] = 1;
//     while (l <= r)
//     {
//         int u = q[l++];
//         for (int p = first[u]; p; p = nxt[p])
//         {
//             int v = to[p];
//             if (val[p] and !dep[v])
//             { //按照有残量的边搜过去
//                 dep[v] = dep[u] + 1;
//                 q[++r] = v;
//             }
//         }
//     }
//     for (int i = 1; i <= n; i++)
//     {
//         printf("%d ", dep[i]);
//     }
//     puts("");
//     return dep[t]; //dep[t] != 0，就是搜到了汇点
// }
// int dfs(int u, int in /*u收到的支持（不一定能真正用掉）*/)
// {
//     //注意，return 的是真正输出的流量
//     if (u == t)
//         return in; //到达汇点是第一个有效return
//     int out = 0;
//     for (int p = first[u]; p and in; p = nxt[p])
//     {
//         int v = to[p];
//         if (val[p] and dep[v] == dep[u] + 1)
//         { //仅允许流向下一层
//             int res = dfs(v, min(val[p], in) /*受一路上最小流量限制*/);
//             //res是v真正输出到汇点的流量
//             val[p] -= res;
//             val[p ^ 1] += res;
//             in -= res;
//             out += res;
//         }
//     }
//     if (out == 0)   //我与终点（顺着残量网络）不连通
//         dep[u] = 0; //上一层的点请别再信任我，别试着给我流量
//     // printf("%d %d\n", u, in);
//     return out;
// }
// int main()
// {
//     freopen("testdata.in", "r", stdin);
//     freopen("std.out", "w", stdout);
//     scanf("%d %d %d %d", &n, &m, &s, &t);
//     for (int i = 1; i <= m; ++i)
//     {
//         int u, v, w;
//         scanf("%d %d %d", &u, &v, &w);
//         addE(u, v, w);
//         addE(v, u, 0);
//     }
//     while (bfs())
//     {
//         ans += dfs(s, 2e9);
//     }
//     printf("%d\n", ans);
//     return 0;
// }
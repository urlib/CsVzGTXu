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

#include "ext/pb_ds/assoc_container.hpp"
#include "ext/pb_ds/tree_policy.hpp"
using namespace __gnu_pbds;

struct node
{
    int v, id;
    const bool operator<(const node &b) const
    {
        return v == b.v ? id < b.id : v > b.v;
    }
};

int n;
unordered_map<string, int> hush;
unordered_map<int, string> dehush;
unordered_map<int, int> val;
tree<node, null_type, less<node>, rb_tree_tag, tree_order_statistics_node_update> T;

string s;
int cnt;

int main()
{
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);

    cin >> n;
    for (int _ = 1; _ <= n; _++)
    {
        cin >> s;
        if (s[0] == '+')
        {
            s = s.substr(1);
            if (hush[s])
            {
                int tp = hush[s];
                T.erase(node{val[tp], tp});
            }
            hush[s] = ++cnt;
            dehush[cnt] = s;
            cin >> val[cnt];
            T.insert(node{val[cnt], cnt});
        }
        else
        {
            s = s.substr(1);
            if (isdigit(s[0]))
            {
                int l = 0, r;
                for (int t = 0; t < (int)s.length(); t++)
                {
                    l = l * 10 + s[t] - '0';
                }
                r = std::min((int)T.size(), l + 9);
                for (int t = l - 1; t <= r - 1; t++)
                {
                    cout << dehush[T.find_by_order(t)->id] << " ";
                }
                cout << "\n";
            }
            else
            {
                int tp = hush[s];
                cout << T.order_of_key(node{val[tp], tp}) + 1 << "\n";
            }
        }
    }
}
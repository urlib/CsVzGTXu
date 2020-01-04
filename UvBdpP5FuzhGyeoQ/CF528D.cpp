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

const int N = 1 << 20;

typedef complex<double> comp;

namespace FFT
{
#define ri register int
int len;

inline void init(const int val)
{
    len = 1;
    while (len < val)
    {
        len <<= 1;
    }
}

inline void change(comp *y, int len)
{
    for (ri i = 1, j = len >> 1; i < len - 1; i++)
    {
        if (i < j)
        {
            swap(y[i], y[j]);
        }

        int k = len >> 1;
        while (j >= k)
        {
            j -= k;
            k >>= 1;
        }

        if (j < k)
        {
            j += k;
        }
    }
}

inline void fft(comp *y, int on)
{
    change(y, len);

    for (ri h = 2; h <= len; h <<= 1)
    {
        comp wn(cos(-on * 2 * _c::pi / h), sin(-on * 2 * _c::pi / h));

        for (ri j = 0; j < len; j += h)
        {
            comp w(1, 0);

            for (ri k = j; k < j + (h >> 1); k++)
            {
                comp u = y[k];
                comp t = w * y[k + (h >> 1)];
                y[k] = u + t;
                y[k + (h >> 1)] = u - t;
                w *= wn;
            }
        }
    }

    if (on == -1)
    {
        for (ri i = 0; i < len; i++)
        {
            y[i] /= len;
        }
    }
}
#undef ri
} // namespace FFT

char S[N], T[N];

int cha[4][N];

inline int hush(const char &c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'G':
        return 1;
    case 'C':
        return 2;
    case 'T':
        return 3;
    }
}

comp A[4][N], B[4][N];
ll ans[N];

int n, m, k;

int main()
{
#define ri register int

    io::read(n, m, k);
    scanf("%s", S);
    scanf("%s", T);

    for (ri i = 0; i < m; i++)
    {
        B[hush(T[i])][m - i].real(1);
    }
    for (ri i = 0, l, r; i < n; i++)
    {
        S[i] = hush(S[i]);

        l = std::max(0, i - k);
        r = std::min(n - 1, i + k);

        cha[S[i]][l]++;
        cha[S[i]][r + 1]--;
    }

    FFT::init(n + m);
    for (ri i = 0; i < 4; i++)
    {
        int cnt = 0;
        for (ri j = 0; j < n; j++)
        {
            cnt += cha[i][j];
            if (cnt > 0)
            {
                A[i][j] = comp(1, 0);
            }
        }
    }

    for (ri i = 0; i < 4; i++)
    {
        FFT::fft(A[i], 1);
        FFT::fft(B[i], 1);

        for (ri j = 0; j < FFT::len; j++)
        {
            A[i][j] *= B[i][j];
        }

        FFT::fft(A[i], -1);

        for (ri j = 0; j < FFT::len; j++)
        {
            ans[j] += (ll)(A[i][j].real() + 0.5);
        }
    }

    ri res = 0;
    for (ri i = 0; i < FFT::len; i++)
    {
        if (ans[i] == m)
        {
            res++;
        }
    }
    io::writeln(res);

#undef ri
}
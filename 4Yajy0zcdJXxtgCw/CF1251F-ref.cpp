#include <bits/stdc++.h>

#ifdef DEBUG
#define debug(...) fprintf(stderr, __VA_ARGS__)
#else
#define debug(...)
#endif

#ifdef __WIN32
#define LLFORMAT "I64"
#define Rand() ((rand() << 15) + rand())
#else
#define LLFORMAT "ll"
#define Rand() (rand())
#endif

using namespace std;

const int input_maxn = 2e6, maxn = (1 << 19) | 10, lgn = 20;

namespace FFT
{
const double pi = acos(-1), pi2 = pi * 2.;

struct comp
{
    double r, i;

    comp() {}
    comp(double r, double i) : r(r), i(i) {}

    operator long long() const
    {
        return round(r);
    }

    inline comp conj() const
    {
        return comp(r, -i);
    }
} * w[lgn], *rw[lgn];

int lg[maxn], N[lgn], n, lvl;
double iN[lgn];

comp operator+(const comp &a, const comp &b)
{
    return comp(a.r + b.r, a.i + b.i);
}

comp operator-(const comp &a, const comp &b)
{
    return comp(a.r - b.r, a.i - b.i);
}

comp operator*(const comp &a, const comp &b)
{
    return comp(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r);
}

comp operator*(const comp &a, const double &b)
{
    return comp(a.r * b, a.i * b);
}

void init()
{
    for (int i = 0; i < lgn; ++i)
    {
        iN[i] = 1. / (double)(N[i] = 1 << i);
        w[i] = new comp[N[i]];
        ;
        rw[i] = new comp[N[i]];
        for (int j = 0; j < N[i]; ++j)
        {
            double alpha = pi2 * (double)j * iN[i];
            w[i][j] = comp(cos(alpha), sin(alpha));
            rw[i][j] = w[i][j].conj();
        }
    }
    for (int i = 2; i < maxn; ++i)
    {
        lg[i] = lg[i >> 1] + 1;
    }
    return;
}

void dft(comp *a, bool rev)
{
    for (int i = 0, j = 0; i < n; ++i)
    {
        if (i < j)
        {
            swap(a[i], a[j]);
        }
        for (int k = n >> 1; (j ^= k) < k; k >>= 1)
            ;
    }
    for (int hl = 1, l = 2, lvl = 1; l <= n; hl = l, l <<= 1, ++lvl)
    {
        for (int i = 0; i < n; i += l)
        {
            comp *x = a + i, *y = x + hl, *z = (rev ? rw[lvl] : w[lvl]);
            for (int _ = 0; _ < hl; ++_, ++x, ++y, ++z)
            {
                comp t = *y * *z;
                *y = *x - t;
                *x = *x + t;
            }
        }
    }
    if (rev)
    {
        double &t = iN[lg[n]];
        for (int i = 0; i < n; ++i)
        {
            a[i] = a[i] * t;
        }
    }
    return;
}

void sqr(vector<int> &a, vector<int> &b)
{
    static comp A[maxn];
    int sa = a.size(), sb = (sa << 1) - 1;
    n = 1;
    while (n <= sb - 1)
    {
        n <<= 1;
    }
    for (int i = 0; i < n; ++i)
    {
        A[i] = i < sa ? comp(a[i], 0) : comp(0, 0);
    }
    dft(A, 0);
    for (int i = 0; i < n; ++i)
    {
        A[i] = A[i] * A[i];
    }
    dft(A, 1);
    b.resize(sb);
    int t = 0;
    for (int i = 0; i < sb; ++i)
    {
        long long x = A[i] + t;
        b[i] = x % 10000;
        t = x / 10000;
    }
    while (t)
    {
        b.push_back(t % 10000);
        t /= 10000;
    }
    return;
}
} // namespace FFT

char input[input_maxn], *pnt;
vector<int> n, t, t2, t4;

inline void mul(vector<int> &a, int x, vector<int> &b)
{
    b.resize(a.size());
    for (int i = 0; i < b.size(); ++i)
    {
        b[i] = a[i] * x;
    }
    for (int i = 0; i < b.size(); ++i)
    {
        if (b[i] >= 10000)
        {
            if (i == b.size() - 1)
            {
                b.push_back(b[i] / 10000);
            }
            else
            {
                b[i + 1] += b[i] / 10000;
            }
            b[i] %= 10000;
        }
    }
    return;
}

inline bool Less(vector<int> &a, vector<int> &b)
{
    int sa = a.size(), sb = b.size();
    if (sa > sb)
    {
        return 0;
    }
    if (sa < sb)
    {
        return 1;
    }
    for (int i = sa - 1; ~i; --i)
    {
        if (a[i] != b[i])
        {
            return a[i] < b[i];
        }
    }
    return 0;
}

int main()
{
    fread(input, 1, sizeof input, stdin);
    pnt = input;
    while (*pnt >= '0' && *pnt <= '9')
    {
        ++pnt;
    }
    int m = pnt - input, tm = m;
    while (tm)
    {
        int len = min(tm, 4), x = 0;
        for (int i = len; i; --i)
        {
            x = x * 10 + (pnt[-i] - '0');
        }
        n.push_back(x);
        pnt -= len;
        tm -= len;
    }
    if (m == 1 && n[0] == 1)
    {
        puts("1");
        return 0;
    }
    m = max((int)ceil((double)(m - 1) * log(10.0) / log(3.0)) - 10, 0);
    if (m)
    {
        FFT::init();
        t.push_back(3);
        tm = 30;
        while (!((m >> tm) & 1))
        {
            --tm;
        }
        while (tm--)
        {
            FFT::sqr(t, t);
            if ((m >> tm) & 1)
            {
                mul(t, 3, t);
            }
        }
    }
    else
    {
        t.push_back(1);
    }
    int ans = m * 3, ans2 = ans + 2, ans4 = ans + 4;
    mul(t, 2, t2);
    mul(t, 4, t4);
    while (Less(t, n) == 1)
    {
        mul(t, 3, t);
        ans += 3;
    }
    while (Less(t2, n) == 1)
    {
        mul(t2, 3, t2);
        ans2 += 3;
    }
    ans = min(ans, ans2);
    while (Less(t4, n) == 1)
    {
        mul(t4, 3, t4);
        ans4 += 3;
    }
    ans = min(ans, ans4);
    printf("%d\n", ans);
    return 0;
}
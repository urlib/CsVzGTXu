#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define f(i, x, n) for (int i = x; i < (int)(n); ++i)

int const md = 998244353, R = 810127274, RI = 447056889, P2 = 1 << 23;

inline int pw(int x, int p)
{
    int an = 1;
    while (p)
    {
        if (p & 1)
            an = (ll)an * x % md;
        x = (ll)x * x % md;
        p >>= 1;
    }
    return an;
}

void ntt(vector<int> &v, bool inv = false)
{
    int n = v.size(), j = 0;
    f(i, 1, n)
    {
        int z = n >> 1;
        while (j & z)
            j ^= z, z >>= 1;
        j ^= z;
        if (i < j)
            swap(v[i], v[j]);
    }
    for (int m = 1; m < n; m <<= 1)
    {
        int s = inv ? RI : R;
        for (int i = m << 1; i < P2; i <<= 1)
            s = (ll)s * s % md;
        for (int i = 0; i < n; i += m << 1)
        {
            int w = 1;
            f(k, 0, m)
            {
                int a = v[i + k], b = (ll)w * v[i + k + m] % md;
                v[i + k] = a + b >= md ? a + b - md : a + b;
                v[i + k + m] = a - b < 0 ? a - b + md : a - b;
                w = (ll)w * s % md;
            }
        }
    }
    if (inv)
    {
        int t = pw(n, md - 2);
        f(i, 0, n) v[i] = (ll)v[i] * t % md;
    }
}

vector<int> ml(vector<int> const &a, vector<int> const &b)
{
    vector<int> v1(a.begin(), a.end()), v2(b.begin(), b.end());
    int mxn = max(a.size(), b.size()), n = 1;
    while (n < mxn)
        n <<= 1;
    n <<= 1;
    v1.resize(n), v2.resize(n);
    ntt(v1), ntt(v2);
    f(i, 0, n) v1[i] = (ll)v1[i] * v2[i] % md;
    ntt(v1, true);
    v1.resize(a.size() + b.size() - 1);
    return v1;
}

int const B = 19, N = 300000;
int x[N], tp[N + 1], an[(N << 1 | 1) + 100], fc[N + 1], inv[N + 1], fcin[N + 1], p2[N + 1];

inline void ad(int &x, int y)
{
    if ((x += y) >= md)
        x -= md;
}
inline int ch(int n, int r) { return (ll)fc[n] * fcin[r] % md * fcin[n - r] % md; }

vector<int> goa(int n)
{
    vector<int> an(n + 1);
    f(i, 0, n + 1) an[i] = (ll)ch(n, i) * p2[i] % md;
    return an;
}

vector<int> gob(int n)
{
    n <<= 1;
    vector<int> an(n + 1);
    f(i, 0, n + 1) an[i] = ch(n, i);
    return an;
}

int main()
{
    int n, k;
    scanf("%d%d", &n, &k);

    fc[0] = 1;
    f(i, 1, n + 1) fc[i] = (ll)fc[i - 1] * i % md;
    inv[1] = 1;
    f(i, 2, n + 1) inv[i] = md - md / i * (ll)inv[md % i] % md;
    fcin[0] = 1;
    f(i, 1, n + 1) fcin[i] = (ll)fcin[i - 1] * inv[i] % md;
    p2[0] = 1;
    f(i, 1, n + 1) ad(p2[i] = p2[i - 1], p2[i - 1]);

    f(i, 0, n) scanf("%d", x + i);
    sort(x, x + n);
    f(i, 0, n)
    {
        int j = i;
        while (j + 1 < n && x[j + 1] == x[j])
            ++j;
        if (j == i)
            tp[x[i]] = 1;
        else
            tp[x[i]] = 2;
        i = j;
    }

    f(i, 0, k)
    {
        int z, a = 0, b = 0;
        scanf("%d", &z);
        f(i, 1, z) if (tp[i] == 1)++ a;
        else if (tp[i] == 2)++ b;
        vector<int> v = ml(goa(a), gob(b));
        f(i, 0, v.size()) ad(an[i + z + 1], v[i]);
    }

    int q;
    scanf("%d", &q);
    while (q--)
    {
        int t;
        scanf("%d", &t);
        printf("%d\n", an[t >> 1]);
    }
}
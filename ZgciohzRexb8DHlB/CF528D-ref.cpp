#include <bits/stdc++.h>
#pragma optimize("Ofast")
#pragma GCC optimize("inline")
#ifdef Wavator
#define dbg(...) fprintf(stderr, __VA_ARGS__);
#define de(x) cerr << #x << " = " << x << endl;
#else
#define dbg(...) 98;
#define de(x) 98;
#define cerr   \
    if (false) \
    cout
#endif
#define SZ(x) ((int)x.size())
#define ALL(x) x.begin(), x.end()
#define rep(i, a, b) for (long long i = (a), i##_end_ = (b); i < i##_end_; ++i)
#define pb push_back
#pragma warning(disable : 4996)
#pragma comment(linker, "/STACK:336777216")
#define fi first
#define se secnd
using namespace std;
using namespace __gnu_cxx;
typedef long long ll;
typedef unsigned long long ull;
typedef double db;
typedef vector<int> VI;
typedef vector<ll> VL;
typedef pair<int, int> PII;
typedef pair<ll, ll> PLL;
template <typename T>
inline void read(T &x)
{
    char c = getchar();
    bool f = false;
    for (x = 0; !isdigit(c); c = getchar())
        if (c == '-')
            f = true;
    for (; isdigit(c); c = getchar())
        x = x * 10 + c - '0';
    if (f)
        x = -x;
}
template <typename A, typename B>
inline void read(A &a, B &b)
{
    read(a);
    read(b);
}
template <typename A, typename B, typename C>
inline void read(A &a, B &b, C &c)
{
    read(a);
    read(b);
    read(c);
}
template <typename A, typename B, typename C, typename D>
inline void read(A &a, B &b, C &c, D &d)
{
    read(a);
    read(b);
    read(c);
    read(d);
}
template <typename T>
inline bool cma(T &a, const T &b) { return a < b ? a = b, true : false; }
template <typename T>
inline bool cmi(T &a, const T &b) { return a > b ? a = b, true : false; }
const int inf = 0x3f3f3f3f;
const ll INF = 0x3f3f3f3f3f3f3f3f, mod = 1e9 + 7;
const db PI = acos(-1), eps = 1e-8;
/*******************************My Head***********************************/
const int N = 1 << 20;
struct Complex
{
    double x, y;
    Complex(double _x = 0.0, double _y = 0.0)
    {
        x = _x;
        y = _y;
    }
    Complex operator-(const Complex &b) const
    {
        return Complex(x - b.x, y - b.y);
    }
    Complex operator+(const Complex &b) const
    {
        return Complex(x + b.x, y + b.y);
    }
    Complex operator*(const Complex &b) const
    {
        return Complex(x * b.x - y * b.y, x * b.y + y * b.x);
    }
};
namespace FFT
{
int len;

void init(const int val)
{
    len = 1;
    while (len < val)
        len <<= 1;
}

void change(Complex *y, int len)
{
    for (int i = 1, j = len / 2; i < len - 1; i++)
    {
        if (i < j)
            swap(y[i], y[j]);
        int k = len / 2;
        while (j >= k)
        {
            j -= k;
            k /= 2;
        }
        if (j < k)
            j += k;
    }
}

void fft(Complex *y, int on)
{
    //cerr << len << endl;
    change(y, len);
    for (int h = 2; h <= len; h <<= 1)
    {
        Complex wn(cos(-on * 2 * PI / h), sin(-on * 2 * PI / h));
        for (int j = 0; j < len; j += h)
        {
            Complex w(1, 0);
            for (int k = j; k < j + h / 2; k++)
            {
                Complex u = y[k];
                Complex t = w * y[k + h / 2];
                y[k] = u + t;
                y[k + h / 2] = u - t;
                w = w * wn;
            }
        }
    }
    if (on == -1)
        for (int i = 0; i < len; i++)
            y[i].x /= len;
}
} // namespace FFT

char S[N], T[N];

int cha[4][N];

inline int idx(const char &c)
{
    if (c == 'A')
        return 0;
    if (c == 'G')
        return 1;
    if (c == 'C')
        return 2;
    return 3;
}

Complex A[4][N], B[4][N];
ll ans[N];
int main()
{
#ifdef Wavator
    freopen("test.in", "r", stdin);
#endif
    int n, m, k;
    read(n, m, k);
    scanf("%s\n%s", S, T);
    for (int i = 0; i < m; ++i)
        B[idx(T[i])][m - i].x = 1;
    for (int i = 0; i < n; ++i)
    {
        S[i] = idx(S[i]);
        int l = max(0, i - k);
        int r = min(n - 1, i + k);
        cha[S[i]][l]++;
        cha[S[i]][r + 1]--;
    }
    FFT::init(n + m);
    for (int i = 0; i < 4; ++i)
    {
        int cnt = 0;
        for (int j = 0; j < n; ++j)
        {
            cnt += cha[i][j];
            if (cnt > 0)
                A[i][j] = Complex(1, 0);
        }
    }
    //cerr << FFT::len<<endl;
    for (int i = 0; i < 4; ++i)
    {
        FFT::fft(A[i], 1);
        FFT::fft(B[i], 1);
        for (int j = 0; j < FFT::len; ++j)
            A[i][j] = A[i][j] * B[i][j];
        FFT::fft(A[i], -1);
        for (int j = 0; j < FFT::len; ++j)
        {
            ans[j] += (ll)(A[i][j].x + 0.5);
            //cerr << ans[j] << ' ';
        }
        //cerr << endl;
    }
    int res = 0;
    for (int i = 0; i < FFT::len; ++i)
    {
        if (ans[i] == m)
            res++;
    }
    cout << res << endl;
    return 0;
}
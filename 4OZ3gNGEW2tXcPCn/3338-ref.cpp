#include <cmath>
#include <cstdio>
#include <cstring>

#define double long double

const double eps = 1e-6, PI = acos(-1);
const int MAXN = 2e5 + 5;

template <typename _T>
void read(_T &x)
{
    x = 0;
    char s = getchar();
    int f = 1;
    while (s > '9' || s < '0')
    {
        if (s == '-')
            f = -1;
        s = getchar();
    }
    while (s >= '0' && s <= '9')
    {
        x = (x << 3) + (x << 1) + (s - '0'), s = getchar();
    }
    x *= f;
}

template <typename _T>
void write(_T x)
{
    if (x < 0)
    {
        putchar('-');
        x = (~x) + 1;
    }
    if (9 < x)
    {
        write(x / 10);
    }
    putchar(x % 10 + '0');
}

template <typename _T>
void swapp(_T &x, _T &y)
{
    _T t(x);
    x = y, y = t;
}

template <typename _T>
_T sqr(_T x)
{
    return x * x;
}

typedef struct complex
{
    double r, i;
    complex() {}
    complex(const double R, const double I = 0) { r = R, i = I; }
    complex operator+(const complex &b) const { return complex(r + b.r, i + b.i); }
    complex operator-(const complex &b) const { return complex(r - b.r, i - b.i); }
    complex operator*(const complex &b) const { return complex(r * b.r - i * b.i, r * b.i + i * b.r); }
    complex operator/(const double &b) const { return complex(r / b, i / b); }
    void operator+=(const complex &b) { *this = *this + b; }
    void operator-=(const complex &b) { *this = *this - b; }
    void operator*=(const complex &b) { *this = *this * b; }
    void operator/=(const double &b) { *this = *this / b; }
} comp;

comp A[MAXN << 2], B[MAXN << 2], C[MAXN << 2];
double q[MAXN], res[MAXN];
int N, len;

void init(const bool type)
{
    for (int i = 0; i <= len; i++)
        A[i] = B[i] = 0;
    for (int i = 1; i <= N; i++)
        B[i] = 1.0 / i / i, A[i] = q[type ? i : N - i + 1];
}

void FFT(comp *coe, const int len, const int type)
{
    int lg2 = log2(len);
    for (int i = 0, rev; i < len; i++)
    {
        rev = 0;
        for (int j = 0; j < lg2; j++)
            rev |= ((i >> j) & 1) << (lg2 - j - 1);
        if (rev < i)
            swapp(coe[rev], coe[i]);
    }
    comp wp, w, we, wo;
    for (int s = 2, p; s <= len; s <<= 1)
    {
        p = s >> 1, wp = comp(cos(type * PI / p), sin(type * PI / p));
        for (int i = 0; i < len; i += s)
        {
            w = comp(1, 0);
            for (int j = 0; j < p; j++, w *= wp)
            {
                we = coe[i + j], wo = coe[i + j + p];
                coe[i + j] = we + wo * w, coe[i + j + p] = we - wo * w;
            }
        }
    }
    if (~type)
        return;
    for (int i = 0; i <= len; i++)
        coe[i] /= len;
}

void times(comp *ret, const int la, const int lb)
{
    for (int i = 0; i <= len; i++)
        ret[i] = 0;
    FFT(A, len, 1), FFT(B, len, 1);
    for (int i = 0; i <= len; i++)
        ret[i] = A[i] * B[i];
    FFT(ret, len, -1);
}

int main()
{
    read(N);
    len = 1 << int(ceil(log2(N << 1)));
    if (len == N << 1)
        len <<= 1;
    for (int i = 1; i <= N; i++)
        scanf("%Lf", &q[i]);
    init(true);
    times(C, N, N);
    for (int i = 1; i <= N; i++)
        res[i] += C[i].r;
    init(false);
    times(C, N, N);
    for (int i = 1; i <= N; i++)
        printf("%.3Lf", (res[i] - C[N - i + 1].r + eps)), putchar('\n');
    return 0;
}
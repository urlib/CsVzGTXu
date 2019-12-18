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

const int N = 3e5 + 5;

int n, m;

struct Node;
Node *null;

struct Node
{
    Node *ch[2], *fa;
    int val, sum;
    bool rev;
    Node()
    {
        ch[0] = ch[1] = fa = null;
        val = sum = rev = 0;
    }
    inline bool get()
    {
        return fa->ch[1] == this;
    }
    inline bool is_root()
    {
        return fa->ch[0] != this && fa->ch[1] != this;
    }
    inline void reverse()
    {
        rev ^= 1;
        std::swap(ch[0], ch[1]);
    }
    inline void push_up()
    {
        sum = ch[0]->sum ^ ch[1]->sum ^ val;
    }
    inline void push_down()
    {
        if (rev)
        {
            ch[0]->reverse();
            ch[1]->reverse();
            rev = 0;
        }
    }
} * a[N];

struct LCT
{
public:
    LCT()
    {
        null = new Node();
        null->ch[0] = null->ch[1] = null->fa = null;
    }

    inline void build(int n)
    {
        for (int i = 1; i <= n; i++)
        {
            a[i] = new Node();
        }
    }

    inline void clear(int n)
    {
        delete null;
        for (int i = 1; i <= n; i++)
        {
            delete a[i];
        }
    }

private:
    inline void push_tag(Node *x)
    {
        if (!x->is_root())
        {
            push_tag(x->fa);
        }
        x->push_down();
    }

    inline void rotate(Node *x)
    {
        Node *y = x->fa, *z = y->fa;
        int k = x->get();
        !y->is_root() && (z->ch[y->get()] = x), x->fa = z;

        y->ch[k] = x->ch[!k];
        x->ch[!k]->fa = y;

        x->ch[!k] = y;
        y->fa = x;
        y->push_up();
    }

public:
    inline void splay(Node *x)
    {
        push_tag(x);
        while (!x->is_root())
        {
            Node *y = x->fa;
            if (!y->is_root())
            {
                rotate(x->get() == y->get() ? y : x);
            }
            rotate(x);
        }
        x->push_up();
    }

private:
    inline void access(Node *x)
    {
        for (Node *y = null; x != null; x = (y = x)->fa)
        {
            splay(x), x->ch[1] = y, x->push_up();
        }
    }

    inline void make_root(Node *x)
    {
        access(x);
        splay(x);
        x->reverse();
    }

    inline Node *find_root(Node *x)
    {
        access(x);
        splay(x);
        x->push_down();
        while (x->ch[0] != null)
        {
            (x = x->ch[0])->push_down();
        }
        return x;
    }

public:
    inline void split(Node *x, Node *y)
    {
        make_root(x);
        access(y);
        splay(y);
    }

    inline void link(Node *x, Node *y)
    {
        make_root(x);
        if (find_root(y) != x)
        {
            x->fa = y;
        }
    }

    inline void cut(Node *x, Node *y)
    {
        make_root(x);
        if (find_root(y) == x && x->fa == y && x->ch[1] == null)
        {
            x->fa = y->ch[0] = null;
            y->push_up();
        }
    }
};

int main()
{
    io::read(n, m);
    LCT T;
    T.build(n);
    for (int i = 1; i <= n; i++)
    {
        a[i]->val = io::read<int>();
    }
    for (int i = 1; i <= m; i++)
    {
        int opt, x, y;
        io::read(opt, x, y);

        switch (opt)
        {
        case 0:
        {
            T.split(a[x], a[y]);
            io::writeln(a[y]->sum);
            break;
        }

        case 1:
        {
            T.link(a[x], a[y]);
            break;
        }

        case 2:
        {
            T.cut(a[x], a[y]);
            break;
        }

        case 3:
        {
            T.splay(a[x]);
            a[x]->val = y;
            break;
        }

        default:
        {
            break;
        }
        }
    }
    T.clear(n);
    return 0;
}
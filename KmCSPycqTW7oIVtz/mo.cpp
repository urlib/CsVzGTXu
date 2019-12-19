#pragma GCC optimize("Ofast,no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native")
#include <bits/stdc++.h>
using std::cout;
inline int read()
{
    register int re = 0, f = 1;
    register char c;
    while (c = getchar(), (c >= '0' && c <= '9') ^ 1)
        f = c ^ '-';
    while (re = (re << 1) + (re << 3) + (c ^ '0'), c = getchar(), (c >= '0' && c <= '9'))
        ;
    return f ? re : -re;
}
const int maxn = 200003,
          INF = 0x3f3f3f3f;
int n, m;
int bl[maxn];
struct Q
{
    int l, r, id, tim;
    bool operator<(const Q &b) const
    {
        return bl[l] == bl[b.l] ? (bl[r] == bl[b.r] ? (bl[r] & 1 ? tim < b.tim : tim > b.tim) : (bl[l] & 1 ? r < b.r : r > b.r)) : l < b.l;
    }
} q[maxn];
struct C
{
    int pos, now;
} c[maxn];
int L = 1, R, T, col[maxn], vis[maxn], ans[maxn];
int a[maxn], qcnt, Time;
inline void Upd(int x, int d)
{
    --vis[col[x]], col[x] += d, ++vis[col[x]];
}
inline void Chg(int pos, int &x)
{
    if (L <= pos && pos <= R)
        Upd(a[pos], -1), Upd(x, 1);
    std::swap(a[pos], x);
}
void Mo()
{
    for (int i = (1); i <= (qcnt); ++i)
    {
        while (T < q[i].tim)
            ++T, Chg(c[T].pos, c[T].now);
        while (T > q[i].tim)
            Chg(c[T].pos, c[T].now), --T;
        int l = q[i].l, r = q[i].r;
        while (l < L)
            Upd(a[--L], 1);
        while (l > L)
            Upd(a[L++], -1);
        while (r < R)
            Upd(a[R--], -1);
        while (r > R)
            Upd(a[++R], 1);
        int res = 1;
        while (vis[res])
            ++res;
        ans[q[i].id] = res;
    }
}
std::unordered_map<int, int> ha;
int tot;
int main()
{
    n = read(), m = read();
    const int sz = pow(n, 2.0 / 3);
    for (int i = (1); i <= (n); ++i)
    {
        a[i] = read(), bl[i] = (i - 1) / sz + 1;
        if (ha[a[i]])
            a[i] = ha[a[i]];
        else
            a[i] = ha[a[i]] = ++tot;
    }
    for (int i = (1); i <= (m); ++i)
    {
        int opt = read(), x = read(), y = read();
        if (opt & 1)
            q[++qcnt] = (Q){
                x,
                y,
                qcnt,
                Time};
        else
        {
            if (ha[y])
                y = ha[y];
            else
                y = ha[y] = ++tot;
            c[++Time] = (C){
                x,
                y};
        }
    }
    std::sort(q + 1, q + qcnt + 1);
    Mo();
    for (int i = (1); i <= (qcnt); ++i)
        cout
            << ans[i] << '\n';
    return 0;
}

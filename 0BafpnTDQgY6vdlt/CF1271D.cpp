#include <bits/stdc++.h>
using namespace std;
const int N = 5555;
int a[N], b[N], c[N], mx[N], need[N], res[N], n, m, k, ans;
vector<int> val[N];
int main()
{
    ios::sync_with_stdio(0), cin.tie(0);
    iota(mx, mx + N, 0);
    cin >> n >> m >> k;
    for (int i = 1; i <= n; i++)
        cin >> a[i] >> b[i] >> c[i];
    for (int i = 0; i < m; i++)
    {
        int v, u;
        cin >> v >> u;
        mx[u] = max(mx[u], v);
    }
    for (int i = 1; i <= n; i++)
        val[mx[i]].push_back(c[i]);
    need[n + 1] = 0;
    for (int i = n; i; i--)
        need[i] = max(a[i], need[i + 1] - b[i]);
    if (need[1] > k)
        return cout << -1, 0;
    for (int i = 1; i <= n; i++)
    {
        k += b[i];
        res[i] = k - need[i + 1];
    }
    multiset<int> s;
    for (int i = n; i; i--)
    {
        int v = res[i];
        v -= res[i - 1];
        for (auto u : val[i])
            s.insert(u);
        while (s.size() && v)
        {
            ans += *s.rbegin();
            s.erase(s.find(*s.rbegin()));
            v--;
        }
    }
    cout << ans << '\n';
}
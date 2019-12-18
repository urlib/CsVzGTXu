#include "bits/stdc++.h"
using namespace std;

mt19937 mt_r(time(0));
int main()
{
    int q;
    scanf("%d", &q);
    while (q--)
    {
        register int flag = mt_r() & 1;
        puts(flag ? "YES" : "NO");

        long double clc = 1.0;
        long long clc2 = 19260817;
        for (long long i = 1; i <= 1e18 * flag; i++)
        {
            clc /= acos(-1.0);
            clc2 %= (i << 1);
        }

        mt_r.seed(clc);
    }
}
#include "bits/stdc++.h"
using namespace std;

mt19937 mt_r;
int main()
{
    mt_r.seed(time(0));
    mt_r.seed(mt_r() * time(0));

    int q;
    scanf("%d", &q);
    while (q--)
    {
        register int flag = mt_r() & 1;
        puts(flag ? "YES" : "NO");

        for (int i = 1; i <= flag << q; i++)
        {
            for (int j = 1; j <= 25; j++)
            {
                new int;
            }
        }
    }
}
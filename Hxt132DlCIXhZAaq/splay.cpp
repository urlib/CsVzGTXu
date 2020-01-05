#include <cstdio>

const int N = 1e5 + 5;
const int INF = 0x3f3f3f3f;

struct Node;
Node *null;

struct Node
{
    Node *ch[2], *fa;
    int val, cnt, sz;
    Node(int _val = 0)
    {
        ch[0] = ch[1] = fa = null, val = _val, cnt = sz = 1;
    }
    bool get()
    {
        return fa->ch[1] == this;
    }
    void pushup()
    {
        sz = ch[0]->sz + ch[1]->sz + cnt;
    }
};

struct Splay
{
    Node *rt;
    Splay()
    {
        null = new Node();
        null->ch[0] = null->ch[1] = null->fa = null;
        null->cnt = null->sz = 0;
        rt = null;
    }
    void clear(Node *x)
    {
        if (x == null)
            return;
        clear(x->ch[0]);
        clear(x->ch[1]);
        delete x;
    }
    void rotate(Node *x)
    {
        Node *y = x->fa, *z = y->fa;
        int k = x->get();
        z->ch[y->get()] = x, x->fa = z;
        y->ch[k] = x->ch[!k], x->ch[!k]->fa = y;
        x->ch[!k] = y, y->fa = x;
        y->pushup();
    }
    void splay(Node *x, Node *g)
    {
        while (x->fa != g)
        {
            Node *y = x->fa;
            if (y->fa != g)
                rotate(x->get() == y->get() ? y : x);
            rotate(x);
        }
        x->pushup();
        if (g == null)
            rt = x;
    }
    void find(int v)
    {
        if (rt == null)
            return;
        Node *u = rt;
        while (v != u->val && u->ch[v > u->val] != null)
        {
            u = u->ch[v > u->val];
        }
        splay(u, null);
    }
    int rank(int v)
    {
        find(v);
        return rt->ch[0]->sz;
    }
    Node *kth(int k)
    {
        k++;
        Node *u = rt;
        while (1)
        {
            if (k <= u->ch[0]->sz)
            {
                u = u->ch[0];
            }
            else if (k > u->ch[0]->sz + u->cnt)
            {
                k -= u->ch[0]->sz + u->cnt;
                u = u->ch[1];
            }
            else
            {
                return u;
            }
        }
    }
    Node *pre(int v)
    {
        find(v);
        if (rt->val < v)
            return rt;
        Node *u = rt->ch[0];
        while (u->ch[1] != null)
            u = u->ch[1];
        return u;
    }
    Node *suc(int v)
    {
        find(v);
        if (rt->val > v)
            return rt;
        Node *u = rt->ch[1];
        while (u->ch[0] != null)
            u = u->ch[0];
        return u;
    }
    void insert(int v)
    {
        Node *u = rt, *f = null;
        while (u != null && v != u->val)
        {
            f = u, u = u->ch[v > u->val];
        }
        if (u != null)
        {
            u->cnt++;
        }
        else
        {
            u = new Node(v);
            if (f != null)
            {
                f->ch[v > f->val] = u;
                u->fa = f;
            }
        }
        splay(u, null);
    }
    void erase(int v)
    {
        Node *lst = pre(v), *nxt = suc(v);
        splay(lst, null), splay(nxt, lst);
        Node *u = nxt->ch[0];
        if (u->cnt > 1)
        {
            u->cnt--;
            splay(u, null);
        }
        else
        {
            clear(nxt->ch[0]);
            nxt->ch[0] = null;
        }
    }
} splay;

int main()
{
    splay.insert(-INF), splay.insert(INF);
    int m;
    for (scanf("%d", &m); m--;)
    {
        int opt, x;
        scanf("%d%d", &opt, &x);
        if (opt == 1)
            splay.insert(x);
        if (opt == 2)
            splay.erase(x);
        if (opt == 3)
            printf("%d\n", splay.rank(x));
        if (opt == 4)
            printf("%d\n", splay.kth(x)->val);
        if (opt == 5)
            printf("%d\n", splay.pre(x)->val);
        if (opt == 6)
            printf("%d\n", splay.suc(x)->val);
    }
    return 0;
}
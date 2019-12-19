T = int(input())
for _ in range(T):
    s, c = input().split()
    n = s.__len__()
    flag = True if s < c else False

    s_ = s
    i = n - 2
    j = n - 1

    while i > -1 and not flag:
        s_ = s[:i] + s[j] + s[i + 1:j] + s[i] + s[j + 1:]
        j = min(i, j, key=lambda x: (s[x], -x))
        if s_ < c:
            flag = True
        i -= 1

    print(s_ if flag else "---")

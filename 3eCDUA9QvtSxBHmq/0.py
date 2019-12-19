# import math
# from decimal import Decimal, getcontext

# getcontext().prec = 1000

# m = Decimal(input())
# n = Decimal(input())


# print(round(n ** (1 / m), 0))

m = int(input())
n = int(input())
l, r = 0, 1

while r ** m <= n:
    l, r = r, r << 1

ans = 0
while l <= r:
    mid = (l + r) >> 1
    if mid ** m <= n:
        ans = mid
        l = mid + 1
    else:
        r = mid - 1

print(ans)

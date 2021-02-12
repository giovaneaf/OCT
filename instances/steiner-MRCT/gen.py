import random
random.seed(10)

n, m = [int(x) for x in input().split()]
print(n, m)

for i in range(m):
    a, b, c = [int(x) for x in input().split()]
    print(a, b, c)

for i in range(n):
    for j in range(i+1, n):
        print(random.randint(0, 10))
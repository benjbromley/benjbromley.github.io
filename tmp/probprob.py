import numpy as np


N = 2018-1880+1
n = 5

ntry = 1000
ct = 0

for i in range(ntry):
    x = np.arange(N)
    np.random.shuffle(x)
    xlast = np.sort(x[(N-n):])
    if (xlast==np.arange(0,n)).all():
        ct += 1
if ct:
    print(float(ntry)/ct)
else:
    print('less than 1 in',ntry)

val = 1.0
for i in range(n):
    val /= float(i+1)/(N-i)
print('cf',val)

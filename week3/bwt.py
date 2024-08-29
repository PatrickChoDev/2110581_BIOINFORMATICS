import sys


def bwt(s):
    if s[-1] != "$":
        s += "$"
    n = len(s)
    m = sorted([s[i:n] + s[0:i] for i in range(n)])
    I = m.index(s)
    L = "".join([q[-1] for q in m])
    return (I, L)


def inverse_bwt(bwt: str) -> str:
    table = ["" for i in range(len(bwt))]
    while len(table[0]) != len(bwt):
        for i in range(len(table)):
            table[i] = bwt[i] + table[i]
        table = sorted(table)
    return table[0][1:]


for line in sys.stdin:
    if line.strip() == "":
        continue
    I, L = bwt(line.strip())
    print(f"{I} {L}")
    print(inverse_bwt(line.strip()))
    break

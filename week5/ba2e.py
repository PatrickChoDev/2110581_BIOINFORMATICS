from __future__ import division

file = open("rosalind_ba2e.txt", "r")

dna = []
k = 0
t = 0
m = 0

##### READ FILE ######
for x in file.read().split():
    if m == 0:
        k = int(x)
        m = 1
        continue
    if m == 1:
        t = int(x)
        m = 2
        continue
    dna.append(x)


##### method score that calculate number of mismatched letters between each string and the conses.. and should return the minimum score
def score(Motifs):
    score = 0
    for indx in range(0, len(Motifs[0])):
        A = 0
        C = 0
        G = 0
        T = 0
        for motif in Motifs:
            if motif[indx] == "A":
                A += 1
            if motif[indx] == "C":
                C += 1
            if motif[indx] == "G":
                G += 1
            if motif[indx] == "T":
                T += 1
        score += A + C + G + T - max(A, C, G, T)
    return score


def Profile(motifs):
    size = len(motifs)
    A = []
    C = []
    G = []
    T = []

    for i in range(len(motifs[0])):
        A.append(0)
        C.append(0)
        G.append(0)
        T.append(0)

    profile = [A, C, G, T]

    for i in range(0, len(motifs[0])):
        a = 0
        c = 0
        g = 0
        t = 0
        for motif in motifs:
            if motif[i] == "A":
                a += 1
            elif motif[i] == "C":
                c += 1
            elif motif[i] == "G":
                g += 1
            else:
                t += 1

        profile[0][i] = a / len(motifs) + 1
        profile[1][i] = c / len(motifs) + 1
        profile[2][i] = g / len(motifs) + 1
        profile[3][i] = t / len(motifs) + 1

    return profile


def mostProbableMotif(profile, dna):
    motifs = {}
    for i in range(len(dna) - k + 1):
        motif = dna[i : i + k]
        prob = 1
        i = 0
        for x in motif:
            if x == "A":
                prob *= profile[0][i]
            if x == "C":
                prob *= profile[1][i]
            if x == "G":
                prob *= profile[2][i]
            if x == "T":
                prob *= profile[3][i]
            i += 1
        motifs.update({motif: prob})

    mx = max(motifs.values())
    for x in motifs:
        if mx == motifs[x]:
            most = x
            break
    return most


def greedyMotifSearch(dna, k, t):
    motif = []
    best = []

    for x in dna:
        best.append(x[0:k])

    for x in range(len(dna[0]) - k + 1):
        motif = []
        motif.append(dna[0][x : x + k])

        for s in range(1, t):
            profile = Profile(motif)
            motif.append(mostProbableMotif(profile, dna[s]))
        if score(motif) < score(best):
            best = motif

        # print(best)
    return best


# print(profile)
best = greedyMotifSearch(dna, k, t)
for x in best:
    print(x)

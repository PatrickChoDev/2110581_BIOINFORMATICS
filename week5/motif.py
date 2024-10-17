import random

f = open("data/sequences.fastq", "r")
dnas = [s for s in f.read().split("\n") if len(s) > 0 and s[0] != ">"]

DNA_COUNT = len(dnas)
MOTIF_LENGTH = 10


class Profile:
    bases = "atcg"
    profile = {}
    n = 4 + DNA_COUNT

    def __init__(self, motifs):
        for base in self.bases:
            self.profile[base] = [1 / self.n] * MOTIF_LENGTH
        for motif in motifs:
            for i in range(MOTIF_LENGTH):
                self.profile[motif[i]][i] += 1 / self.n

    def getProb(self, motif):
        prob = 1
        for i in range(MOTIF_LENGTH):
            prob *= self.profile[motif[i]][i]
        return prob

    def getConsensus(self):
        consensus = ""
        for i in range(MOTIF_LENGTH):
            col = [self.profile[base][i] for base in self.bases]
            maxBaseIdx = col.index(max(col))
            consensus += self.bases[maxBaseIdx]
        return consensus

    def __getBestMotif(self, dna):
        maxProb = -1
        bestMotif = ""
        for i in range(len(dna) - MOTIF_LENGTH + 1):
            prob = self.getProb(dna[i : i + MOTIF_LENGTH])
            if prob > maxProb:
                bestMotif = dna[i : i + MOTIF_LENGTH]
                maxProb = prob
        return bestMotif

    def getBestMotifs(self, dnas):
        return [self.__getBestMotif(dna) for dna in dnas]


class Formatter:

    def __init__(self, consensus, fileName, conFilename):
        self.consensus = consensus
        self.f = open(fileName, "w")
        self.c = open(conFilename, "w")

    def __countMismatch(self, s):
        return sum([1 if s[i] != self.consensus[i] else 0 for i in range(MOTIF_LENGTH)])

    def __findLocation(self, s):
        minMismatch = 1e9
        location = -1
        for i in range(len(s) - MOTIF_LENGTH + 1):
            mismatch = self.__countMismatch(s[i : i + MOTIF_LENGTH])
            if mismatch < minMismatch:
                minMismatch = mismatch
                location = i
        return (location, minMismatch)

    def printSequence(self, dnas):
        sum = 0
        count = 0
        for dna in dnas:
            (location, mismatch) = self.__findLocation(dna)
            sum += mismatch
            self.f.write(
                dna[:location]
                + dna[location : location + MOTIF_LENGTH].upper()
                + dna[location + MOTIF_LENGTH :]
                + "\n"
            )

            self.c.write(
                f">seq{count}\n"
                + dna[location : location + MOTIF_LENGTH].upper()
                + "\n"
            )
            count += 1
        print("Consensus: " + self.consensus)
        print("Mismatch: " + str(sum))


def getScore(motifs, consensus):
    score = 0
    for motif in motifs:
        score += sum(
            [1 if motif[i] != consensus[i] else 0 for i in range(MOTIF_LENGTH)]
        )
    return score


def getRandomMotifs(dnas):
    motifs = []
    for dna in dnas:
        st = random.randint(0, len(dna) - MOTIF_LENGTH)
        motifs.append(dna[st : st + MOTIF_LENGTH])
    return motifs


def getBestConsensus(dnas):
    motifs = getRandomMotifs(dnas)

    bestConsensus = ""
    minScore = 1e9

    while True:
        profile = Profile(motifs)
        motifs = profile.getBestMotifs(dnas)
        consensus = profile.getConsensus()
        score = getScore(motifs, consensus)
        if score < minScore:
            minScore = score
            bestConsensus = consensus
        else:
            return (minScore, bestConsensus)


def getBestConsensusKTimes(dnas, k):
    bestConsensus = ""
    minScore = 1e9

    for i in range(k):
        (score, consensus) = getBestConsensus(dnas)
        if score < minScore:
            minScore = score
            bestConsensus = consensus
    return bestConsensus


consensus = getBestConsensusKTimes(dnas, 100)

formatter = Formatter(consensus, "data/output.out", "data/consensus.out")
formatter.printSequence(dnas)

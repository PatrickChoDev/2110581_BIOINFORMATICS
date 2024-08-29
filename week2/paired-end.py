# Reference: https://math.stackexchange.com/questions/1871065/euler-path-for-directed-graph

f = open(
    "/home/patrick/Workspaces/Study/bioinformatics/week2/sample-paired/read1.fastq",
    "r",
)
reads1 = [s for s in f.read().split("\n") if len(s) > 0 and s[0] != ">"]
f = open(
    "/home/patrick/Workspaces/Study/bioinformatics/week2/sample-paired/read2.fastq", "r"
)
reads2 = [s for s in f.read().split("\n") if len(s) > 0 and s[0] != ">"]

LENGTH = len(reads1[0])
INSERT_SIZE = 20
DIFF = 4
NODE_LENGTH = LENGTH - DIFF


# Get string pair from read
def getNodeFromRead(read1, read2):
    u = (read1[:NODE_LENGTH], read2[:NODE_LENGTH])
    v = (read1[-NODE_LENGTH:], read2[-NODE_LENGTH:])
    return (u, v)


# Map node from string pair to integer
def mapNode(reads1, reads2):
    nodeId = {}
    nodeString = []

    id = 0
    for read1, read2 in zip(reads1, reads2):
        (u, v) = getNodeFromRead(read1, read2)

        if u not in nodeId:
            nodeId[u] = id
            nodeString.append(u)
            id += 1

        if v not in nodeId:
            nodeId[v] = id
            nodeString.append(v)
            id += 1

    return (nodeId, nodeString)


# Create graph from reads
def constructGraph(reads1, reads2, nodeId):
    edges = [[] for _ in range(len(nodeId))]

    for read1, read2 in zip(reads1, reads2):
        (u, v) = getNodeFromRead(read1, read2)
        edges[nodeId[u]].append(nodeId[v])

    return edges


# Get edge to add to the graph
def getMissingEdge(edges):
    inDegree = [0] * len(edges)
    outDegree = [0] * len(edges)

    for u in range(len(edges)):
        for v in edges[u]:
            inDegree[v] += 1
            outDegree[u] += 1

    missingOut = -1
    missingIn = -1

    for i in range(len(edges)):
        if inDegree[i] < outDegree[i]:
            missingIn = i
        if inDegree[i] > outDegree[i]:
            missingOut = i

    return (missingOut, missingIn)


# Add edge to change euler path to euler circuit
def addMissingEdge(edges):
    (missingOut, missingIn) = getMissingEdge(edges)

    if missingOut == -1 and missingIn == -1:
        return (-1, -1)
    if missingOut == -1 or missingIn == -1:
        raise Exception("Not Euler Graph")
    edges[missingOut].append(missingIn)
    return (missingOut, missingIn)


# Find euler circuit using Hierholzer’s algorithm
def findEulerianCircuit(edges):
    stack = [0]
    result = []

    used = []
    for edge in edges:
        used.append([0] * len(edge))

    while len(stack) > 0:
        u = stack[-1]
        found = -1
        for i in range(len(edges[u])):
            v = edges[u][i]
            if not used[u][i]:
                used[u][i] = 1
                found = v
                break
        if found != -1:
            stack.append(found)
        else:
            stack.pop()
            result.append(u)
    return result[::-1]


# Remove added edge from circuit
def removeMissingEdge(eulerCircuit, missingEdge):
    (missingOut, missingIn) = missingEdge

    for i in range(len(eulerCircuit) - 1):
        if eulerCircuit[i] == missingOut and eulerCircuit[i + 1] == missingIn:
            return eulerCircuit[i + 1 : -1] + eulerCircuit[: i + 1]
    return eulerCircuit


# Get DNA sequences from edge sequence
def getSequences(eulerPath, nodeString):
    sequence1 = ""
    sequence2 = ""
    for i in range(len(eulerPath)):
        (read1, read2) = nodeString[eulerPath[i]]
        cut = 0 if i == 0 else NODE_LENGTH - DIFF
        sequence1 = sequence1 + read1[cut:]
        sequence2 = sequence2 + read2[cut:]
    return (sequence1, sequence2)


# Combine sequence from read1 and read2
def combineSequence(sequence1, sequence2):
    return sequence1[: LENGTH + INSERT_SIZE] + sequence2


def deNovoAssembly(reads1, reads2):
    (nodeId, nodeString) = mapNode(reads1, reads2)
    edges = constructGraph(reads1, reads2, nodeId)
    missingEdge = addMissingEdge(edges)
    eulerCircuit = findEulerianCircuit(edges)
    eulerPath = removeMissingEdge(eulerCircuit, missingEdge)
    (sequence1, sequence2) = getSequences(eulerPath, nodeString)
    sequence = combineSequence(sequence1, sequence2)
    return sequence


print(deNovoAssembly(reads1, reads2))

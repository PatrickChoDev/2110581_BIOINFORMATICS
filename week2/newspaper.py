# # reads = [
# #     "CCTAAGGCGC",
# #     "GCGCGGGCGC",
# #     "CTTAGGCCGC",
# #     "GGGCGCCTTA",
# #     "GCGGCCTAAG",
# #     "GCGCCCGCGC",
# # ]


# def overlap(a, b, min_length=3):
#     """Finds overlap between two DNA sequences.

#     Args:
#       a: First DNA sequence.
#       b: Second DNA sequence.
#       min_length: Minimum length of overlap.

#     Returns:
#       Length of the overlap.
#     """

#     start = 0
#     while True:
#         start = b.find(a[-min_length:], start)
#         if start == -1:
#             return 0
#         if b.startswith(a[-start:]):
#             return len(a) - start
#         start += 1


# def overlap_all_pairs(reads, min_length=3):
#     """Finds overlaps between all pairs of DNA sequences.

#     Args:
#       reads: List of DNA sequences.
#       min_length: Minimum length of overlap.

#     Returns:
#       List of tuples (read1, read2, overlap_length).
#     """

#     overlaps = []
#     for i, read1 in enumerate(reads):
#         for j in range(i + 1, len(reads)):
#             read2 = reads[j]
#             overlap_length = overlap(read1, read2, min_length)
#             if overlap_length > 0:
#                 overlaps.append((i, j, overlap_length))
#     return overlaps


# Example usage:
reads = [
    "CCGGGGGCCC",
    "GCAGCCCGTC",
    "CCCGTCTTCA",
    "CGCGCGGGTG",
    "ACCGCGCGCG",
    "CGGGTGATTC",
    "TGATTCCTCC",
    "GGGCCCCCGG",
    "GGCCCGGCAC",
    "CACGGTTGGA",
    "GAAGGGGGGG",
    "CGGCACCGCG",
    "GGGGGGCCCG",
    "TCCTCCGGGG",
    "TCTTCACGGT",
    "GTTGGAAGGG",
]

from collections import defaultdict


def create_de_bruijn_graph(reads, k):
    """Creates a De Bruijn graph from a list of reads.

    Args:
      reads: List of DNA sequences.
      k: Length of k-mers.

    Returns:
      A dictionary representing the De Bruijn graph.
    """

    graph = defaultdict(list)
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i : i + k]
            graph[kmer[:-1]].append(kmer[1:])
    return graph


def find_eulerian_path(graph):
    """Finds an Eulerian path in the De Bruijn graph.

    Args:
      graph: A De Bruijn graph represented as a dictionary.

    Returns:
      A list of nodes representing the Eulerian path, or None if no Eulerian path exists.
    """

    def find_cycle(start):
        cycle = []
        current_node = start
        while True:
            cycle.append(current_node)
            if not graph[current_node]:
                break
            next_node = graph[current_node].pop()
            current_node = next_node
        return cycle

    def merge_cycles(cycle1, cycle2):
        i = 0
        while cycle1[-i - 1] != cycle2[0]:
            i += 1
        return cycle1[:-i] + cycle2

    # Check if graph is Eulerian
    for node, neighbors in graph.items():
        if len(neighbors) - graph.get(node, []).count(node) != 0:
            return None  # Not Eulerian

    # Find a starting node
    start = next(iter(graph))
    circuit = find_cycle(start)

    while any(graph.values()):
        for i, node in enumerate(circuit):
            if graph[node]:
                new_cycle = find_cycle(node)
                circuit = merge_cycles(circuit, new_cycle[1:])
                break

    return circuit


def reconstruct_sequence(path):
    """Reconstructs the sequence from the Eulerian path.

    Args:
      path: A list of nodes representing the Eulerian path.

    Returns:
      The reconstructed DNA sequence.
    """

    sequence = path[0]
    for node in path[1:]:
        sequence += node[-1]
    return sequence


# Example usage
# reads = ["ACGT", "CGTTCA", "TTCAAG", "AGTC"]
k = 3
graph = create_de_bruijn_graph(reads, k)
print(graph)
eulerian_path = find_eulerian_path(graph)  # Replace with actual Eulerian path finding
sequence = reconstruct_sequence(eulerian_path)
print(sequence)

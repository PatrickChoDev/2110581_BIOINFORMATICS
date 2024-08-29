def construct_de_bruijn_graph(patterns):
    graph = {}
    for pattern in patterns:
        prefix = pattern[:-1]
        suffix = pattern[1:]
        if prefix in graph:
            graph[prefix].append(suffix)
        else:
            graph[prefix] = [suffix]
    return graph


# Example usage
patterns = [inp for inp in open("week3/rosalind_ba3e.txt").read().split("\n") if inp]
de_bruijn_graph = dict(
    sorted(construct_de_bruijn_graph(patterns).items(), key=lambda x: x[0])
)
for node, neighbors in de_bruijn_graph.items():
    print(f"{node} -> {','.join(sorted(neighbors))}")

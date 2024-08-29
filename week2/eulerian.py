def find_eulerian_path(graph):
    # Calculate in-degree and out-degree
    in_degree = {}
    out_degree = {}
    for key in graph:
        if key not in out_degree:
            out_degree[key] = 0
        for value in graph[key]:
            if value not in in_degree:
                in_degree[value] = 0
            out_degree[key] += 1
            in_degree[value] += 1

    # Find the starting node
    start_node = None
    for node in graph:
        if node in out_degree and out_degree[node] > in_degree.get(node, 0):
            start_node = node
            break

    # Perform Hierholzer's algorithm
    path = []
    stack = [start_node]
    while stack:
        current_node = stack[-1]
        if current_node in graph and graph[current_node]:
            next_node = graph[current_node].pop(0)
            stack.append(next_node)
        else:
            path.append(stack.pop())

    # Reverse the path to get the Eulerian path
    path.reverse()
    return path

# sample = '''
# 0 -> 2
# 1 -> 3
# 2 -> 1
# 3 -> 0,4
# 6 -> 3,7
# 7 -> 8
# 8 -> 9
# 9 -> 6
# '''

with open('rosalind_ba3g.txt') as f:
    sample = f.read()

# Parse the input
graph = {}
for line in sample.strip().split('\n'):
    node, edges = line.split(' -> ')
    graph[node] = edges.split(',')
    for i in range(len(graph[node])):
        graph[node][i] = graph[node][i].strip()

# Test the function
eulerian_path = find_eulerian_path(graph)
# print(eulerian_path)
print('->'.join(eulerian_path))

# Save the output
with open('rosalind_ba3g.out', 'w') as f:
    f.write('->'.join(eulerian_path))

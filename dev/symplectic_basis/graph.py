import networkx as nx
import matplotlib.pyplot as plt


def debug_to_graph(lines):
    graph_segment = False
    graph_index = -1
    index = 0
    retval = [[]]

    for line in lines:
        list = line.split(' ')

        if graph_segment is False:
            graph_segment = list[0] == "Graph"
        elif list[0] == "Boundary":
            graph_segment = True
            graph_index += 1
            retval[index].append([])
        elif len(list) < 4 or list[4] != "Vertex":
            graph_segment = False
            index += 1
            graph_index = -1
            retval.append([])
        else:
            vertex = int(list[5])
            edges = list[12:]
            for edge in edges[:-1]:
                retval[index][graph_index].append((vertex, int(edge)))

    retval.pop()
    return retval


def file_to_graph(filename: str):
    with open(filename, "r") as file:
        graphs = debug_to_graph(file.readlines())
        g = [[nx.Graph() for _ in range(len(graphs[0]))] for _ in range(len(graphs))]

        for i, cusp in enumerate(graphs):
            for j, graph in enumerate(cusp):
                for v1, v2 in graph:
                    g[i][j].add_node(v1)
                    g[i][j].add_node(v2)
                    g[i][j].add_edge(v1, v2)

    return g


if __name__ == "__main__":
    graphs = file_to_graph("logs/link-116084.log")
    pos = nx.spring_layout(graphs[0][0])
    for i, graph in enumerate(graphs):
        plt.figure(i)
        nx.draw(graph[0], pos=nx.spring_layout(graph[0]), with_labels=True, font_weight="bold")

    plt.show()

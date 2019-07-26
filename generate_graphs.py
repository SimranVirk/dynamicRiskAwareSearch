#Generate large graphs
import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np
# import rags_test as rags 

class GraphGenerator:
	def __init__(self, mean_max, var_max, radius_connect):
		self.graphnum = 0
		self.mean_max = mean_max
		self.var_max = var_max
		self.radius_connect = radius_connect
		self.g = nx.DiGraph()
		return

	def gen_graph(self, num_verts, max_w, max_h):
		g = nx.DiGraph()
		g.add_nodes_from(range(num_verts))
		node_locs = {}
		attr = {}
		seen_locs = [[0,0], [max_w, max_h]]

		
		i = 0
		while i < num_verts:
			if i == 0:
				location = np.array([0,0])
			elif i == num_verts - 1:
				location = np.array([max_w, max_h])

			else:
				loc = [random.randint(0, max_w), random.randint(0, max_h)]
				if loc in seen_locs:
					continue

				location = np.array(loc)
				seen_locs.append(loc)

			node_locs[i] = {'loc':location}
			attr[i] = location
			i += 1

		nx.set_node_attributes(g, node_locs)


		edges = []
		for i in g.nodes:
			for j in g.nodes:
				dist = np.linalg.norm(g.nodes[i]['loc'] - g.nodes[j]['loc'])

				if dist < self.radius_connect and i != j:
					mean = dist + random.randint(0, self.mean_max)
					var = random.randint(0, self.var_max)
					g.add_edge(i, j, mean = mean, var = var)
		return g, attr



test = GraphGenerator(20, 5, 15)
graph, ns = test.gen_graph(20, 30, 30)
# nx.draw_networkx(graph, ns, edgelist = graph.edges)
# nx.draw_networkx(graph, with_labels = True)

# plt.axis('off')
# plt.show(block=True)  



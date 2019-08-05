#Generate large graphs
import networkx as nx
import random
# import matplotlib.pyplot as plt
import numpy as np
import random
# import rags_test as rags 
# from test import ragspath
# from rags_test import find_path

class GraphGenerator:
	def __init__(self, mean_max, var_max, radius_connect):
		self.graphnum = 0
		self.mean_max = mean_max
		self.var_max = var_max
		self.radius_connect = radius_connect
		self.g = nx.DiGraph()
		self.changes = {}
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


		self.g = g
		return g, attr


	def getEdgeChanges(self, num_changes):
		"""
			Create a set of edges to change and their final mean and variances
		"""
		edgesToChange = random.sample(self.g.edges, num_changes)
		for edge in edgesToChange:
			self.changes[edge] = []
			self.changes[edge].append(random.randint(0, self.mean_max))
			self.changes[edge].append(random.randint(0, self.var_max))


	def changeEdge(self, edge):
		"""
			Change edge u,v according to the stored changes
		"""
		meanvar = self.changes[edge]
		self.g.edges[edge[0], edge[1]]["mean"] = meanvar[0]
		self.g.edges[edge[0], edge[1]]["var"] = meanvar[1]



# test = GraphGenerator(20, 5, 15)
# graph, ns = test.gen_graph(20, 30, 30)
# nx.draw_networkx(graph, ns, edgelist = graph.edges)
# nx.draw_networkx(graph, with_labels = False)

# # plt.axis('off')
# plt.show(block=True)  
# print(find_path(graph, 0, 19, 0.6))
# print(rags(graph, 0, 19, 0.9))


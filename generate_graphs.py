#Generate large graphs
import networkx as nx
import random
import numpy as np
import random
import copy


class GraphGenerator:
	def __init__(self, mean_max, var_max):
		self.g = None
		self.changes = None
		self.changes1 = None

		self.mean_max = mean_max
		self.var_max = var_max
		return

	def gen_graphs(self, num_graphs, num_verts, mean_max, var_max, radius_connect, save = False):
		for i in range(num_graphs):
			g, ns = self.gen_graph(num_verts, mean_max, var_max, radius_connect)
			# print(g.edges(data = True))

			if save:
				self.save_graph(g, num_verts, i)

	def gen_graph(self, num_verts, mean_max, var_max, radius_connect):
		g = nx.DiGraph()
		g.add_nodes_from(range(num_verts))
		node_locs = {}
		attr = {}
		max_w = num_verts + 10
		max_h = num_verts + 10
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

				if dist < radius_connect and i != j:
					mean = int(dist)  + random.randint(0, mean_max)
					var = random.randint(1, var_max)
					g.add_edge(i, j, mean = mean, var = var)


		self.g = g
		return g, attr


	def getEdgeChanges(self, num_changes):
		"""
			Create a set of edges to change and their final mean and variances
		"""
		self.changes = {}
		self.changes1 = {}

		edgesToChange = random.sample(self.g.edges, num_changes)
		for edge in edgesToChange:
			self.changes[edge] = []
			self.changes[edge].append(random.randint(0, self.mean_max))
			self.changes[edge].append(random.randint(0, self.var_max))

		self.changes1 = copy.deepcopy(self.changes)


	def changedEdge(self, edge):
		"""
			Change edge u,v according to the stored changes
		"""
		x = self.changes.pop(edge, None)
		if x == None:
			print("ERROR removing change")


	def changedEdge1(self, edge):
		"""
			Change edge u,v according to the stored changes
		"""
		x = self.changes1.pop(edge, None)
		if x == None:
			print("ERROR removing change")

	def save_graph(self, graph, num_verts, graph_number):
		# print("gen ", num_verts)
		nx.write_edgelist(graph, "graph_files/%d_nodes/graph%d.edgelist" % (num_verts, graph_number), data = ['mean', 'var'])

	def load_graph(self, num_verts, graph_number):
		graph_read = nx.read_edgelist("graph_files/%d_nodes/graph%d.edgelist" % (num_verts, graph_number), create_using = nx.DiGraph(), nodetype = int, data = (('mean', int), ('var', int)))
		self.g = graph_read
		return graph_read

if __name__ == "__main__":
	gg = GraphGenerator(10, 10)
	for i in [60,70,80,90,100]:
		gg.gen_graphs(20, i, 10, 10, 20, save = True)
	# graph = gg.load_graph(20, 1)
	# print(graph.edges(data = True))


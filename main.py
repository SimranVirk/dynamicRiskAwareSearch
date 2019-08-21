import copy
import csv
import generate_graphs as gen
import networkx as nx
import drags

num_nodes = 60
total_num_edges = 0
num_changes = 20
iters = 5
start_goal_disconnected = 0
threshold = 0.75

print("Running comparison ", iters, " times for graph size ", num_nodes, " nodes and update frequency of one in ", num_changes)


# 1 -- dynamic rags, 2 -- replan from scratch- 
#total time, replan time, number of updates, cost
t1 = 0
rt1 = 0
u1 = 0
c1 = 0


t2 = 0
rt2 = 0
u2 = 0
c2 = 0


better = 0
bettercost  = 0

with open('results/update2/test%d_%d.csv' % (num_nodes, num_changes), 'a', newline = '') as csvfile:

	fieldnames = ['iter','num_nodes', 'num_edges', 'threshold', 'num_edges_changed', 
'd_totaltime', 'd_replantime', 'd_cost', 'd_meancost', 'd_varcost', 'd_numupdates',
'd_paths_expanded', 'd_heap_ops', 'd_paths_considered', 'd_closed_accesses',
'r_totaltime', 'r_replantime', 'r_cost', 'r_meancost', 'r_varcost', 'r_numupdates',
'r_paths_expanded', 'r_heap_ops', 'r_paths_considered', 'r_closed_accesses']
	writer = csv.DictWriter(csvfile, fieldnames = fieldnames)

	writer.writeheader()



	for i in range(iters):
		print("iteration ", i)
		gg = gen.GraphGenerator(10, 10)

		# Generate a new graph
		# graph, ns = gg.gen_graph(num_nodes, 10, 10, 20)

		# Load an existing graph
		graph = gg.load_graph(num_nodes, i)
		ns = []

		# # Uncomment to print graph
		# nx.draw_networkx_edges(graph, ns, edgelist = graph.edges, arrowsize = 2)
		# plt.axis('off')
		# plt.show()


		num_edges = len(graph.edges)
		total_num_edges += num_edges
		print(num_edges, " edges")

		graph_copy = copy.deepcopy(graph)
		d = drags.DRAGS(graph_copy, ns, 0, num_nodes - 1, threshold)
		gg.getEdgeChanges(int(num_edges / num_changes))
		print(num_edges/ num_changes, " edge changes")

		l2 = d.run_replan(gg, graph)

		# the start and goal were not connected in this graph so cancel the run
		if l2[0] == -1:
			start_goal_disconnected += 1
			continue

		l1 = d.run(gg)

		t1 += l1[0]
		rt1 += l1[1]
		u1 += l1[5]

		t2 += l2[0]
		rt2 += l2[1]
		u2 += l2[5]

		c1 += l1[2]
		c2 += l2[2]


		writer.writerow({'iter': i,'num_nodes': num_nodes, 'num_edges' : num_edges, 
				'threshold': threshold, 'num_edges_changed': num_changes, 
				'd_totaltime': l1[0], 'd_replantime': l1[1], 'd_cost': l1[2], 'd_meancost': l1[3], 'd_varcost': l1[4], 'd_numupdates': l1[5],
		'd_paths_expanded': l1[6], 'd_heap_ops': l1[7], 'd_paths_considered': l1[8], 'd_closed_accesses': l1[9],
		'r_totaltime': l2[0], 'r_replantime': l2[1], 'r_cost': l2[2], 'r_meancost': l2[3], 'r_varcost': l2[4], 'r_numupdates': l2[5],
		'r_paths_expanded': l2[6], 'r_heap_ops': l2[7], 'r_paths_considered': l2[8], 'r_closed_accesses': l2[9]})



		if l1[1] < l2[1]:
			better += 1

		if l1[3] < l2[3]:
			bettercost += 1

print("============")
print("stats - rt1: ", rt1, "t1", t1 , "u1: ",u1, "c1 ", c1)
print("rt2: ", rt2, "t2", t2, "u2: ", u2, "c2: ", c2)
print("better ", better)
print("better cost ", bettercost)
print("cancelled runs ", start_goal_disconnected)
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import csv

points = {20:[], 25:[], 30:[], 35:[], 40:[]}
rpoints = {20:[], 25:[], 30:[], 35:[], 40:[]}
p_points = {}

# num_changes = 20
sizes = [20, 25, 30, 35, 40]


for j in [5]:
	num_changes = j
	for i in sizes:
		num_nodes = i
		p_points[(num_nodes, num_changes)] = []


		with open('results/updated/updated_%d_%d.csv' % (num_nodes, num_changes), newline = '') as csvfile:
			r = csv.reader(csvfile)
			first_iter = True

			d_heap_ops = []
			d_path_expansions = []
			d_cost = []

			r_heap_ops = []
			r_path_expansions = []
			r_cost = []

			p_heap_ops = []
			p_path_expansions = []
			p_cost = []


			for row in r:
				if first_iter:
					first_iter = False
				else:
					d_heap = float(row[11])
					d_path_exp =  float(row[10])
					d_c = float(row[6])

					r_heap = float(row[20])
					r_path_exp =  float(row[19])
					r_c = float(row[15])

					d_heap_ops.append(d_heap)
					d_path_expansions.append(d_path_exp)
					d_cost.append(d_c)

					r_heap_ops.append(r_heap)
					r_path_expansions.append(r_path_exp)
					r_cost.append(r_c)

					p_path_expansions.append((r_path_exp - d_path_exp) * 100.0/d_path_exp)
					p_heap_ops.append((r_heap - d_heap) * 100.0/d_heap)
					p_cost.append((r_c - d_c) * 100.0/d_c)



			d_dict = d_heap_ops
			r_dict = r_heap_ops
			p_dict = p_heap_ops

			points[num_nodes].append(np.mean(np.array(d_dict).astype(float)))
			points[num_nodes].append(np.std(np.array(d_dict).astype(float)))

			rpoints[num_nodes].append(np.mean(np.array(r_dict).astype(float)))
			rpoints[num_nodes].append(np.std(np.array(r_dict).astype(float)))

			p_points[num_nodes, num_changes].append(np.mean(np.array(p_dict).astype(float)))
			p_points[num_nodes, num_changes].append(np.std(np.array(p_dict).astype(float)))


	x = np.array([x for x in points.keys()])
	y = np.array([v[0] for k,v in points.items()])
	e = np.array([v[1] for k,v in points.items()])

	plt.errorbar(x,y,e, label = "d-RAGS")

	x1 = np.array([x for x in rpoints.keys()])
	y1 = np.array([v[0] for k,v in rpoints.items()])
	e1 = np.array([v[1] for k,v in rpoints.items()])
	plt.errorbar(x1,y1,e1, color = 'r', label = "RAGS")
	plt.title("Number of Path Expansions vs Graph Size")
	plt.ylabel("# Path Expansions")
	plt.xlabel("Graph Size")
	plt.xticks(sizes)

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))


	plt.savefig("diagrams/errorbarplot/single_edgechangenum/heap_ops_legend%d" % (num_changes), bbox_inches = "tight")
	# plt.show()


# x = np.array([x for x in sizes])

# y5 = np.array([v[0] for k,v in p_points.items() if k[1] == 5])
# e5 = np.array([v[1] for k,v in p_points.items() if k[1] == 5])
# plt.errorbar(x,y5,e5, color = 'blue', label = '1 in 5 edges changed')

# y10 = np.array([v[0] for k,v in p_points.items() if k[1] == 10])
# e10 = np.array([v[1] for k,v in p_points.items() if k[1] == 10])
# plt.errorbar(x,y10,e10, color = 'red', label = '1 in 10 edges changed')

# y20 = np.array([v[0] for k,v in p_points.items() if k[1] == 20])
# e20 = np.array([v[1] for k,v in p_points.items() if k[1] == 20])
# plt.errorbar(x,y20,e20, color = 'orange', label = '1 in 20 edges changed')

# plt.title("Percentage Extra Heap Operations")
# plt.ylabel("Heap Operations")
# plt.xlabel("Graph Size")
# plt.xticks(sizes)

# # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# plt.savefig("diagrams/errorbarplot/percentage_difference/heap_ops", bbox_inches = "tight")
# # plt.show()



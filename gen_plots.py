import numpy as np
import matplotlib.pyplot as plt
import csv

points = {20:[], 25:[], 30:[], 35:[], 40:[]}
rpoints = {20:[], 25:[], 30:[], 35:[], 40:[]}

num_changes = 20
sizes = [20, 25, 30, 35, 40]

for i in sizes:
	num_nodes = i

	with open('results/updated/updated_%d_%d.csv' % (num_nodes, num_changes), newline = '') as csvfile:
		r = csv.reader(csvfile)
		first_iter = True

		d_heap_ops = []
		d_path_expansions = []
		d_cost = []

		r_heap_ops = []
		r_path_expansions = []
		r_cost = []


		for row in r:
			if first_iter:
				first_iter = False
			else:
				d_heap_ops.append(float(row[11]))
				d_path_expansions.append(float(row[10]))
				d_cost.append(float(row[6]))

				r_heap_ops.append(float(row[20]))
				r_path_expansions.append(float(row[19]))
				r_cost.append(float(row[15]))

		d_dict = d_heap_ops
		r_dict = r_heap_ops

		points[num_nodes].append(np.mean(np.array(d_dict).astype(float)))
		points[num_nodes].append(np.std(np.array(d_dict).astype(float)))

		rpoints[num_nodes].append(np.mean(np.array(r_dict).astype(float)))
		rpoints[num_nodes].append(np.std(np.array(r_dict).astype(float)))


x = np.array([x for x in points.keys()])
y = np.array([v[0] for k,v in points.items()])
e = np.array([v[1] for k,v in points.items()])

plt.errorbar(x,y,e)

x1 = np.array([x for x in rpoints.keys()])
y1 = np.array([v[0] for k,v in rpoints.items()])
e1 = np.array([v[1] for k,v in rpoints.items()])
plt.errorbar(x1,y1,e1,)
plt.show()







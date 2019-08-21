import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import csv

points = {60:[], 70:[], 80:[], 90:[]}
rpoints = {60:[], 70:[], 80:[], 90:[]}
p_points = {}

sizes = [60,70,80,90]


cols = []




def read_csvs(graph_sizes, update_frequencies):

	csv_dicts = {}

	for u in update_frequencies:
		num_changes = u
		for i in graph_sizes:
			csv_dict = {}

			num_nodes = i
			p_points[(num_nodes, u)] = []


			with open('results/update2/updated_%d_%d.csv' % (num_nodes, u), newline = '') as csvfile:
				r = csv.reader(csvfile)

				first_iter = True

				for row in r:
					if first_iter:
						first_iter = False
						for colname in row:
							csv_dict[colname] = []
							cols.append(colname)
					else:
						for i in range(len(row)):
							csv_dict[cols[i]].append(float(row[i]))

			csv_dicts[num_nodes, num_changes] = csv_dict

	return csv_dicts




def gen_percent_comparison(csv_dicts, graph_sizes, colname, axisname, update_freqs):

	#map graphsize, update freq to mean, stddev of column
	points  = {}

	x = np.array([x for x in sizes])

	for u in update_freqs:
		for g in graph_sizes:
			file = csv_dicts[g, u]
			diffs = np.array(np.array(file["r_" + colname]) - np.array(file["d_" + colname]))*100.0 / np.array(file["d_" + colname])

			points[g,u] = [np.mean(diffs), np.std(diffs)]



	for u in update_freqs:
		x = []
		y = []
		e = []
		for k,v in points.items():
			if k[1] == u:
				x.append(k[0])
				y.append(v[0])
				e.append(v[1])

		plt.errorbar(x,y,e, label = '1 in %d edges changed'%(u))

	plt.title("Percentage Extra %s"%(axisname))
	plt.ylabel("Percent")
	plt.xlabel("Graph Size")
	plt.xticks(sizes)

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	plt.savefig("diagrams/update2/percent/%s" % (axisname), bbox_inches = "tight")
	plt.show()


def gen_twoline_plot(csv_dicts, graph_sizes, colname, axisname, update_freq):
	dcol = "d_" + colname
	rcol = "r_" + colname

	#map graph size to mean and stddev of column
	drags = {}
	rags = {}

	for k,v in csv_dicts.items():
		if k[1] == update_freq:
			m = np.mean(np.array(v[dcol]))
			s = np.std(np.array(v[dcol]))

			drags[k[0]] = [m, s]

			mr = np.mean(np.array(v[rcol]))
			sr = np.std(np.array(v[rcol]))

			rags[k[0]] = [mr, sr]

	x = np.array([x for x in graph_sizes])
	y = np.array([v[0] for k,v in drags.items()])
	e = np.array([v[1] for k,v in drags.items()])
	plt.errorbar(x,y,e, label = "d-RAGS")

	y1 = np.array([v[0] for k,v in rags.items()])
	e1 = np.array([v[1] for k,v in rags.items()])
	plt.errorbar(x,y1,e1, color = 'r', label = "RAGS")

	plt.title("%s vs Graph Size" % (axisname))
	plt.ylabel(axisname)
	plt.xlabel("Graph Size")
	# # plt.ylim(0)
	plt.xticks(sizes)

	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))


	plt.savefig("diagrams/update2/single/%s_%d" % (axisname, update_freq), bbox_inches = "tight")
	plt.show()



if __name__ == "__main__":
	d = read_csvs(sizes, [5,20])

	# gen_twoline_plot(d, sizes, "cost", "Cost", 5)
	gen_twoline_plot(d, sizes, "replantime", "Time", 5)
	# gen_percent_comparison(d, sizes, "paths_expanded", "Heap Operations", [5,20])





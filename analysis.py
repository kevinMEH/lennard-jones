import numpy as np
import matplotlib.pyplot as plot

def load_file(fileName):
    file = open(fileName, "r")
    contents = file.read().strip()
    contents = contents.split("\n\n\n")[2]
    results = contents.split("\n")
    array = []
    for result in results:
        values = []
        for csv in result.split(","):
            values.append(float(csv.strip()))
        array.append(values)
    array = np.asarray(array)
    return array

all_all_values = [
    load_file("./results/analysis/decreasing_std/5k_record/particle_analysis_1_0.999000_0.999913.txt"),
    load_file("./results/analysis/decreasing_std/5k_record/particle_analysis_2_0.999000_0.999913.txt"),
    load_file("./results/analysis/decreasing_std/5k_record/particle_analysis_3_0.999000_0.999913.txt"),
]

record_every = 5
batch_size = 1000

colors = [ "tab:blue", "tab:green", "tab:orange", "tab:red", "tab:pink" ]

record_thresholds = np.linspace(record_every * batch_size, record_every * batch_size * all_all_values[0][0].size, all_all_values[0][0].size)

figure, (axes1, axes2) = plot.subplots(2, 1, figsize=(13, 9))
plot.subplots_adjust(wspace=0.2, hspace=0.3)
plot.subplots_adjust(left=0.075, right=0.95)
plot.subplots_adjust(top=0.925, bottom=0.075)

for all_values, color in zip(all_all_values, colors):
    best_potentials = np.asarray(all_values[0])
    acceptance_rate = np.asarray(all_values[1])

    axes1.plot(record_thresholds, best_potentials, c=color)
    for x, y in list(zip(record_thresholds, best_potentials))[:3]:
        axes1.annotate("{:.0f}".format(y), (x + 1500, y), va="center", c=color)

    axes2.plot(record_thresholds, acceptance_rate, c=color)

axes1.set_xticks(record_thresholds[::2])
axes1.set_title("Batch potentials over time")
axes1.set_xlabel("Steps")
axes1.set_ylabel("Average batch potential")

axes2.set_xticks(record_thresholds[::2])
axes2.set_title("Batch potentials over time")
axes2.set_xlabel("Steps")
axes2.set_ylabel("Overall acceptance rate")

plot.show()
import numpy as np
import matplotlib.pyplot as plot

fileName = "particle_1_64.txt"
file = open("./results/" + fileName, "r")
contents = file.read().strip()

results = contents.split("\n\n")
results = [ result.split("\n") for result in results ]
all_values = []
for result in results:
    values = []
    for line in result:
        value = float(line.split(": ")[1].strip())
        values.append(value)
    all_values.append(values)

temperatures = []
standard_deviations = []
best_potentials = []
converging_batches = []
converging_acceptance_rates = []

for [_temperature, _standard_deviation, _best_potential, _, _converging_batches, _, _converging_acceptance_rate, _] in all_values:
    temperatures.append(_temperature)
    standard_deviations.append(_standard_deviation)
    best_potentials.append(_best_potential)
    converging_batches.append(_converging_batches)
    converging_acceptance_rates.append(_converging_acceptance_rate)

converging_batches = np.asarray(converging_batches)
real_converging_batches = converging_batches.copy()
min_converging_batches = converging_batches[np.argmin(converging_batches)]
max_converging_batches = converging_batches[np.argmax(converging_batches)]
converging_batches = (converging_batches - min_converging_batches) / (max_converging_batches - min_converging_batches)
# Min size = 20, max size = 200
converging_batches = 200 / (1 + converging_batches * 10)

figure, axes = plot.subplots(figsize=(16, 10));
contour = axes.tricontourf(temperatures, standard_deviations, best_potentials, 100)
# The bigger the circle, the less the batch size
scatter = axes.scatter(temperatures, standard_deviations, s=converging_batches, fc="w", ec="k")
background_color = dict(facecolor="k", alpha=0.25)
for temp, std, batch_count, accept_rate in zip(temperatures, standard_deviations, real_converging_batches, converging_acceptance_rates):
    axes.annotate("({:.0f}, {:.3f})".format(batch_count, accept_rate), (temp + 1, std + 0.0015), rotation=40, c="w", bbox=background_color)
axes.set_xlabel("Temperatures")
axes.set_ylabel("Standard Deviations")
axes.set_xticks([0, 1, 10, 20, 35, 50, 75, 100, 105])
axes.set_yticks([0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 0.9])
figure.colorbar(contour, ax=axes)

plot.show()
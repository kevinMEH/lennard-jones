import numpy as np
import matplotlib.pyplot as plot

fileName = "particle_increasing_temperature.txt"
file = open("./results/" + fileName, "r")
contents = file.read().strip()

contents = contents.split("\n\n\n")[1]

results = contents.split("\n\n")
results = [ result.split("\n") for result in results ]
all_values = []
for result in results:
    values = []
    for line in result:
        values.append(float(line.split(": ")[1].strip()))
    all_values.append(values)

move_particles = np.asarray([ values[0] for values in all_values ])
best_potentials_mean = np.asarray([ values[3] for values in all_values ])
converging_batches = np.asarray([ values[5] for values in all_values ])

batch_size = 1000
limit = 15

move_particles = move_particles[:limit]
best_potentials_mean = best_potentials_mean[:limit]
converging_batches = converging_batches[:limit]

xticks = np.arange(1, 16)

figure, ((axes1, axes2), (axes3, axes4)) = plot.subplots(2, 2, figsize=(16, 9))
plot.subplots_adjust(wspace=0.2, hspace=0.3)
plot.subplots_adjust(left=0.075, right=0.95)
plot.subplots_adjust(top=0.925, bottom=0.075)
axes1.plot(move_particles, best_potentials_mean)
for x, y in zip(move_particles[[0, 1, 2]], best_potentials_mean[[0, 1, 2]]):
    axes1.scatter(x, y, c="tab:blue")
    axes1.annotate("{:.0f}".format(y), (x + 0.6, y), va="center")
axes1.set_title("Best potential vs Particles moved")
axes1.set_xticks(xticks)
axes1.set_xlabel("Particle moved per step")
axes1.set_ylabel("Mean best potential")

axes3.plot(move_particles, converging_batches)
for x, y in zip(move_particles[[0, 1, 2]], converging_batches[[0, 1, 2]]):
    axes3.scatter(x, y, c="tab:blue")
    axes3.annotate("{:.0f}".format(y), (x + 0.6, y), va="center")
axes3.set_title("Total steps vs Particles moved")
axes3.set_xticks(xticks)
axes3.set_xlabel("Particle moved per step")
axes3.set_ylabel("Total steps (in 1000s)")

axes4.plot(move_particles, move_particles * converging_batches)
for x, y in zip(move_particles[[0, 1, 2]], (move_particles * converging_batches)[[0, 1, 2]]):
    axes4.scatter(x, y, c="tab:blue")
    axes4.annotate("{:.0f}".format(y), (x + 0.6, y), va="center")
axes4.set_title("Total particle calculations vs Particles moved")
axes4.set_xticks(xticks)
axes4.set_xlabel("Particle moved per step")
axes4.set_ylabel("Total particle calculations (in 1000s)")

plot.show()

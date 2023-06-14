import numpy as np
import matplotlib.pyplot as plot

fileName = "particle_constant.txt"
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

move_particles_constant = np.asarray([ values[0] for values in all_values ])
best_potentials_mean_constant = np.asarray([ values[2] for values in all_values ])
converging_batches_constant = np.asarray([ values[4] for values in all_values ])

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

move_particles_increasing = np.asarray([ values[0] for values in all_values ])
best_potentials_mean_increasing = np.asarray([ values[3] for values in all_values ])
converging_batches_increasing = np.asarray([ values[5] for values in all_values ])

batch_size = 1000
limit = 15

move_particles_constant = move_particles_constant[:limit]
best_potentials_mean_constant = best_potentials_mean_constant[:limit]
converging_batches_constant = converging_batches_constant[:limit]

move_particles_increasing = move_particles_increasing[:limit]
best_potentials_mean_increasing = best_potentials_mean_increasing[:limit]
converging_batches_increasing = converging_batches_increasing[:limit]

xticks = np.arange(1, 16)

figure, ((axes1, axes2), (axes3, axes4)) = plot.subplots(2, 2, figsize=(16, 9))
plot.subplots_adjust(wspace=0.2, hspace=0.3)
plot.subplots_adjust(left=0.075, right=0.95)
plot.subplots_adjust(top=0.925, bottom=0.075)

axes1.plot(move_particles_constant, best_potentials_mean_constant, c="tab:blue")
axes1.plot(move_particles_increasing, best_potentials_mean_increasing, c="tab:orange", zorder=-10)
for x, y in zip(move_particles_constant[[0, 1, 2]], best_potentials_mean_constant[[0, 1, 2]]):
    axes1.scatter(x, y, c="tab:blue")
    axes1.annotate("{:.0f}".format(y), (x + 0.6, y), c="tab:blue", va="center")
for x, y in zip(move_particles_increasing[[0, 1, 2]], best_potentials_mean_increasing[[0, 1, 2]]):
    axes1.scatter(x, y, c="tab:orange", zorder=-10)
    axes1.annotate("{:.0f}".format(y), (x - 0.4, y), c="tab:orange", va="center", ha="right", zorder=-10)
axes1.set_title("Best potential vs Particles moved")
axes1.set_xticks(xticks)
axes1.set_xlabel("Particle moved per step")
axes1.set_ylabel("Mean best potential")

axes2.annotate("Blue = constant temperature (25)", (0.5, 0.6), c="tab:blue", ha="center", fontsize=14)
axes2.annotate("Orange = increasing temperature (25 * particleCount)", (0.5, 0.3), c="tab:orange", ha="center", fontsize=14)

axes3.plot(move_particles_constant, converging_batches_constant, c="tab:blue")
axes3.plot(move_particles_increasing, converging_batches_increasing, c="tab:orange", zorder=-10)
for x, y in zip(move_particles_constant[[0, 1, 2]], converging_batches_constant[[0, 1, 2]]):
    axes3.scatter(x, y, c="tab:blue")
    axes3.annotate("{:.0f}".format(y), (x + 0.6, y), c="tab:blue", va="center")
for x, y in zip(move_particles_increasing[[0, 1, 2]], converging_batches_increasing[[0, 1, 2]]):
    axes3.scatter(x, y, c="tab:orange", zorder=-10)
    axes3.annotate("{:.0f}".format(y), (x - 0.4, y), c="tab:orange", va="center", ha="right", zorder=-10)
axes3.set_title("Total steps vs Particles moved")
axes3.set_xticks(xticks)
axes3.set_xlabel("Particle moved per step")
axes3.set_ylabel("Total steps (in 1000s)")

axes4.plot(move_particles_constant, move_particles_constant * converging_batches_constant, c="tab:blue")
axes4.plot(move_particles_increasing, move_particles_increasing * converging_batches_increasing, c="tab:orange", zorder=-10)
for x, y in zip(move_particles_constant[[0, 1, 2]], (move_particles_constant * converging_batches_constant)[[0, 1, 2]]):
    axes4.scatter(x, y, c="tab:blue")
    axes4.annotate("{:.0f}".format(y), (x + 0.6, y), c="tab:blue", va="center")
for x, y in zip(move_particles_increasing[[0, 1, 2]], (move_particles_increasing * converging_batches_increasing)[[0, 1, 2]]):
    axes4.scatter(x, y, c="tab:orange", zorder=-10)
    axes4.annotate("{:.0f}".format(y), (x - 0.4, y), c="tab:orange", va="center", ha="right", zorder=-10)
axes4.set_title("Total particle calculations vs Particles moved")
axes4.set_xticks(xticks)
axes4.set_xlabel("Particle moved per step")
axes4.set_ylabel("Total particle calculations (in 1000s)")

plot.show()

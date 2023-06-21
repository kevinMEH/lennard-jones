import numpy as np
import matplotlib.pyplot as plot

def load_values(fileName):
    file = open(fileName, "r")
    contents = file.read().strip()
    contents = contents.split("\n\n\n")[1]
    lines = contents.split("\n")
    line_values = []
    for line in lines:
        split = line.split(": ")
        if(len(split) != 2): continue
        line_values.append(float(split[1]))
    return [ line_values[3], line_values[5], line_values[6], line_values[8] ]

values = [
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999653.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999700.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999768.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999800.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999827.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999850.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999866.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999884.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999913.txt"),
    load_values("./results/analysis/decreasing_std/5k_record/0.999000/particle_analysis_1_0.999000_0.999931.txt"),
]

temperature_factors = [ value[0] for value in values ]
std_factors = [ value[1] for value in values ]
best_potentials = [ value[2] for value in values ]
converging_steps = [ value[3] for value in values ]

record_every = 5
batch_size = 1000

figure, (axes1, axes2) = plot.subplots(2, 1, figsize=(13, 9))
plot.subplots_adjust(wspace=0.2, hspace=0.3)
plot.subplots_adjust(left=0.075, right=0.95)
plot.subplots_adjust(top=0.925, bottom=0.075)

axes1.plot(std_factors, best_potentials)
for factor, potential in zip(std_factors, best_potentials):
    axes1.annotate("{:0f}".format(potential), (factor, potential - 2), ha="center", va="top")
    if(abs(factor - 0.999850) < 0.00001):
        axes1.scatter(factor, potential, c="tab:orange")
axes2.plot(std_factors, converging_steps)
for factor, steps in zip(std_factors, converging_steps):
    axes2.annotate("{:0f}".format(steps), (factor, steps - 10000), ha="center", va="top")
    if(abs(factor - 0.999850) < 0.00001):
        axes2.scatter(factor, steps, c="tab:orange")

axes1.set_title("Best potential vs STD")
axes1.set_xticks(std_factors)
axes1.set_xlabel("STD")
axes1.set_ylabel("Best potential")
axes1.set_ybound(-315)

axes2.set_title("Converging steps vs STD")
axes2.set_xticks(std_factors)
axes2.set_xlabel("STD")
axes2.set_ylabel("Converging steps")
axes2.set_ybound(0)

plot.show()
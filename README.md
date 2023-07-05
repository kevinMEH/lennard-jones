# Lennard Jones

Extremely fast program for optimization of a Lennard Jones system using the
Metropolis Hastings and Simulated Annealing algorithms.

## Considerations

The algorithm is very fast. Each particle is associated with a potential array,
which contains the potential between that particle and every other particle in
the system. Each time the particle is moved, that potential array is updated,
and that particle's entry in every other particle's potential array is updated.

Compared to the naive algorithm of computing the potentials for every single
particle and summing it up at the end of each step, this algorithm is magnitudes
faster. (A complexity of O(1) vs O(n^2))

However, this also means that more memory is used for the simulation. The memory
consumption is approximately `sizeof(Particle) * n + sizeof(double) * n^2`,
where n is the amount of particles in the system. (As opposed to just
`sizeof(Particle) * n` for the naive algorithm.) So be mindful of memory usage
when running large scale simulations.

The algorithm is multi-thread friendly.

Feel free to modify the code for your own usages. A lot of the code is unique
for my project, like for example the moving multiple particles code, the warm
start code, the recording potential values and acceptance rates code, etc.

## Original Research Question

This program was built for a research project. The question was: "Will moving
multiple particles during each step at the beginning of a simulation result in
an improvement in the number of steps needed for convergence?"

The answer is no, not really.

The potential when moving multiple particles at once yields only extremely
marginal benefits for the first few steps of the simulation.

---

Let `f(t)` be the potential vs time function.

Let `f'(t)` be the change in potential per step. It is the derivative of the
potential vs time function, `f(t)`.

Let `f(t)_x` and `f'(t)_x` be `f(t)` and `f'(t)` respectively for simulations
moving x particles per step.

For `f'(t)`, it should ideally be a negative value, as we want our potential to
decreaese. So, the lower the better.

For all simulations, `f(t)_x` starts off at the same potential value, which is
the potential of a bunch of random particles distributed in a bounding box.

---

Running the experiments and recording the average potentials after every 1000
steps, we find that moving 2 particles per step vs 1 particle per step results
in a lower potential for only the first ~2000 steps of the simulation.

After that, moving multiple particles results in worse performance.

This means that initially, `f'(t)_2` starts off at a lower value than `f'(t)_1`,
which results in `f(t)_2` being lower than `f(t)_1`. However, after a while,
`f'(t)_2` starts becoming greater than `f'(t)_1`, making moving 2 particles a
worse option. The goal is to change from moving 2 particles to 1 particle right
when `f'(t)_2` == `f'(t)_1`, so that we maximize the decrease in our potential
for every step.

From our results, `f(2000)_1` ~== `f(2000)_2`.

From this, let's assume that `f'(1000)_1` ~== `f'(1000)_2`, so it takes 1000
steps for the derivates for 1 particle and 2 particles to match. So we should
switch from moving 2 particles to 1 after 1000 steps.

To find the amount of steps we have saved, we need to find t such that
`f(t)_1` == `f(1000)_2`.

Let's assume that t = 1500, so `f(1500)_1 ~== f(1000)_2`. This means that we
have saved 500 steps from our simulation.

It takes about 25000 steps when moving 1 particle to get to a decent potential
of -295.

**1 - (24500 / 25000) = 2%**

So we've decreased the runtime by a whopping 2%!

Yeah, not really anything exciting to write a paper about. Especially
considering the fact that I implemented Lennard Jones in such a way that
particle potentials are recalculated in a very efficient manner every time a
particle is moved instead of recalculating all particles at the end of every
step.

This means that moving 1 particle for 50 steps and moving 50 particles for 1
step takes approximately the same time, so there is virtually no reason to even
consider moving multiple particles per step, except in mathematical theory.

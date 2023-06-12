#ifndef LENNARD_JONES
#define LENNARD_JONES

extern int TOTAL_PARTICLES;

extern int BATCH_SIZE;
extern int ALLOWED_STRIKES;

extern double BOUNDING_BOX_SIZE;
extern double HALF_BOUNDING_BOX_SIZE;

typedef struct Particle {
    double x;
    double y;
    double z;
    double* potentials; // Array of potentials to other particles
} Particle;

typedef struct SimulationResults {
    double bestPotential;
    int convergingBatches;
    double convergingAcceptanceRate;
} SimulationResults;

SimulationResults simulateAnnealing(int moveParticles, double temperature, double coolingFactor, double standardDeviation);

#endif
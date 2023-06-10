#ifndef LENNARD_JONES
#define LENNARD_JONES

extern int TOTAL_PARTICLES;
extern int MOVE_PARTICLES;

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

SimulationResults simulateAnnealing(double temperature, double cooling_factor, double standard_deviation);

#endif
#ifndef LENNARD_JONES
#define LENNARD_JONES

#define RECORD 0

#include "random.h"

extern int TOTAL_PARTICLES;

extern double BOUNDING_BOX_SIZE;
extern double HALF_BOUNDING_BOX_SIZE;

#if RECORD == 1
extern int recordSize;
extern int recordEvery;
#endif

typedef struct Particle {
    double x;
    double y;
    double z;
    double* potentials; // Array of potentials to other particles
} Particle;

typedef struct SimulationResults {
    double bestPotential;
    int convergingSteps;
    int convergingParticleComputations;
    double convergingAcceptanceRate;
} SimulationResults;

typedef struct LennardJonesOptions {
    #if RECORD == 1
    double* bestPotentialRecord;
    double* acceptanceRateRecord;
    #endif

    int warmStart;
    int warmStartTemperature;
    int warmStartCount;
    
    int BATCH_SIZE;
    int ALLOWED_STRIKES;

    double minimumImprovement;
} LennardJonesOptions;

SimulationResults simulateAnnealing(
    XorshiftState* state,
    int moveParticles, double particleFactor,
    double temperature, double temperatureFactor,
    double standardDeviation, double standardDeviationFactor,
    LennardJonesOptions options
);

#endif
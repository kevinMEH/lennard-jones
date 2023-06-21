#ifndef LENNARD_JONES
#define LENNARD_JONES

#define RECORD 0

#include "random.h"

extern int TOTAL_PARTICLES;

extern int BATCH_SIZE;
extern int ALLOWED_STRIKES;

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

SimulationResults simulateAnnealing(XorshiftState* state,
    int moveParticles, double particleFactor,
    double temperature, double temperatureFactor,
    double standardDeviation, double standardDeviationFactor
    #if RECORD == 1
    , double* bestPotentialRecord, double* acceptanceRateRecord
    #endif
);

#endif
#include <stdio.h>
#include <time.h>
#include "src/random.h"
#include "src/lennard_jones.h"

int TOTAL_PARTICLES = 64;

double BOUNDING_BOX_SIZE = 4.0;
double HALF_BOUNDING_BOX_SIZE = 4.0/2;

int main() {
    int moveParticles = 2;
    double temperatureFactor = 0.999;
    double particleFactor = 1;
    double temperature = 25;
    double standard_deviation = 0.1;
    double standardDeviationFactor = 1.0;
    
    XorshiftState state;
    state.x[0] = 11992933392989292238ULL;
    state.x[1] = 995759136711242123ULL + time(0) * time(0) * time(0);
    
    LennardJonesOptions options;
    options.warmStart = 0;
    options.warmStartTemperature = 0;
    options.warmStartCount = 0;
    options.BATCH_SIZE = 512;
    options.ALLOWED_STRIKES = 3;
    options.minimumImprovement = 0;
    
    SimulationResults result = simulateAnnealing(&state, moveParticles, particleFactor, temperature, temperatureFactor, standard_deviation, standardDeviationFactor, options);
    printf("%lf\n", result.bestPotential);
    printf("%d\n", result.convergingSteps);
    printf("%lf\n", result.convergingAcceptanceRate);
    printf("%d\n", result.convergingParticleComputations);
}
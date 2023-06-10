#include <stdio.h>
#include <time.h>
#include "src/random.h"
#include "src/lennard_jones.h"

int TOTAL_PARTICLES = 10;
int MOVE_PARTICLES = 2;

int BATCH_SIZE = 1000;
int ALLOWED_STRIKES = 3;

double BOUNDING_BOX_SIZE = 4.0;
double HALF_BOUNDING_BOX_SIZE = 4.0/2;

/**
 * Geometric cooling: This number will be multiplied with the temperature after
 * every step
*/
double cooling_factor = 0.999;

/**
 * Temperature means, if the increase in potential is equal to ln(2) * temperature,
 * then we have a 1/2 chance of accepting this new state. If the increase in
 * Potential is equal to 2 * ln(2) * temperature, then we have a 1/4 chance of
 * accepting this new state. Etc. ln(2) = 0.693
*/
double temperature = 25;
double standard_deviation = 0.1;

int main() {
    if(TOTAL_PARTICLES > 64) {
        printf("Bitmap cannot represent more than 64 particles. Use an array instead.");
        return 1;
    }
    if(MOVE_PARTICLES > TOTAL_PARTICLES/2 && MOVE_PARTICLES != TOTAL_PARTICLES) {
        printf("Keep the amount of moving particles below 50%% of the total particles.");
        printf("The slowdown is approximately 1/(1 - move/total).");
    }
    
    setXorshiftState((uint64_t) 11992933392989292238ULL, (uint64_t) 995759136711242123ULL + time(0) * time(0) * time(0));
    SimulationResults result = simulateAnnealing(temperature, cooling_factor, standard_deviation);
    printf("%lf\n", result.bestPotential);
    printf("%d\n", result.convergingBatches);
    printf("%lf\n", result.convergingAcceptanceRate);
}
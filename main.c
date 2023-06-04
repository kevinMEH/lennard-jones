#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include "src/random.h"

#define TOTAL_PARTICLES 64
#define MOVE_PARTICLES 1

#define BOUNDING_BOX_SIZE 4.0
#define HALF_BOUNDING_BOX_SIZE BOUNDING_BOX_SIZE/2

/**
 * We are using 2^-(Delta / Temperature) instead of e^-(Delta / Temperature)
 * as the probability function.
 * So temperature means, if the increase in potential is equal to temperature,
 * then we have a 1/2 chance of accepting this new state. If the increase in
 * Potential is double the temperature, then we have a 1/4 chance of accepting
 * this new state. Etc.
*/
#define TEMPERATURE 1
#define STANDARD_DEVIATION 0.25
#define BATCH_SIZE 100000
#define ALLOWED_STRIKES 1

typedef struct Particle {
    double x;
    double y;
    double z;
    double* potentials; // Array of potentials to other particles
} Particle;

// Sigma = 1, epsilon = 1
double lennardJonesPotential(Particle* __restrict first, Particle* __restrict second) {
    double distanceSquared = pow(first->x - second->x, 2) + pow(first->y - second->y, 2) + pow(first->z - second->z, 2);
    return (4 / pow(distanceSquared, 6)) - (4 / pow(distanceSquared, 3));
}

// Recalculate potentials for the particle, changing the corresponding
// potentials of other particles.
void recalculatePotential(Particle* particles, int selfIndex) {
    Particle* self = &particles[selfIndex];
    for(int otherIndex = 0; otherIndex < TOTAL_PARTICLES; otherIndex++) {
        Particle* other = &particles[otherIndex];
        double potential = lennardJonesPotential(self, other);
        self->potentials[otherIndex] = potential;
        other->potentials[selfIndex] = potential;
    }
    self->potentials[selfIndex] = 0;
}

void moveParticle(Particle* particle) {
    particle->x = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->x + xorshiftNormal(0, STANDARD_DEVIATION)));
    particle->y = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->y + xorshiftNormal(0, STANDARD_DEVIATION)));
    particle->z = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->z + xorshiftNormal(0, STANDARD_DEVIATION)));
}

double totalPotential(Particle* particles) {
    double total = 0;
    for(int i = 0; i < TOTAL_PARTICLES; i++) {
        Particle* current = &particles[i];
        for(int j = 0; j < TOTAL_PARTICLES; j++) {
            total += current->potentials[j];
        }
    }
    return total;
}

int shouldMove(double changePotential, double temperature) {
    return xorshiftDouble() < (1 / exp2(changePotential / temperature));
}

// Takes an array of changed particles indices, copies over all the new
// coordinates and all the potentials. The potentials array pointer of the
// destination will be untouched. (It will be changed but resetted.)
void copyOverAllParticles(Particle* __restrict source, Particle* __restrict destination) {
    double* destinationBlockPotentials = destination[0].potentials;
    memcpy(destination, source, sizeof(Particle) * TOTAL_PARTICLES);
    memcpy(destinationBlockPotentials, source[0].potentials, sizeof(double) * TOTAL_PARTICLES * TOTAL_PARTICLES);
    for(int i = 0; i < TOTAL_PARTICLES; i++) {
        destination[i].potentials = &destinationBlockPotentials[i * TOTAL_PARTICLES];
    }
}

int main() {
    setXorshiftState((uint64_t) 11992933392989292238ULL, (uint64_t) 995759136711242123ULL + time(0) * time(0) * time(0));

    Particle particles[TOTAL_PARTICLES];
    double* potentialsBlock = (double*) calloc(TOTAL_PARTICLES * TOTAL_PARTICLES, sizeof(double));
    
    // Initialize particles randomly in initial box
    for(int i = 0; i < TOTAL_PARTICLES; i++) {
        particles[i] = (Particle) {
            (xorshiftDouble() * BOUNDING_BOX_SIZE - HALF_BOUNDING_BOX_SIZE),
            (xorshiftDouble() * BOUNDING_BOX_SIZE - HALF_BOUNDING_BOX_SIZE),
            (xorshiftDouble() * BOUNDING_BOX_SIZE - HALF_BOUNDING_BOX_SIZE),
            &potentialsBlock[i * TOTAL_PARTICLES],
        };
    }
    // Initial potential calculation
    for(int i = 0; i < TOTAL_PARTICLES; i++) {
        recalculatePotential(particles, i);
    }
    
    // Modified particles: These are the particles we will modify, and then copy
    // to the source if the SA process succeeds
    Particle modifiedParticles[TOTAL_PARTICLES];
    double* modifiedPotentialsBlock = (double*) calloc(TOTAL_PARTICLES * TOTAL_PARTICLES, sizeof(double));
    modifiedParticles[0].potentials = modifiedPotentialsBlock;
    copyOverAllParticles(particles, modifiedParticles);
    
    // Best particles
    Particle bestParticles[TOTAL_PARTICLES];
    double* bestPotentialsBlock = (double*) calloc(TOTAL_PARTICLES * TOTAL_PARTICLES, sizeof(double));
    bestParticles[0].potentials = bestPotentialsBlock;
    copyOverAllParticles(particles, bestParticles);
    

    // Start of Simulated Annealing
    /**
     * How convengence is determined: A variable n, called the BATCH_SIZE, is
     * defined. Every n iterations, the average potential is calculated. As long
     * as that potential is less than the best batch's potential, no strike is
     * added. If the potential is greater, a strike is added. If on the next
     * iteration, the potential becomes lower again, the strikes are removed. If
     * the amount of strikes reaches the amount of ALLOWED_STRIKES, then the
     * program ends.
    */
    
    double lastPotential = totalPotential(particles);
    double bestPotential = lastPotential;
    double bestBatchPotential = bestPotential * 1.1;
    int strikes = 0;
    
    while(1) {
        double currentBatchTotalPotential = 0.0;
        for(int i = 0; i < BATCH_SIZE; i++) {
            #if TOTAL_PARTICLES > 64
            printf("Bitmap cannot represent more than 64 particles. Use an array instead.");
            exit(1);
            #endif
            unsigned long long changed = 0;
            for(int i = 0; i < MOVE_PARTICLES; i++) {
                #if TOTAL_PARTICLES != MOVE_PARTICLES
                int indexToMove = xorshiftDouble() * TOTAL_PARTICLES;
                while(changed & (1ULL << indexToMove)) {
                    indexToMove = xorshiftDouble() * TOTAL_PARTICLES;
                }
                changed = changed | (1ULL << indexToMove);
                #else
                int indexToMove = i;
                #endif
                moveParticle(&modifiedParticles[indexToMove]);
                recalculatePotential(modifiedParticles, indexToMove);
            }
            double modifiedTotalPotential = totalPotential(modifiedParticles);
            if(shouldMove(modifiedTotalPotential - lastPotential, TEMPERATURE)) {
                copyOverAllParticles(modifiedParticles, particles);
                lastPotential = modifiedTotalPotential;
                if(lastPotential < bestPotential) {
                    copyOverAllParticles(modifiedParticles, bestParticles);
                    bestPotential = lastPotential;
                }
            } else {
                copyOverAllParticles(particles, modifiedParticles);
            }
            currentBatchTotalPotential += lastPotential;
        }
        double currentBatchPotential = currentBatchTotalPotential / BATCH_SIZE;
        printf("Current batch potential: %lf\n", currentBatchPotential);
        if(currentBatchPotential >= bestBatchPotential) {
            strikes++;
            if(strikes == ALLOWED_STRIKES) break;
        } else {
            strikes = 0;
            bestBatchPotential = currentBatchPotential;
        }
    }
    
    printf("Best potential: %lf\n", bestPotential);
    // printf("Best potential particles:\n");
    // for(int i = 0; i < TOTAL_PARTICLES; i++) {
    //     Particle bestParticle = bestParticles[i];
    //     printf("%.9lf,%.9lf,%.9lf\n", bestParticle.x, bestParticle.y, bestParticle.z);
    // }
    
    free(potentialsBlock);
    free(modifiedPotentialsBlock);
    free(bestPotentialsBlock);
}
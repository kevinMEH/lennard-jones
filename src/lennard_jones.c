#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "random.h"
#include "lennard_jones.h"

#define LARGE_MOVE_SIZE 0
#define PARTICLE_GEOMETRIC_COOLING 1
#define MULTI_THREADED 1

double changeMoveParticles(double previousMoveParticlesDouble, double coolingFactor) {
    return previousMoveParticlesDouble * coolingFactor;
}

double changeTemperature(double previousTemperature, double coolingFactor) {
    return previousTemperature * coolingFactor;
}

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

void moveParticle(Particle* particle, double standardDeviation) {
    particle->x = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->x + xorshiftNormal(0, standardDeviation)));
    particle->y = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->y + xorshiftNormal(0, standardDeviation)));
    particle->z = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->z + xorshiftNormal(0, standardDeviation)));
}

void moveParticleState(XorshiftState* state, Particle* particle, double standardDeviation) {
    particle->x = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->x + xorshiftNormalState(state, 0, standardDeviation)));
    particle->y = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->y + xorshiftNormalState(state, 0, standardDeviation)));
    particle->z = fmin(HALF_BOUNDING_BOX_SIZE,
        fmax(-HALF_BOUNDING_BOX_SIZE,
            particle->z + xorshiftNormalState(state, 0, standardDeviation)));
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
    return xorshiftDouble() < (1 / exp(changePotential / temperature));
}

int shouldMoveState(XorshiftState* state, double changePotential, double temperature) {
    return xorshiftDoubleState(state) < (1 / exp(changePotential / temperature));
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

void randomizeParticles(Particle* particles) {
    double* potentialsBlock = particles[0].potentials;
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
}

void randomizeParticlesState(XorshiftState* state, Particle* particles) {
    double* potentialsBlock = particles[0].potentials;
    // Initialize particles randomly in initial box
    for(int i = 0; i < TOTAL_PARTICLES; i++) {
        particles[i] = (Particle) {
            (xorshiftDoubleState(state) * BOUNDING_BOX_SIZE - HALF_BOUNDING_BOX_SIZE),
            (xorshiftDoubleState(state) * BOUNDING_BOX_SIZE - HALF_BOUNDING_BOX_SIZE),
            (xorshiftDoubleState(state) * BOUNDING_BOX_SIZE - HALF_BOUNDING_BOX_SIZE),
            &potentialsBlock[i * TOTAL_PARTICLES],
        };
    }
    // Initial potential calculation
    for(int i = 0; i < TOTAL_PARTICLES; i++) {
        recalculatePotential(particles, i);
    }
}

#if MULTI_THREADED == 1
SimulationResults simulateAnnealing(XorshiftState* state, int moveParticles, double temperature, const double temperatureCoolingFactor, const double particleCoolingFactor, const double standardDeviation) {
#else
SimulationResults simulateAnnealing(int moveParticles, double temperature, const double temperatureCoolingFactor, const double particleCoolingFactor, const double standardDeviation) {
#endif
    Particle particles[TOTAL_PARTICLES];
    double* potentialsBlock = (double*) calloc(TOTAL_PARTICLES * TOTAL_PARTICLES, sizeof(double));
    particles[0].potentials = potentialsBlock;

    #if MULTI_THREADED == 1
    randomizeParticlesState(state, particles);
    #else
    randomizeParticles(particles);
    #endif
    
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
    double bestBatchPotential = lastPotential;
    
    int totalBatches = 0;

    int acceptCount = 0;
    int convergingAcceptCount = 0; // Acceptance count during the converging process

    int particleComputations = 0;
    int convergingParticleComputations = 0;

    int strikes = 0;
    
    #if PARTICLE_GEOMETRIC_COOLING == 1
    double moveParticlesReal = particleCoolingFactor == 1.0 ? moveParticles - 1 : changeMoveParticles(moveParticles, particleCoolingFactor);
    #endif
    
    while(1) {
        totalBatches++;

        double currentBatchTotalPotential = 0.0;
        for(int i = 0; i < BATCH_SIZE; i++) {

            #if LARGE_MOVE_SIZE == 0
            unsigned long long changed = 0;
            for(int i = 0; i < moveParticles; i++) {
                int indexToMove;
                if(TOTAL_PARTICLES != moveParticles) {

                    #if MULTI_THREADED == 1
                    indexToMove = xorshiftDoubleState(state) * TOTAL_PARTICLES;
                    while(changed & (1ULL << indexToMove)) {
                        indexToMove = xorshiftDoubleState(state) * TOTAL_PARTICLES;
                    }
                    #else
                    indexToMove = xorshiftDouble() * TOTAL_PARTICLES;
                    while(changed & (1ULL << indexToMove)) {
                        indexToMove = xorshiftDouble() * TOTAL_PARTICLES;
                    }
                    #endif

                    changed = changed | (1ULL << indexToMove);
                } else {
                    indexToMove = i;
                }
                #if MULTI_THREADED == 1
                moveParticleState(state, &modifiedParticles[indexToMove], standardDeviation);
                #else
                moveParticle(&modifiedParticles[indexToMove], standardDeviation);
                #endif
                recalculatePotential(modifiedParticles, indexToMove);
            }
            #else
            unsigned long long changed = 0; // 1 = don't move, 0 = move
            for(int i = 0; i < TOTAL_PARTICLES - moveParticles; i++) {

                #if MULTI_THREADED == 1
                int indexToMove = xorshiftDoubleState(state) * TOTAL_PARTICLES;
                while(changed & (1ULL << indexToMove)) {
                    indexToMove = xorshiftDoubleState(state) * TOTAL_PARTICLES;
                }
                #else
                int indexToMove = xorshiftDouble() * TOTAL_PARTICLES;
                while(changed & (1ULL << indexToMove)) {
                    indexToMove = xorshiftDouble() * TOTAL_PARTICLES;
                }
                #endif
                
                changed = changed | (1ULL << indexToMove);
            }
            for(int i = 0; i < TOTAL_PARTICLES; i++) {
                if((changed & (1ULL << i)) == 0) {
                    #if MULTI_THREADED == 1
                    moveParticleState(state, &modifiedParticles[i], standardDeviation);
                    #else
                    moveParticle(&modifiedParticles[i], standardDeviation);
                    #endif
                    recalculatePotential(modifiedParticles, i);
                }
            }
            #endif
            
            particleComputations += moveParticles;

            double modifiedPotential = totalPotential(modifiedParticles);
            #if MULTI_THREADED == 1
            if(shouldMoveState(state, modifiedPotential - lastPotential, temperature)) {
            #else
            if(shouldMove(modifiedPotential - lastPotential, temperature)) {
            #endif
                copyOverAllParticles(modifiedParticles, particles);
                lastPotential = modifiedPotential;
                if(modifiedPotential < bestPotential) {
                    copyOverAllParticles(modifiedParticles, bestParticles);
                    bestPotential = modifiedPotential;
                }
                acceptCount++;
            } else {
                copyOverAllParticles(particles, modifiedParticles);
            }
            temperature = changeTemperature(temperature, temperatureCoolingFactor);

            #if PARTICLE_GEOMETRIC_COOLING == 1
            moveParticlesReal = changeMoveParticles(moveParticlesReal, particleCoolingFactor);
            moveParticles = 1 + floor(moveParticlesReal);
            #endif

            currentBatchTotalPotential += lastPotential;
        }
        double currentBatchPotential = currentBatchTotalPotential / BATCH_SIZE;
        if(currentBatchPotential >= bestBatchPotential) {
            strikes++;
            if(strikes == ALLOWED_STRIKES) break;
        } else {
            strikes = 0;
            bestBatchPotential = currentBatchPotential;
            convergingAcceptCount = acceptCount;
            convergingParticleComputations = particleComputations;
        }
    }
    
    free(potentialsBlock);
    free(modifiedPotentialsBlock);
    free(bestPotentialsBlock);
    return (SimulationResults) {
        bestPotential / 2, // We count each potential twice
        totalBatches - ALLOWED_STRIKES,
        (double) convergingAcceptCount / ((totalBatches - ALLOWED_STRIKES) * BATCH_SIZE),
        convergingParticleComputations
    };
}
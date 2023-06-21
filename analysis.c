#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include "src/random.h"
#include "src/lennard_jones.h"

int TOTAL_PARTICLES = 64;

int BATCH_SIZE = 1000;
int ALLOWED_STRIKES = 2;

double BOUNDING_BOX_SIZE = 4.0;
double HALF_BOUNDING_BOX_SIZE = 4.0/2;

int recordEvery = 5;
int recordSize = 20;

FILE* file;

typedef struct ConfigurationParameters {
    int count;
    XorshiftState state;
    int moveParticles;
    double particleFactor;
    double temperature;
    double temperatureFactor;
    double standardDeviation;
    double standardDeviationFactor;
    
    double* bestPotentialArray;
    double* convergingStepsArray;
    int* convergingParticleComputationsArray;
    double* convergingAcceptanceRateArray;
    
    double* bestPotentialRecord;
    double* acceptanceRateRecord;
} ConfigurationParameters;

void* simulateConfiguration(void* parametersP) {
    ConfigurationParameters parameters = *((ConfigurationParameters*) parametersP);

    int count = parameters.count;
    XorshiftState* state = &parameters.state;
    int moveParticles = parameters.moveParticles;
    double particleFactor = parameters.particleFactor;
    double temperature = parameters.temperature;
    double temperatureFactor = parameters.temperatureFactor;
    double standardDeviation = parameters.standardDeviation;
    double standardDeviationFactor = parameters.standardDeviationFactor;

    double* bestPotentialArray = parameters.bestPotentialArray;
    double* convergingStepsArray = parameters.convergingStepsArray;
    int* convergingParticleComputationsArray = parameters.convergingParticleComputationsArray;
    double* convergingAcceptanceRateArray = parameters.convergingAcceptanceRateArray;
    
    double* bestPotentialRecord = parameters.bestPotentialRecord;
    double* acceptanceRateRecord = parameters.acceptanceRateRecord;

    double* localBestPotentialRecord = (double*) calloc(recordSize, sizeof(double));
    double* localAcceptanceRateRecord = (double*) calloc(recordSize, sizeof(double));

    for(int i = 0; i < count; i++) {
        SimulationResults result = simulateAnnealing(state, moveParticles, particleFactor, temperature, temperatureFactor, standardDeviation, standardDeviationFactor, localBestPotentialRecord, localAcceptanceRateRecord);

        bestPotentialArray[i] = result.bestPotential;
        convergingStepsArray[i] = result.convergingSteps;
        convergingAcceptanceRateArray[i] = result.convergingAcceptanceRate;
        convergingParticleComputationsArray[i] = result.convergingParticleComputations;
        for(int i = 0; i < recordSize; i++) {
            bestPotentialRecord[i] += localBestPotentialRecord[i];
            acceptanceRateRecord[i] += localAcceptanceRateRecord[i];
        }

        memset(localBestPotentialRecord, 0, recordSize * sizeof(double));
        memset(localAcceptanceRateRecord, 0, recordSize * sizeof(double));
    }

    free(localBestPotentialRecord);
    free(localAcceptanceRateRecord);

    for(int i = 0; i < recordSize; i++) {
        bestPotentialRecord[i] /= count;
        acceptanceRateRecord[i] /= count;
    }

    pthread_exit(NULL);
}

int analysis(int moveParticles, double particleFactor, double temperature, double temperatureFactor, double standardDeviation, double standardDeviationFactor) {
    const int threadCount = 8;
    const int countPerThread = 250;
    const int totalCount = threadCount * countPerThread;

    setXorshiftState((uint64_t) 11992933392989292238ULL, (uint64_t) 995759136711242123ULL + time(0) * time(0) * time(0));
    
    char* fileName;
    asprintf(&fileName, "./results/analysis/decreasing_std/5k_record/%lf/particle_analysis_%d_%lf_%lf.txt", temperatureFactor, moveParticles, temperatureFactor, standardDeviationFactor);
    file = fopen(fileName, "a+");
    free(fileName);
    
    // Trying out a set of parameters
    fprintf(file, "Box size: %lf\n", BOUNDING_BOX_SIZE);
    fprintf(file, "Sample size: %d\n", totalCount);
    fprintf(file, "Batch size: %d\n", BATCH_SIZE);
    fprintf(file, "Allowed strikes: %d\n", ALLOWED_STRIKES);
    fprintf(file, "Total count: %d\n", totalCount);
    fprintf(file, "\n\n");
    
    pthread_t threads[threadCount];
    ConfigurationParameters* resources[threadCount];
    
    // Finding statistics
    double* bestPotentialArray = (double*) calloc(totalCount, sizeof(double));
    double* convergingStepsArray = (double*) calloc(totalCount, sizeof(double));
    int* convergingParticleComputationsArray = (int*) calloc(totalCount, sizeof(int));
    double* convergingAcceptanceRateArray = (double*) calloc(totalCount, sizeof(double));
    
    double* bestPotentialRecord = (double*) calloc(recordSize, sizeof(double));
    double* acceptanceRateRecord = (double*) calloc(recordSize, sizeof(double));
    
    for(int i = 0; i < threadCount; i++) {
        ConfigurationParameters* parametersP = (ConfigurationParameters*) malloc(sizeof(ConfigurationParameters));
        for(int z = 0; z < 13; z++) xorshift64();
        parametersP->state.x[0] = xorshift64();
        for(int z = 0; z < 17; z++) xorshift64();
        parametersP->state.x[1] = xorshift64();

        parametersP->count = countPerThread;
        parametersP->moveParticles = moveParticles;
        parametersP->particleFactor = particleFactor;
        parametersP->temperature = temperature;
        parametersP->temperatureFactor = temperatureFactor;
        parametersP->standardDeviation = standardDeviation;
        parametersP->standardDeviationFactor = standardDeviationFactor;
    
        parametersP->bestPotentialArray = &bestPotentialArray[i * countPerThread];
        parametersP->convergingStepsArray = &convergingStepsArray[i * countPerThread];
        parametersP->convergingParticleComputationsArray = &convergingParticleComputationsArray[i * countPerThread];
        parametersP->convergingAcceptanceRateArray = &convergingAcceptanceRateArray[i * countPerThread];

        parametersP->bestPotentialRecord = (double*) calloc(recordSize, sizeof(double));
        parametersP->acceptanceRateRecord = (double*) calloc(recordSize, sizeof(double));

        pthread_create(&threads[i], NULL, &simulateConfiguration, parametersP);
        resources[i] = parametersP;
    }

    for(int i = 0; i < threadCount; i++) {
        pthread_join(threads[i], NULL);
        for(int j = 0; j < recordSize; j++) {
            bestPotentialRecord[j] += resources[i]->bestPotentialRecord[j];
            acceptanceRateRecord[j] += resources[i]->acceptanceRateRecord[j];
        }
        free(resources[i]->bestPotentialRecord);
        free(resources[i]->acceptanceRateRecord);
        free(resources[i]);
    }
    
    for(int i = 0; i < recordSize; i++) {
        bestPotentialRecord[i] /= threadCount;
        acceptanceRateRecord[i] /= threadCount;
    }
    
    double bestPotentialMean = 0.0;
    double convergingStepsMean = 0.0;
    double convergingParticleComputationsMean = 0.0;
    double convergingAcceptanceRateMean = 0.0;
    for(int i = 0; i < totalCount; i++) {
        bestPotentialMean += bestPotentialArray[i];
        convergingStepsMean += convergingStepsArray[i];
        convergingParticleComputationsMean += convergingParticleComputationsArray[i];
        convergingAcceptanceRateMean += convergingAcceptanceRateArray[i];
    }
    bestPotentialMean /= totalCount;
    convergingStepsMean /= totalCount;
    convergingParticleComputationsMean /= totalCount;
    convergingAcceptanceRateMean /= totalCount;

    double bestPotentialSampleSTD = 0.0;
    double convergingStepsSampleSTD = 0.0;
    double convergingAcceptanceRateSampleSTD = 0.0;
    double convergingParticleComputationsSampleSTD = 0.0;
    for(int i = 0; i < totalCount; i++) {
        bestPotentialSampleSTD += pow(bestPotentialArray[i] - bestPotentialMean, 2);
        convergingStepsSampleSTD += pow(convergingStepsArray[i] - convergingStepsMean, 2);
        convergingParticleComputationsSampleSTD += pow(convergingParticleComputationsArray[i] - convergingParticleComputationsMean, 2);
        convergingAcceptanceRateSampleSTD += pow(convergingAcceptanceRateArray[i] - convergingAcceptanceRateMean, 2);
    }
    bestPotentialSampleSTD = sqrt(bestPotentialSampleSTD / (totalCount - 1));
    convergingStepsSampleSTD = sqrt(convergingStepsSampleSTD / (totalCount - 1));
    convergingParticleComputationsSampleSTD = sqrt(convergingParticleComputationsSampleSTD / (totalCount - 1));
    convergingAcceptanceRateSampleSTD = sqrt(convergingAcceptanceRateSampleSTD / (totalCount - 1));

    fprintf(file, "Starting move particles count: %d\n", moveParticles);
    fprintf(file, "Move particles cooling factor: %.9lf\n", particleFactor);
    fprintf(file, "----------\n");
    fprintf(file, "Starting temperature: %.9lf\n", temperature);
    fprintf(file, "Temperature cooling factor: %.9lf\n", temperatureFactor);
    fprintf(file, "----------\n");
    fprintf(file, "Starting standard deviation: %.9lf\n", standardDeviation);
    fprintf(file, "Standard deviation cooling factor: %.9lf\n", standardDeviationFactor);
    fprintf(file, "----------\n");
    fprintf(file, "Best potentials mean: %.9lf\n", bestPotentialMean);
    fprintf(file, "Best potentials sample STD: %.9lf\n", bestPotentialSampleSTD);
    fprintf(file, "----------\n");
    fprintf(file, "Converging steps mean: %.9lf\n", convergingStepsMean);
    fprintf(file, "Converging steps sample STD: %.9lf\n", convergingStepsSampleSTD);
    fprintf(file, "----------\n");
    fprintf(file, "Converging particle computations mean: %.9lf\n", convergingParticleComputationsMean);
    fprintf(file, "Converging particle computations sample STD: %.9lf\n", convergingParticleComputationsSampleSTD);
    fprintf(file, "----------\n");
    fprintf(file, "Converging acceptance rate mean: %.9lf\n", convergingAcceptanceRateMean);
    fprintf(file, "Converging acceptance rate sample STD: %.9lf\n", convergingAcceptanceRateSampleSTD);
    fprintf(file, "\n\n");
    
    free(bestPotentialArray);
    free(convergingStepsArray);
    free(convergingParticleComputationsArray);
    free(convergingAcceptanceRateArray);
    
    // Intermediate information every 10,000 steps. CSV format.
    // Best potential, acceptance rate
    for(int i = 0; i < recordSize - 1; i++) {
        fprintf(file, "%lf,", bestPotentialRecord[i]);
    }
    fprintf(file, "%lf\n", bestPotentialRecord[recordSize - 1]);
    for(int i = 0; i < recordSize - 1; i++) {
        fprintf(file, "%lf,", acceptanceRateRecord[i]);
    }
    fprintf(file, "%lf\n", acceptanceRateRecord[recordSize - 1]);
    
    free(bestPotentialRecord);
    free(acceptanceRateRecord);

    fclose(file);
    return 0;
}

int main() {
    int moveParticlesTests[] = { 1 };
    size_t moveParticlesTestsSize = sizeof(moveParticlesTests) / sizeof(int);
    double particleFactor = 1.0;
    double temperature = 25;
    double temperatureFactorTests[] = { 0.998800 };
    size_t temperatureFactorTestsSize = sizeof(temperatureFactorTests) / sizeof(double);
    double standardDeviation = 0.1;
    double standardDeviationFactorTests[] = { 0.999653, 0.999700, 0.999768, 0.999800, 0.999827, 0.999850, 0.999866, 0.999884, 0.999913, 0.999931 };
    size_t standardDeviationFactorTestsSize = sizeof(standardDeviationFactorTests) / sizeof(double);
    for(int i = 0; i < moveParticlesTestsSize; i++) {
        int moveParticles = moveParticlesTests[i];
        for(int j = 0; j < temperatureFactorTestsSize; j++) {
            double temperatureFactor = temperatureFactorTests[j];
            for(int k = 0; k < standardDeviationFactorTestsSize; k++) {
                double standardDeviationFactor = standardDeviationFactorTests[k];
                analysis(moveParticles, particleFactor, temperature, temperatureFactor, standardDeviation, standardDeviationFactor);
            }
        }
    }
}
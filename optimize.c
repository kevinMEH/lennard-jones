#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "src/random.h"
#include "src/lennard_jones.h"

int TOTAL_PARTICLES = 64;
int MOVE_PARTICLES = 1;

int BATCH_SIZE = 1000;
int ALLOWED_STRIKES = 3;

double BOUNDING_BOX_SIZE = 4.0;
double HALF_BOUNDING_BOX_SIZE = 4.0/2;

double temperature = 25;
double standardDeviation = 0.1;
double cooling_factor = 0.999;

int main() {
    setXorshiftState((uint64_t) 11992933392989292238ULL, (uint64_t) 995759136711242123ULL + time(0) * time(0) * time(0));
    const int sampleSize = 1000;
    
    char* fileName;
    asprintf(&fileName, "./results/particle.txt");
    FILE* file = fopen(fileName, "a+");
    free(fileName);
    
    // // Trying out a set of parameters
    // fprintf(file, "Box size: %lf\n", BOUNDING_BOX_SIZE);
    // fprintf(file, "Temperature: %lf\n", temperature);
    // fprintf(file, "Standard deviation: %lf\n", standardDeviation);
    // fprintf(file, "Sample size: %d\n", sampleSize);
    // fprintf(file, "\n\n");

    for(int i = 41; i <= 64; i++) {
        MOVE_PARTICLES = i;

        // Finding statistics
        double* bestPotentialArray = (double*) calloc(sampleSize, sizeof(double));
        double* convergingBatchesArray = (double*) calloc(sampleSize, sizeof(double));
        double* convergingAcceptanceRateArray = (double*) calloc(sampleSize, sizeof(double));

        double bestBestPotential = 9999999999;
        double bestPotentialMean = 0.0;
        double convergingBatchesMean = 0.0;
        double convergingAcceptanceRateMean = 0.0;
        for(int i = 0; i < sampleSize; i++) {
            SimulationResults result = simulateAnnealing(temperature, cooling_factor, standardDeviation);
            double bestPotential = result.bestPotential;
            double convergingBatches = result.convergingBatches;
            double convergingAcceptanceRate = result.convergingAcceptanceRate;
            bestPotentialArray[i] = bestPotential;
            convergingBatchesArray[i] = convergingBatches;
            convergingAcceptanceRateArray[i] = convergingAcceptanceRate;
            bestPotentialMean += bestPotential;
            convergingBatchesMean += convergingBatches;
            convergingAcceptanceRateMean += convergingAcceptanceRate;
            if(bestPotential < bestBestPotential) {
                bestBestPotential = bestPotential;
            }
        }
        bestPotentialMean /= sampleSize;
        convergingBatchesMean /= sampleSize;
        convergingAcceptanceRateMean /= sampleSize;
        
        double bestPotentialSampleSTD = 0.0;
        double convergingBatchesSampleSTD = 0.0;
        double convergingAcceptanceRateSampleSTD = 0.0;
        for(int i = 0; i < sampleSize; i++) {
            bestPotentialSampleSTD += pow(bestPotentialArray[i] - bestPotentialMean, 2);
            convergingBatchesSampleSTD += pow(convergingBatchesArray[i] - convergingBatchesMean, 2);
            convergingAcceptanceRateSampleSTD += pow(convergingAcceptanceRateArray[i] - convergingAcceptanceRateMean, 2);
        }
        bestPotentialSampleSTD = sqrt(bestPotentialSampleSTD / (sampleSize - 1));
        convergingBatchesSampleSTD = sqrt(convergingBatchesSampleSTD / (sampleSize - 1));
        convergingAcceptanceRateSampleSTD = sqrt(convergingAcceptanceRateSampleSTD / (sampleSize - 1));

        fprintf(file, "Move particles: %d\n", MOVE_PARTICLES);
        fprintf(file, "Best best potential: %lf\n", bestBestPotential);
        fprintf(file, "Best potentials mean: %lf\n", bestPotentialMean);
        fprintf(file, "Best potentials sample STD: %lf\n", bestPotentialSampleSTD);
        fprintf(file, "Converging batches mean: %lf\n", convergingBatchesMean);
        fprintf(file, "Converging batches sample STD: %lf\n", convergingBatchesSampleSTD);
        fprintf(file, "Converging acceptance rate mean: %lf\n", convergingAcceptanceRateMean);
        fprintf(file, "Best potentials sample STD: %lf\n", convergingAcceptanceRateSampleSTD);
        fprintf(file, "\n");
        
        fflush(file);
        
        free(bestPotentialArray);
        free(convergingBatchesArray);
        free(convergingAcceptanceRateArray);
    }

    fclose(file);

    return 0;
}
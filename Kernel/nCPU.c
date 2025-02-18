#include <stdio.h>
#include <stdlib.h>

#include "WolframLibrary.h"
#include "WolframNumericArrayLibrary.h"


#include <math.h>

#define GROUPSIZE 256

void solveFP( //cancel reflections
  double *dest, 
  double *src, 
  double *meta,
  mint size,
  int group_id,
  int cycles

) {
    for (int k=0; k<cycles; ++k) {
        const int work_item_id = group_id + GROUPSIZE * k;
        if (work_item_id >= size || work_item_id < 2) continue;

        const double CONSTF = 6.28331 * meta[1];    //thikcness * 6.28331
    
        const double freq = src[3 * work_item_id ]; //freqs
        const double t    = src[3 * work_item_id + 1]; //abs of transmission
        const double ph   = src[3 * work_item_id + 2]; //phase
    
        const double n = dest[5 * work_item_id ]; //n (working)
        const double k = dest[5 * work_item_id + 1]; //k (working)
    
        // Perform the optimized mathematical operations
        double var40 = k * k;
        double var41 = 1.0f + n;
        double var54 = var41 * var41;
        double var57 = var40 + var54;
        double var85 = n * n;
        double var90 = 2.0f * CONSTF * freq * n;
        double var81 = -2.0f * CONSTF * freq * k;
        double var82 = exp(var81);  // Exp[var81] -> exp(var81)
        double var58 = 1.0f / (var57 * var57);
        double var93 = -1.0f + var40 + var85;
        double var91 = cos(var90);  // Cos[var90] -> cos(var90)
        double var83 = -2.0f + k;
        double var84 = var83 * k;
        double var86 = -1.0f + var84 + var85;
        double var87 = 2.0f + k;
        double var88 = k * var87;
        double var89 = -1.0f + var88 + var85;
        double var94 = sin(var90);  // Sin[var90] -> sin(var90)
        double var80 = var57 * var57;
        double var99 = var40 + (-1.0f + n) * (-1.0f + n);
        double var98 = exp(2.0f * CONSTF * freq * k);  // Exp[...] -> exp(...)
        double var92 = 4.0f * k * var93 * var94;
    
        // Compute abs value (magnitude)
        double abs_val = sqrt(
            var58 * ((var99 / var98) * (var99 / var98) + var80 + var82 * (-2.0f * var86 * var89 * var91 + 2.0f * var92))
        );
    
        // Compute arg value (phase)
        double term1 = var82 * var58 * (4.0f * k * var93 * var91 + var86 * var89 * var94);
        double term2 = var82 * var58 * (var98 * var80 - var86 * var89 * var91 + var92);
        double arg_val = atan2(term1, term2);  // Use atan2 for phase calculation
    
        // Handle non-finite values (OpenCL doesn't have a direct equivalent to NumericQ or FiniteQ)
        if (!isfinite(abs_val) || !isfinite(arg_val)) {
            abs_val = 1.0f;
            arg_val = 0.0f;
        }    
    
        if (t * abs_val > 1.2) {
            dest[work_item_id*5 + 2] = 1.0f;
            dest[work_item_id*5 + 3] = 0.0f; 
        } else {
            dest[work_item_id*5 + 2] = t * abs_val; //modified amplitude of transmission
            dest[work_item_id*5 + 3] = ph - arg_val;  //modified phase
        }
    }
}

void savekValues(
   double *dest, 
  mint size,
  int group_id,
  int cycles
) {

    for (int k=0; k<cycles; ++k) {
        const int work_item_id = group_id + GROUPSIZE * k;
        if (work_item_id >= size) continue;

        dest[5*work_item_id + 4] = dest[5*work_item_id + 1];
    }    
}

void  movingAverage(
  double *dest, 
  mint kernelsize,
  mint size,
  int group_id,
  int cycles  
) { //Average n and k-values
    
    for (int k=0; k<cycles; ++k) {
        int work_item_id = group_id + GROUPSIZE * k;
        if (work_item_id >= size - kernelsize  || work_item_id < 4) continue;

        double accumulatorN = 0.0;
        double accumulatorK = 0.0;

        for (int u=0; u<kernelsize; u++) {
            accumulatorN = accumulatorN + dest[5*(work_item_id+u) + 0];
            accumulatorK = accumulatorK + dest[5*(work_item_id+u) + 1];
        }

        dest[5*work_item_id] = (accumulatorN) / ((double)kernelsize);
        dest[5*work_item_id + 1] = (accumulatorK) / ((double)kernelsize);    
    }
}

void  solveNK( //high-persision n and k values extraction
   double *dest, 
   double *src, 
   double *meta,
  mint size,
  mint iterations,
  int group_id,
  int cycles  
) {

    for (int k=0; k<cycles; ++k) {
        const int work_item_id = group_id + GROUPSIZE * k;
        if (work_item_id >= size) continue;

        const double ft = 1.0f / meta[1]; // 1/thickness (in cm-1)

        const double logT = log(dest[5*work_item_id + 2] * meta[2]); //amplitude
        const double ph = dest[5*work_item_id + 3] + meta[3]; //phase
        const double freq = src[3*work_item_id]; //freqs

        const double fT = ft / freq;

        const double n0 = 1.0f + (0.029979f*meta[0]*ft); //DC limit for n

        double np = dest[5*work_item_id]; //current n value
        double kp = dest[5*work_item_id + 1]; //current k value
        double n=np;
        double k=kp; 
        double modulus, arg, denominator, n2, re, im;

        for(int j = 0; j < iterations; j++) {
            n2 = 1.0f + np;
            n2 = n2 * n2;

            denominator = 1.0f / (kp * kp + n2);
            denominator = denominator * denominator;

            re = denominator * (np * n2 + kp * kp * (2.0f + np));
            im = denominator * (kp * (kp * kp + np * np - 1.0f));

            modulus = sqrt(re * re + im * im);
            arg = atan2(im, re);  // Equivalent to Arg[I im + re]

            n = 1.0f + 0.159152f * fT * (ph - arg);
            k = -0.159152f * fT * (logT - log(4.0f * modulus));
            //n = 1.0f + 0.159152f * fT * (ph - arg);
            //k = -0.159152f * fT * (logT - log(4.0f * modulus));

            np = n;
            kp = k;
        }

        // Validity check for n and k
        if (isnan(k) || isnan(n) || n < 1.0f) {
            k = 0.0f;
            n = n0;
        }

        dest[5*work_item_id] = n; //update n
        dest[5*work_item_id + 1] = k; //update k

    }
}

void  initialize( //rough approximation of n,k using transmission data
   double *dest, 
   double *src, 
   double *meta,
  mint    size,
  int group_id,
  int cycles  
) {

    const double ft = 1.0f / meta[1]; // 1/thickness in (cm-1)
    
    for (int j=0; j<cycles; ++j) {
        const int work_item_id = group_id + GROUPSIZE * j;
        if (work_item_id >= size) continue;

        if (work_item_id < 2) { // low-freq cutoff
            const double n0 = 1.0f + (0.029979f*meta[0] / meta[1]);

            dest[5*work_item_id] = n0;
            dest[5*work_item_id + 1] = 0.0f;
            dest[5*work_item_id + 2] = src[3*work_item_id+1];
            dest[5*work_item_id + 3] = src[3*work_item_id+2];
            dest[5*work_item_id + 4] = 0.0f;

            return;
        } 

        const double t = src[3*work_item_id + 1]; //amplitude
        const double p = src[3*work_item_id + 2]; //phase
        const double f = src[3*work_item_id]; //freqs

        dest[5*work_item_id] = 1.0f + 0.159152f * (p + meta[3]) * ft / f; //n
        dest[5*work_item_id + 1] = - 0.159152f * log(t * meta[2]) * ft / f; //kappa (k)
        dest[5*work_item_id + 2] = t; //amplitude of transmission
        dest[5*work_item_id + 3] = p; //phase
        dest[5*work_item_id + 4] = 0.0f; //save for later
    }
}


void clRun(
   double *dest, 
   double *src, 
   double *meta,
  int itemSize,
  int groupSize,
  int iterationsNK,
  int movAvg,
  int iterationsFP,
  int sid,
  int group_id
) {
    //int sid = get_group_id(0);

    //int group_id = get_local_id(0);
    int cycles = (int)ceil(((float)itemSize)/((float)GROUPSIZE));


     double *localDest = &dest[sid * (itemSize * 5)];
     double *localMeta = &meta[sid * 5];

    for (int lid=0; lid < GROUPSIZE; ++lid) initialize(localDest, src, localMeta, itemSize, lid, cycles);
    for (int lid=0; lid < GROUPSIZE; ++lid) solveNK(localDest, src, localMeta, itemSize, iterationsNK, lid, cycles);
    for (int lid=0; lid < GROUPSIZE; ++lid) movingAverage(localDest, movAvg+1, itemSize, lid, cycles);
    for (int lid=0; lid < GROUPSIZE; ++lid) savekValues(localDest, itemSize, lid, cycles);

    for (int h=0; h<iterationsFP; ++h) {
        for (int lid=0; lid < GROUPSIZE; ++lid) solveFP(localDest, src, localMeta, itemSize, lid, cycles);
        for (int lid=0; lid < GROUPSIZE; ++lid) solveNK(localDest, src, localMeta, itemSize, iterationsNK, lid, cycles);
        for (int lid=0; lid < GROUPSIZE; ++lid) movingAverage(localDest, movAvg+1, itemSize, lid, cycles);
    }

} 




DLLEXPORT mint WolframLibrary_getVersion() {
    return WolframLibraryVersion;
}

DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
    return LIBRARY_NO_ERROR;
}

DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) {
    return;
}


DLLEXPORT int run(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {    

    int err = LIBRARY_NO_ERROR;
   
    mint address_dest = MArgument_getInteger(Args[0]);
    double *dest = libData->MTensor_getRealData(*((MTensor*)address_dest));

    mint address_src = MArgument_getInteger(Args[1]);
    double *src = libData->MTensor_getRealData(*((MTensor*)address_src));    

    mint address_meta = MArgument_getInteger(Args[2]);
    double *meta = libData->MTensor_getRealData(*((MTensor*)address_meta));    
    
    mint itemSize = MArgument_getInteger(Args[3]);
    mint groupSize = MArgument_getInteger(Args[4]);
    mint iterationsNK = MArgument_getInteger(Args[5]);
    mint movAvg = MArgument_getInteger(Args[6]);
    mint iterationsFP = MArgument_getInteger(Args[7]);

    mint totalRunners = MArgument_getInteger(Args[8]);

    for (int sid=0; sid < groupSize; ++sid) {
        // {
            
            clRun(dest, src, meta, itemSize, groupSize, iterationsNK, movAvg, iterationsFP, sid, 0);
        //}
    }

    MArgument_setInteger(Res, 0);

    return err;
}


DLLEXPORT int load(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {    

    int err = LIBRARY_NO_ERROR;

    MTensor *ptr = (MTensor*)malloc(sizeof(MTensor));
    *ptr = MArgument_getMTensor(Args[0]);

    mint address = (mint)ptr;

    MArgument_setInteger(Res, address);
    return err;
}

DLLEXPORT int get(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {    

    int err = LIBRARY_NO_ERROR;
   
    mint address = MArgument_getInteger(Args[0]);
    MTensor *ptr = (MTensor*)address;
    MArgument_setMTensor(Res, *ptr);

    return err;
}

DLLEXPORT int unload(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {    

    int err = LIBRARY_NO_ERROR;
   
    mint address = MArgument_getInteger(Args[0]);
    MTensor *ptr = (MTensor*)address;
    //libData->MTensor_disownAll(*ptr);
    libData->MTensor_free(*ptr);
    free(ptr);

    MArgument_setInteger(Res, 0);

    return err;
}
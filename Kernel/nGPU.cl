#define GROUPSIZE 256

void __attribute__((always_inline)) solveFP( //cancel reflections
  __global float *dest, 
  __global float *src, 
  __global float *meta,
  mint size,
  int group_id,
  int cycles

) {
    for (int k=0; k<cycles; ++k) {
        const int work_item_id = group_id + GROUPSIZE * k;
        if (work_item_id >= size || work_item_id < 2) continue;

        const float CONSTF = 6.28331 * meta[1];    //thikcness * 6.28331
    
        const float freq = src[3 * work_item_id ]; //freqs
        const float t    = src[3 * work_item_id + 1]; //abs of transmission
        const float ph   = src[3 * work_item_id + 2]; //phase
    
        const float n = dest[5 * work_item_id ]; //n (working)
        const float k = dest[5 * work_item_id + 1]; //k (working)
    
        // Perform the optimized mathematical operations
        float var40 = k * k;
        float var41 = 1.0f + n;
        float var54 = var41 * var41;
        float var57 = var40 + var54;
        float var85 = n * n;
        float var90 = 2.0f * CONSTF * freq * n;
        float var81 = -2.0f * CONSTF * freq * k;
        float var82 = exp(var81);  // Exp[var81] -> exp(var81)
        float var58 = 1.0f / (var57 * var57);
        float var93 = -1.0f + var40 + var85;
        float var91 = cos(var90);  // Cos[var90] -> cos(var90)
        float var83 = -2.0f + k;
        float var84 = var83 * k;
        float var86 = -1.0f + var84 + var85;
        float var87 = 2.0f + k;
        float var88 = k * var87;
        float var89 = -1.0f + var88 + var85;
        float var94 = sin(var90);  // Sin[var90] -> sin(var90)
        float var80 = var57 * var57;
        float var99 = var40 + (-1.0f + n) * (-1.0f + n);
        float var98 = exp(2.0f * CONSTF * freq * k);  // Exp[...] -> exp(...)
        float var92 = 4.0f * k * var93 * var94;
    
        // Compute abs value (magnitude)
        float abs_val = sqrt(
            var58 * ((var99 / var98) * (var99 / var98) + var80 + var82 * (-2.0f * var86 * var89 * var91 + 2.0f * var92))
        );
    
        // Compute arg value (phase)
        float term1 = var82 * var58 * (4.0f * k * var93 * var91 + var86 * var89 * var94);
        float term2 = var82 * var58 * (var98 * var80 - var86 * var89 * var91 + var92);
        float arg_val = atan2(term1, term2);  // Use atan2 for phase calculation
    
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

void __attribute__((always_inline)) savekValues(
  __global float *dest, 
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

void __attribute__((always_inline)) movingAverage(
  __global float *dest, 
  mint kernelsize,
  mint size,
  int group_id,
  int cycles  
) { //Average n and k-values
    
    for (int k=0; k<cycles; ++k) {
        const int work_item_id = group_id + GROUPSIZE * k;
        if (work_item_id >= size - kernelsize) continue;

        float accumulatorN = 0.0f;
        float accumulatorK = 0.0f;
        for (int u=0; u<kernelsize; ++u) {
            accumulatorN = accumulatorN + dest[5*(work_item_id+u) + 0];
            accumulatorK = accumulatorK + dest[5*(work_item_id+u) + 1];
        }

        dest[5*work_item_id] = (accumulatorN) / ((float)kernelsize);
        dest[5*work_item_id + 1] = (accumulatorK) / ((float)kernelsize);        
    }
}

void __attribute__((always_inline)) solveNK( //high-persision n and k values extraction
  __global float *dest, 
  __global float *src, 
  __global float *meta,
  mint size,
  mint iterations,
  int group_id,
  int cycles  
) {

    for (int k=0; k<cycles; ++k) {
        const int work_item_id = group_id + GROUPSIZE * k;
        if (work_item_id >= size) continue;

        const float ft = 1.0f / meta[1]; // 1/thickness (in cm-1)

        const float logT = log(dest[5*work_item_id + 2] * meta[2]); //amplitude
        const float ph = dest[5*work_item_id + 3] + meta[3]; //phase
        const float freq = src[3*work_item_id]; //freqs

        const float fT = ft / freq;

        const float n0 = 1.0f + (0.029979f*meta[0]*ft); //DC limit for n

        float np = dest[5*work_item_id]; //current n value
        float kp = dest[5*work_item_id + 1]; //current k value
        float n=np;
        float k=kp; 
        float modulus, arg, denominator, n2, re, im;

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

void __attribute__((always_inline)) initialize( //rough approximation of n,k using transmission data
  __global float *dest, 
  __global float *src, 
  __global float *meta,
  mint    size,
  int group_id,
  int cycles  
) {

    const float ft = 1.0f / meta[1]; // 1/thickness in (cm-1)
    
    for (int j=0; j<cycles; ++j) {
        const int work_item_id = group_id + GROUPSIZE * j;
        if (work_item_id >= size) continue;

        if (work_item_id < 2) { // low-freq cutoff
            const float n0 = 1.0f + (0.029979f*meta[0] / meta[1]);

            dest[5*work_item_id] = n0;
            dest[5*work_item_id + 1] = 0.0f;
            dest[5*work_item_id + 2] = src[3*work_item_id+1];
            dest[5*work_item_id + 3] = src[3*work_item_id+2];
            dest[5*work_item_id + 4] = 0.0f;

            return;
        } 

        const float t = src[3*work_item_id + 1]; //amplitude
        const float p = src[3*work_item_id + 2]; //phase
        const float f = src[3*work_item_id]; //freqs

        dest[5*work_item_id] = 1.0f + 0.159152f * (p + meta[3]) * ft / f; //n
        dest[5*work_item_id + 1] = - 0.159152f * log(t * meta[2]) * ft / f; //kappa (k)
        dest[5*work_item_id + 2] = t; //amplitude of transmission
        dest[5*work_item_id + 3] = p; //phase
        dest[5*work_item_id + 4] = 0.0f; //save for later
    }
}

__kernel void clRun(
  __global float *dest, 
  __global float *src, 
  __global float *meta,
  mint itemSize,
  mint groupSize,
  mint iterationsNK,
  mint movAvg,
  mint iterationsFP
) {
    int sid = get_group_id(0);

    int group_id = get_local_id(0);
    int cycles = (int)ceil(((float)itemSize)/((float)GROUPSIZE));


    __global float *localDest = &dest[sid * (itemSize * 5)];
    __global float *localMeta = &meta[sid * 5];

    initialize(localDest, src, localMeta, itemSize, group_id, cycles);
    solveNK(localDest, src, localMeta, itemSize, iterationsNK, group_id, cycles);
    movingAverage(localDest, movAvg+1, itemSize, group_id, cycles);
    savekValues(localDest, itemSize, group_id, cycles);

    for (int h=0; h<iterationsFP; ++h) {
        solveFP(localDest, src, localMeta, itemSize, group_id, cycles);
        solveNK(localDest, src, localMeta, itemSize, iterationsNK, group_id, cycles);
        movingAverage(localDest, movAvg+1, itemSize, group_id, cycles);
    }

}





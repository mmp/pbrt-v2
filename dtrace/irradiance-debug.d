/* -*- mode: c++; -*- */
/*
  Prints a substantial amount of information about the execution of the
  irradiance caching integrator; can be useful for debugging.  Example
  output:

Starting interpolation @ p = (19.203182,0.000000,8.111011), n = (-0.000000,1.000000,0.000000)
  Checked sample at p=(15.200317,2.194881,9.090816), n=(0.034456,0.964192,-0.262958), E=(0.000000,0.000000,0.000000), max dist 2.438947, dir=(0.000000,0.000000,0.000000)
    dist 0.000000, perror 0.000000 nerror 0.000000
  Checked sample at p=(18.060753,0.000000,9.018047), n=(-0.000000,1.000000,0.000000), E=(0.878971,0.878232,0.882887), max dist 0.768360, dir=(-0.185814,4.260122,3.494656)
    dist 0.000000, perror 0.000000 nerror 0.000000
  Checked sample at p=(18.107153,-0.000000,7.707138), n=(-0.000000,1.000000,0.000000), E=(0.282019,0.282019,0.282019), max dist 1.104246, dir=(2.138046,1.835254,-0.559076)
    dist 0.000000, perror 0.000000 nerror 0.000000
  Checked sample at p=(18.618122,0.000000,7.757750), n=(-0.000000,1.000000,0.000000), E=(0.027848,0.027848,0.027848), max dist 0.735851, dir=(-0.061622,0.273127,-0.045451)
    dist 0.000000, perror 0.000000 nerror 0.000000
  Checked sample at p=(17.750380,0.000000,8.032465), n=(-0.000000,1.000000,0.000000), E=(0.141541,0.141541,0.141541), max dist 1.094631, dir=(0.463679,1.303884,-0.404309)
    dist 0.000000, perror 0.000000 nerror 0.000000
  Checked sample at p=(17.712950,0.000000,7.254876), n=(-0.000000,1.000000,0.000000), E=(0.129228,0.129196,0.131477), max dist 1.307988, dir=(0.922071,0.705974,0.622719)
    dist 0.000000, perror 0.000000 nerror 0.000000
  Checked sample at p=(19.046423,-0.000000,8.228955), n=(-0.000000,1.000000,0.000000), E=(0.850349,0.850349,0.850349), max dist 0.549654, dir=(-2.606634,5.070159,5.431772)
    dist 0.000000, perror 0.000000 nerror 0.000000

*/

#pragma D option quiet

struct Point { float x, y, z; };
struct Vector { float x, y, z; };
struct Normal { float x, y, z; };
struct Spectrum { float r, g, b; };
struct IrradianceSample {
    struct Spectrum E;
    struct Normal n;
    struct Point p;
    struct Vector wi;
    float maxDist;
};

/* const struct Point *, const struct Normal *, float maxDist, const struct Spectrum *E, const struct Vector *primaryDir, float pixel spacing */
:::irradiance_cache_added_new_sample {
    this->p = copyin(arg0, sizeof(struct Point));
    this->n = copyin(arg1, sizeof(struct Normal));
    /*    float this->maxDist; = arg2;*/
    this->E = copyin(arg3, sizeof(struct Spectrum));
    this->wi = copyin(arg4, sizeof(struct Vector));
    printf("Added sample at p=(%f,%f,%f), n=(%f,%f,%f), E=(%f,%f,%f), max dist %f, dir=(%f,%f,%f), pixel spacing %f\n",
           ((struct Point *)this->p)->x,
           ((struct Point *)this->p)->y, 
           ((struct Point *)this->p)->z,    
           ((struct Normal *)this->n)->x,
           ((struct Normal *)this->n)->y,
           ((struct Normal *)this->n)->z,  
           ((struct Spectrum *)this->E)->r,
           ((struct Spectrum *)this->E)->g, 
           ((struct Spectrum *)this->E)->b,    
           (float)arg2, /*this->maxDist,*/
           ((struct Vector *)this->wi)->x,
           ((struct Vector *)this->wi)->y, 
           ((struct Vector *)this->wi)->z,
           (float)arg5);
}

/* (const struct Point *p, const struct Normal *n) */
:::irradiance_cache_started_interpolation {
    this->p = copyin(arg0, sizeof(struct Point));
    this->n = copyin(arg1, sizeof(struct Normal));
    printf("Starting interpolation @ p = (%f,%f,%f), n = (%f,%f,%f)\n",
           ((struct Point *)this->p)->x,
           ((struct Point *)this->p)->y, 
           ((struct Point *)this->p)->z,    
           ((struct Normal *)this->n)->x,
           ((struct Normal *)this->n)->y,
           ((struct Normal *)this->n)->z);
}

/* (const struct Point *p, const struct Normal *n, int successful, int nfound) */
:::irradiance_cache_finished_interpolation {
}

/* (const struct IrradianceSample *, float dist, float perr, float nerr) */
:::irradiance_cache_checked_sample {
    this->samp = copyin(arg0, sizeof(struct IrradianceSample));
    printf("  Checked sample at p=(%f,%f,%f), n=(%f,%f,%f), E=(%f,%f,%f), max dist %f, dir=(%f,%f,%f)\n",
           ((struct IrradianceSample *)this->samp)->p.x,
           ((struct IrradianceSample *)this->samp)->p.y, 
           ((struct IrradianceSample *)this->samp)->p.z,    
           ((struct IrradianceSample *)this->samp)->n.x,
           ((struct IrradianceSample *)this->samp)->n.y,
           ((struct IrradianceSample *)this->samp)->n.z,  
           ((struct IrradianceSample *)this->samp)->E.r,
           ((struct IrradianceSample *)this->samp)->E.g, 
           ((struct IrradianceSample *)this->samp)->E.b,    
           ((struct IrradianceSample *)this->samp)->maxDist,
           ((struct IrradianceSample *)this->samp)->wi.x,
           ((struct IrradianceSample *)this->samp)->wi.y, 
           ((struct IrradianceSample *)this->samp)->wi.z);
    printf("    dist %f, perror %f nerror %f\n", (float)arg1, (float)arg2, (float)arg3);
}

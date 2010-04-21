/* -*- mode: c++; -*- */
/*
  Prints information about all of the irradiance samples created, in a
  format suitable for inclusion in a pbrt input file (e.g. for
  visualization).

  Example:

  AttributeBegin Translate 20.000000 6.807435 18.278511 Shape "sphere" "float radius" [.025] AttributeEnd
  AttributeBegin Translate 19.616219 5.935906 20.000000 Shape "sphere" "float radius" [.025] AttributeEnd

 */

#pragma D option quiet

struct Point { float x, y, z; };

:::irradiance_cache_added_new_sample {
    this->p = copyin(arg0, sizeof(struct Point));
    printf("\nAttributeBegin Translate ");
    printf("%f %f %f ", 
           ((struct Point *)this->p)->x,
           ((struct Point *)this->p)->y,
           ((struct Point *)this->p)->z);
    printf("Shape \"sphere\" \"float radius\" [.025] AttributeEnd\n");
}


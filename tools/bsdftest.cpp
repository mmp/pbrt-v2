// Kevin Egan

#include <stdio.h>
#include <stdlib.h>

#include "pbrt.h"
#include "reflection.h"
#include "montecarlo.h"
#include "memory.h"
#include "api.h"
#include "paramset.h"
#include "shapes/disk.h"

static MemoryArena arena;
static RNG rng;

// extract the red channel from a Spectrum class
double spectrumRedValue(const Spectrum & s)
{
    return ((float*) (& s))[0];
}

typedef void (*CreateBSDFFunc)(BSDF* bsdf);

void createBlinn0(BSDF* bsdf);
void createBlinn05(BSDF* bsdf);
void createBlinn2(BSDF* bsdf);
void createBlinn30and0(BSDF* bsdf);
void createAniso0_0(BSDF* bsdf);
void createAniso10_10(BSDF* bsdf);
void createAniso30_30(BSDF* bsdf);
void createLambertian(BSDF* bsdf);
void createOrenNayar0(BSDF* bsdf);
void createOrenNayar20(BSDF* bsdf);
void createFresnelBlend0(BSDF* bsdf);
void createFresnelBlend30(BSDF* bsdf);
void createPlastic(BSDF* bsdf);
void createSubstrate(BSDF* bsdf);


typedef void (*GenSampleFunc)(BSDF* bsdf,
    const Vector & wo, Vector* wi, float* pdf, Spectrum* f);

void Gen_Sample_f(BSDF* bsdf,
    const Vector & wo, Vector* wi, float* pdf, Spectrum* f);
void Gen_CosHemisphere(BSDF* bsdf,
    const Vector & wo, Vector* wi, float* pdf, Spectrum* f);
void Gen_UniformHemisphere(BSDF* bsdf,
    const Vector & wo, Vector* wi, float* pdf, Spectrum* f);


int main(int argc, char *argv[]) 
{
    Options opt;
    pbrtInit(opt);

    // number of monte carlo estimates
    //const int estimates = 1;
    const int estimates = 100000000;

    // radiance of uniform environment map
    const double environmentRadiance = 1.0;


    fprintf(stderr, "outgoing radiance from a surface viewed\n"
        "straight on with uniform lighting\n\n"
        "    uniform incoming radiance = %.3f\n"
        "    monte carlo samples = %d\n\n\n",
        environmentRadiance, estimates);


    CreateBSDFFunc BSDFFuncArray[] = {
        createBlinn0,
        createBlinn05,
        createBlinn2,
        createBlinn30and0,
        createAniso0_0,
        createAniso10_10,
        createAniso30_30,
//CO        createLambertian,
//CO        createOrenNayar0,
//CO        createOrenNayar20,
//CO        createFresnelBlend0,
//CO        createFresnelBlend30,
//CO        createPlastic,
//CO        createSubstrate,
    };

    const char* BSDFFuncDescripArray[] = {
        "Blinn (exponent 0)",
        "Blinn (exponent 0.5)",
        "Blinn (exponent 2)",
        "Blinn (exponent 30 and 0)",
        "Anisotropic (exponent 0, 0)",
        "Anisotropic (exponent 10, 10)",
        "Anisotropic (exponent 30, 30)",
//CO        "Lambertian",
//CO        "Oren Nayar (sigma 0)",
//CO        "Oren Nayar (sigma 20)",
//CO        "FresnelBlend (Blinn exponent 0)",
//CO        "FresnelBlend (Blinn exponent 30)",
//CO        "Plastic",
//CO        "Substrate",
    };

    GenSampleFunc SampleFuncArray[] = {
        Gen_Sample_f,
        Gen_CosHemisphere,
//CO        Gen_UniformHemisphere,
    };

    const char* SampleFuncDescripArray[] = {
        "BSDF Importance Sampling",
        "Cos Hemisphere",
//CO        "Uniform Hemisphere",
    };

    int numModels = sizeof(BSDFFuncArray) / sizeof(BSDFFuncArray[0]);
    int numModelsDescrip = sizeof(BSDFFuncDescripArray) /
        sizeof(BSDFFuncDescripArray[0]);
    int numGenerators = sizeof(SampleFuncArray) / sizeof(SampleFuncArray[0]);
    int numGeneratorsDescrip = sizeof(SampleFuncDescripArray) /
        sizeof(SampleFuncDescripArray[0]);

    if (numModels != numModelsDescrip) {
        fprintf(stderr, "BSDFFuncArray and BSDFFuncDescripArray out of sync!\n");
        exit(1);
    }

    if (numGenerators != numGeneratorsDescrip) {
        fprintf(stderr, "SampleFuncArray and SampleFuncDescripArray out of sync!\n");
        exit(1);
    }

    // for each bsdf model
    for (int model = 0; model < numModels; model++) {

        BSDF* bsdf;

        // create BSDF which requires creating a Shape, casting a Ray
        // that hits the shape to get a DifferentialGeometry object,
        // and passing the DifferentialGeometry object into the BSDF
        {
            Transform t = RotateX(-90);
            bool reverseOrientation = false;
            ParamSet p;

            Reference<Shape> disk = new Disk(new Transform(t), new Transform(Inverse(t)),
                                             reverseOrientation, 0., 1., 0, 360.);
            if (!disk) {
                fprintf(stderr, "Could not load disk plugin\n"
                    "  make sure the PBRT_SEARCHPATH environment variable is set\n");
                exit(1);
            }

            Point origin(0.1, 1, 0);  // offset slightly so we don't hit center of disk
            Vector direction(0, -1, 0);
            float tHit, rayEps;
            Ray r(origin, direction, 1e-3, INFINITY);
            DifferentialGeometry* dg = BSDF_ALLOC(arena, DifferentialGeometry)();
            disk->Intersect(r, &tHit, &rayEps, dg);

            bsdf = BSDF_ALLOC(arena, BSDF)(*dg, dg->nn);
            (BSDFFuncArray[model])(bsdf);
        }


        // facing directly at normal
        Vector woL = Normalize(Vector(0, 0, 1));
        Vector wo = bsdf->LocalToWorld(woL);
        const Normal &n = bsdf->dgShading.nn;

        // for each method of generating samples over the hemisphere
        for (int gen = 0; gen < numGenerators; gen++) {
            double redSum = 0.0;

            const int numHistoBins = 10;
            double histogram[numHistoBins][numHistoBins];
            for (int i = 0; i < numHistoBins; i++) {
                for (int j = 0; j < numHistoBins; j++) {
                    histogram[i][j] = 0;
                }
            }
            int badSamples = 0;
            int outsideSamples = 0;

            int warningTarget = 1;
            for (int sample = 0; sample < estimates; sample++) {
                Vector wi;
                float pdf;
                Spectrum f;

                // sample hemisphere around bsdf, wo is fixed
                (SampleFuncArray[gen])(bsdf, wo, & wi, & pdf, & f);

                double redF = spectrumRedValue(f);

                // add hemisphere sample to histogram
                Vector wiL = bsdf->WorldToLocal(wi);
                float x = Clamp(wiL.x, -1.f, 1.f);
                float y = Clamp(wiL.y, -1.f, 1.f);
                float wiPhi = (atan2(y, x) + M_PI) / (2.0 * M_PI);
                float wiCosTheta = wiL.z;
                bool validSample = (wiCosTheta > 1e-7);
                if (wiPhi < -0.0001 || wiPhi > 1.0001 || wiCosTheta > 1.0001) {
                    // wiCosTheta can be less than 0
                    fprintf(stderr, "bad wi! %.3f %.3f %.3f, (%.3f %.3f)\n",
                        wiL[0], wiL[1], wiL[2], wiPhi, wiCosTheta);
                } else if (validSample) {
                    int histoPhi      = (int) (wiPhi * numHistoBins);
                    int histoCosTheta = (int) (wiCosTheta * numHistoBins);
                    histogram[histoCosTheta][histoPhi] += 1.0 / pdf;
                }

                if (!validSample) {
                    outsideSamples++;
                } else if (pdf == 0.f || isnan(pdf) || redF < 0 || isnan(redF)) {
                    if (badSamples == warningTarget) {
                        fprintf(stderr, "warning %d, bad sample %d! "
                            "pdf: %.3f, redF: %.3f\n",
                            warningTarget, sample, pdf, redF);
                        warningTarget *= 10;
                    }
                    badSamples++;
                } else {
                    // outgoing radiance estimate =
                    //   bsdf * incomingRadiance * cos(wi) / pdf
                    redSum += redF * environmentRadiance * AbsDot(wi, n) / pdf;
                }
            }
            int goodSamples = estimates - badSamples;

            // print results
            fprintf(stderr, "*** BRDF: '%s', Samples: '%s'\n\n"
                "wi histogram showing the relative weight in each bin\n"
                "  all entries should be close to 2pi = %.5f:\n"
                "  (%d bad samples, %d outside samples)\n\n"
                "                          cos(theta) bins\n",
                BSDFFuncDescripArray[model], SampleFuncDescripArray[gen],
                M_PI * 2.0, badSamples, outsideSamples);
            double totalSum = 0.0;
            for (int i = 0; i < numHistoBins; i++) {
                fprintf(stderr, "  phi bin %02d:", i);
                for (int j = 0; j < numHistoBins; j++) {
                    fprintf(stderr, " %5.2f", histogram[i][j] *
                        numHistoBins * numHistoBins / goodSamples);
                    totalSum += histogram[i][j];
                }
                fprintf(stderr, "\n");
            }
            fprintf(stderr, "\n  final average :  %.5f (error %.5f)\n\n"
                "  radiance = %.5f\n\n",
                totalSum / goodSamples, totalSum / goodSamples - M_PI * 2.0,
                redSum / goodSamples);
        }
    }
    
    pbrtCleanup();
    return 0;
}









void Gen_Sample_f(BSDF* bsdf,
    const Vector & wo, Vector* wi, float* pdf, Spectrum* f)
{
    // only glossy or diffuse reflections (no specular reflections)
    BxDFType inflags = BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE | BSDF_GLOSSY);
    BxDFType outflags;
    BSDFSample sample(rng);
    *f = bsdf->Sample_f(wo, wi, sample, pdf, inflags, &outflags);

    // double check bsdf->Pdf() gives us the same answer
    Vector wiL = bsdf->WorldToLocal(*wi);
    float wiCosTheta = wiL.z;
    bool validSample = (wiCosTheta > 1e-7);
    
    if (validSample) {
        float verifyPdf = bsdf->Pdf(wo, *wi, inflags);
        if (fabs(verifyPdf - *pdf) > 1e-4) {
            fprintf(stderr, "BSDF::Pdf() doesn't match BSDF::Sample_f() !\n"
                "  Sample_f pdf %.3f, Pdf pdf %.3f\n"
                "  wo %.3f %.3f %.3f, wi %.3f %.3f %.3f\n",
                *pdf, verifyPdf, wo[0], wo[1], wo[2], (*wi)[0], (*wi)[1], (*wi)[2]);
            fprintf(stderr, "blah! validSample %d, wiCosTheta %.3f, wiL %.3f %.3f %.3f\n",
                validSample, wiCosTheta, wiL[0], wiL[1], wiL[2]);
        }
    }
}

void Gen_CosHemisphere(BSDF* bsdf,
    const Vector & wo, Vector* wi, float* pdf, Spectrum* f)
{
    float u1 = rng.RandomFloat();
    float u2 = rng.RandomFloat();
    Vector wiL = CosineSampleHemisphere(u1, u2);
    *wi = bsdf->LocalToWorld(wiL);
    float cosTheta = wiL.z;
    *pdf = CosineHemispherePdf(cosTheta, u2);

    *f = bsdf->f(wo, *wi);
}

void Gen_UniformHemisphere(BSDF* bsdf,
    const Vector & wo, Vector* wi, float* pdf, Spectrum* f)
{
    float u1 = rng.RandomFloat();
    float u2 = rng.RandomFloat();
    Vector wiL = UniformSampleHemisphere(u1, u2);
    *wi = bsdf->LocalToWorld(wiL);
    *pdf = UniformHemispherePdf();

    *f = bsdf->f(wo, *wi);
}









void createBlinn0(BSDF* bsdf)
{
    const float blinnExponent = 0.0;
    Spectrum Ks(1);
    MicrofacetDistribution* distribution = BSDF_ALLOC(arena, Blinn)(blinnExponent);
    Fresnel* fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
    BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution);
    bsdf->Add(bxdf);
}


void createBlinn05(BSDF* bsdf)
{
    const float blinnExponent = 0.5;
    Spectrum Ks(1);
    MicrofacetDistribution* distribution = BSDF_ALLOC(arena, Blinn)(blinnExponent);
    Fresnel* fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
    BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution);
    bsdf->Add(bxdf);
}


void createBlinn2(BSDF* bsdf)
{
    const float blinnExponent = 2.0;
    Spectrum Ks(1);
    MicrofacetDistribution* distribution = BSDF_ALLOC(arena, Blinn)(blinnExponent);
    Fresnel* fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
    BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution);
    bsdf->Add(bxdf);
}


void createBlinn30and0(BSDF* bsdf)
{
    const float blinnExponent1 = 30.0;
    Spectrum Ks(0.5);
    MicrofacetDistribution* distribution1 = BSDF_ALLOC(arena, Blinn)(blinnExponent1);
    Fresnel* fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
    BxDF* bxdf1 = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution1);
    bsdf->Add(bxdf1);

    const float blinnExponent2 = 0.0;
    MicrofacetDistribution* distribution2 = BSDF_ALLOC(arena, Blinn)(blinnExponent2);
    BxDF* bxdf2 = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution2);
    bsdf->Add(bxdf2);
}


void createAniso0_0(BSDF* bsdf)
{
    const float aniso1 = 0.0;
    const float aniso2 = 0.0;
    Spectrum Ks(1);
    MicrofacetDistribution* distribution =
        BSDF_ALLOC(arena, Anisotropic(aniso1, aniso2));
    Fresnel* fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
    BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution);
    bsdf->Add(bxdf);
}


void createAniso10_10(BSDF* bsdf)
{
    const float aniso1 = 10.0;
    const float aniso2 = 10.0;
    Spectrum Ks(1);
    MicrofacetDistribution* distribution =
        BSDF_ALLOC(arena, Anisotropic(aniso1, aniso2));
    Fresnel* fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
    BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution);
    bsdf->Add(bxdf);
}

void createAniso30_30(BSDF* bsdf)
{
    const float aniso1 = 30.0;
    const float aniso2 = 30.0;
    Spectrum Ks(1);
    MicrofacetDistribution* distribution =
        BSDF_ALLOC(arena, Anisotropic(aniso1, aniso2));
    Fresnel* fresnel = BSDF_ALLOC(arena, FresnelNoOp)();
    BxDF* bxdf = BSDF_ALLOC(arena, Microfacet)(Ks, fresnel, distribution);
    bsdf->Add(bxdf);
}


void createLambertian(BSDF* bsdf)
{
    Spectrum Kd(1);
    BxDF* bxdf = BSDF_ALLOC(arena, Lambertian)(Kd);
    bsdf->Add(bxdf);
}


void createOrenNayar0(BSDF* bsdf)
{
    Spectrum Kd(1);
    float sigma = 20.0;
    BxDF* bxdf = BSDF_ALLOC(arena, OrenNayar)(Kd, sigma);
    bsdf->Add(bxdf);
}


void createOrenNayar20(BSDF* bsdf)
{
    Spectrum Kd(1);
    float sigma = 0.0;
    BxDF* bxdf = BSDF_ALLOC(arena, OrenNayar)(Kd, sigma);
    bsdf->Add(bxdf);
}


void createFresnelBlend0(BSDF* bsdf)
{
    Spectrum d(0.5);
    Spectrum s(0.5);
    float exponent = 0.0;

    MicrofacetDistribution* distribution = BSDF_ALLOC(arena, Blinn)(exponent);
    BxDF* bxdf = BSDF_ALLOC(arena, FresnelBlend)(d, s, distribution);
    bsdf->Add(bxdf);
}


void createFresnelBlend30(BSDF* bsdf)
{
    Spectrum d(0.5);
    Spectrum s(0.5);
    float exponent = 30.0;

    MicrofacetDistribution* distribution = BSDF_ALLOC(arena, Blinn)(exponent);
    BxDF* bxdf = BSDF_ALLOC(arena, FresnelBlend)(d, s, distribution);
    bsdf->Add(bxdf);
}


void createPlastic(BSDF* bsdf)
{
    // Taken from plastic.cpp

    Spectrum kd(0.5);
    BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
    Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
    Spectrum ks = (0.5);
    float rough = 0.1;
    BxDF *spec = BSDF_ALLOC(arena, Microfacet)(ks, fresnel,
                                               BSDF_ALLOC(arena, Blinn)(1.f / rough));
    bsdf->Add(diff);
    bsdf->Add(spec);
}



void createSubstrate(BSDF* bsdf)
{
    // Taken from substrate.cpp

    Spectrum d(0.5);
    Spectrum s(0.5);
    float u = 0.1;
    float v = 0.1;

    bsdf->Add(BSDF_ALLOC(arena, FresnelBlend)(d, s,
                                              BSDF_ALLOC(arena, Anisotropic)(1.f/u, 1.f/v)));
}




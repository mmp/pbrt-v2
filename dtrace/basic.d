/* -*- mode: c++; -*- */
/*
  Prints basic set of statistics about rendering:
  Time on parsing, preprocessing, and rendering
  Number of shapes and transforms in the scene
  Number of Image Maps and memory used for them
  Total number of rays of various types (camera, shadow, etc) traced
*/

#pragma D option quiet

#define NS 1000000000
#define MS 1000000
#define NS_TO_SECS(x) ((x)/NS)
#define NS_TO_SECS_FRAC3(x) (((x) - NS * NS_TO_SECS(x)) / MS)
#define NS_TO_SECS_PARTS(x)  NS_TO_SECS(x), NS_TO_SECS_FRAC3(x)
#define NS_SECS_FMT "%8d.%03d"

#define SAFE_ZERO(x) ((x) > 0 ? (x) : 0xffffffffffffffff)

#define PCT_INT(x, y) ((100*(x))/SAFE_ZERO(y))
#define PCT_FRAC2(x, y) (((10000*(x))/SAFE_ZERO(y)) - 100 * PCT_INT(x,y))
#define PCT_PARTS(x, y) PCT_INT(x, y), PCT_FRAC2(x, y)
#define PCT_FMT "%3d.%02d%%"

#define RATIO_INT(x, y) ((x)/SAFE_ZERO(y))
#define RATIO_FRAC2(x, y) (((100*(x))/SAFE_ZERO(y)) - 100 * RATIO_INT(x,y))
#define RATIO_PARTS(x, y) RATIO_INT(x, y), RATIO_FRAC2(x, y)
#define RATIO_FMT "%d.%02dx"

uint64_t total_time;
uint64_t total_parse_time;
uint64_t total_preprocessing_time;
uint64_t total_rendering_time;

#define PARSING 0
#define PREPROCESSING 1
#define RENDERING 2
#define NUM_PHASES 3
int phase;

uint64_t phase_tasks_time[NUM_PHASES];
uint64_t phase_tasks_count[NUM_PHASES];

:::started_task {
    ++phase_tasks_count[phase];
    self->task_start_time = vtimestamp;
}

:::finished_task {
    phase_tasks_time[phase] += vtimestamp - self->task_start_time;
}

:::started_parsing {
    phase = PARSING;
    parsing = 1;
    parse_start_time = vtimestamp;
}

:::finished_parsing {
    parsing = 0;
    total_parse_time = vtimestamp - parse_start_time;
}

:::allocated_cached_transform {
    ++allocated_cached_transforms;
}

:::found_cached_transform {
    ++found_cached_transforms;
}

:::started_preprocessing {
    phase = PREPROCESSING;
    preprocessing = 1;
    preprocessing_start_time = vtimestamp;
}

:::finished_preprocessing {
    preprocessing = 0;
    total_preprocessing_time = vtimestamp - preprocessing_start_time;
}

:::started_rendering {
    phase = RENDERING;
    rendering = 1;
    rendering_start_time = vtimestamp;
}

:::finished_rendering {
    rendering = 0;
    total_rendering_time = vtimestamp - rendering_start_time;
}

:::started_rendertask {
    self->render_task_start_time = vtimestamp;
}

uint64_t imagemap_memory;

:::loaded_image_map {
    ++num_imagemaps_loaded;
    imagemap_memory += (arg1 * arg2 * arg3) / 1024;;
    @imagemap_memory_histogram = quantize((arg1 * arg2 * arg3) / 1024);
}

uint64_t total_rendertask_time;

:::finished_rendertask {
    this->end_time = vtimestamp;
    total_rendertask_time += this->end_time - self->render_task_start_time;
}

:::created_shape {
    ++num_shapes;
}

::started_camera_ray_integration {
    ++num_camera_rays;
}

:::started_ray_intersectionp {
    ++num_shadow_rays;
}

:::started_ray_intersection {
    ++num_general_rays;
}

:::started_specular_reflection_ray {
    ++num_reflection_rays;
}

:::started_specular_refraction_ray {
    ++num_refraction_rays;
}

:::started_bsdf_shading {
    ++num_points_shaded;
}

:::started_trilinear_texture_lookup {
    ++num_trilinear_texture_lookups;
}

:::started_ewa_texture_lookup {
    ++num_ewa_texture_lookups;
}

/* print it */

dtrace:::END {
    printf("Time breakdown:\n");
    printf("  Parsing total                       %12d.%03ds\n", 
           NS_TO_SECS_PARTS(total_parse_time));
    total_preprocessing_time += phase_tasks_time[PREPROCESSING];
    printf("  Preprocessing total                 %12d.%03ds\n",
           NS_TO_SECS_PARTS(total_preprocessing_time));
    total_rendering_time += total_rendertask_time;
    printf("  Rendering total                     %12d.%03ds\n", 
           NS_TO_SECS_PARTS(total_rendering_time));
    printf("\n");

    printf("Scene Details\n");
    printf("  Total shapes in scene                    %12d\n", num_shapes);
    printf("  Transforms allocated                     %12d (%d kB)\n",
           allocated_cached_transforms, (allocated_cached_transforms * 16 * 2 * 4)/1024);
    printf("  Transform cache hit rate                       %3d.%02d%%\n",
           PCT_PARTS(found_cached_transforms, found_cached_transforms + allocated_cached_transforms));
    printf("\n");

    printf("Rays Traced\n");
    printf("  Camera rays                              %12d\n", num_camera_rays);
    printf("  Shadow rays                              %12d\n", num_shadow_rays);
    printf("  Specular reflection rays                 %12d\n", num_reflection_rays);
    printf("  Specular refraction rays                 %12d\n", num_refraction_rays);
    printf("  Other general                            %12d\n", num_general_rays-num_camera_rays-num_reflection_rays-num_refraction_rays);
    printf("\n");

    printf("Image Maps\n");
    printf("  Number loaded                            %12d\n", num_imagemaps_loaded);
    printf("  Total memory used                        %12d kB\n", imagemap_memory);
    printa("  Memory used (kB)                         %@12u\n", @imagemap_memory_histogram);
}

/* -*- mode: c++; -*- */

#pragma D option quiet

:::started_camera_ray_integration {
    self->camera_ray_start = vtimestamp;
}

:::finished_camera_ray_integration {
    @camera_ray_time = quantize((vtimestamp - self->camera_ray_start));
}

dtrace:::END {
    printa("Time to trace camera rays (ns) %@u", @camera_ray_time);
}

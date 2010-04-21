/* -*- mode: c++; -*- */
/*
  prints the time that each renderng task takes to execute, in nanoseconds.
  also computes a histogram of the time each render task required
*/

#pragma D option quiet

:::started_rendertask,
:::mlt_started_rendering {
    self->render_task_start_time = vtimestamp;
}

:::finished_rendertask,
:::mlt_finished_rendering {
    printf("Render task time %d ns\n", vtimestamp - self->render_task_start_time);
    @task_time = quantize(vtimestamp - self->render_task_start_time);
}

//
// pbrtdisplay.cpp
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <GL/glfw.h>

int main() {
    if (glfwInit() != GL_TRUE) {
        fprintf(stderr, "Error initializing glfw!\n");
        return 1;
    }

#define MAX_RES 4096

    int xres = -1, yres = -1;

    int fd;
    for (;;) {
        fd = open("pbrt.display", O_RDONLY);
        if (fd <= 0) {
            usleep(30000);
            continue;
        }

        int imageSize = 4 * (2 + MAX_RES * MAX_RES);
        int *image = (int *)mmap(NULL, imageSize, PROT_READ, MAP_SHARED, fd, 0);
        if (image == MAP_FAILED) {
            perror("mmap");
            return 1;
        }

        if (image[0] != xres || image[1] != yres) {
            if (xres != -1)
                glfwCloseWindow();

            xres = image[0];
            yres = image[1];
            if (glfwOpenWindow(xres, yres, 8, 8, 8, 0, 0, 0, GLFW_WINDOW) != GL_TRUE) {
                fprintf(stderr, "Error opening window\n");
                return 1;
            }

            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            glOrtho(0, image[0]-1, 0, image[1]-1, 0., 1);
            glPixelZoom(1, -1);
            glRasterPos3f(0, image[1]-1, 0);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
        }

        glDrawPixels(image[0], image[1], GL_RGB, GL_UNSIGNED_BYTE, (void *)&image[2]);
        glfwSwapBuffers();

        usleep(30000);

        close(fd);
        munmap((void *)image, imageSize);
    }

    glfwTerminate();
    return 0;
}

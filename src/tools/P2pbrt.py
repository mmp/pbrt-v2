# Copyright Silvano Galliani 2014
# Compute camera parameters for pbrt given a camera matrix P
# Requires numpy and scipy
# Example:  python P2pbrt.py 0005.png.camera 3072 2048
# Where camera contains P matrix as follows
# 2246.17 -2208.64 -62.134 24477.4
# -316.161 -1091.08 2713.64 -8334.18
# -0.269944 -0.961723 -0.0471142 -7.01217

import numpy as np
import math
from scipy.linalg import rq
from numpy.linalg import det
from numpy.linalg import norm
import sys

if len(sys.argv) < 3:
    print("Usage: {} <camera_matrix.txt> <xresolution> <yresolution>".format(sys.argv[0]))
    sys.exit()

# Read input matrix
P = np.genfromtxt('0005.png.P')

# Save image resolution
nx = float(sys.argv[2])
ny = float(sys.argv[3])
#####


# Extract variable names to be consistent with Zisserman
M = P[:, :3]
p1 = P[:, 0]
p2 = P[:, 1]
p3 = P[:, 2]
p4 = P[:, 3]
m3 = M[2, :]

# Extract K and R (Zisserman p163)
K, R = rq(M)

# make diagonal of K positive
signK = np.sign(np.diag(K))

K = signK*K
R = signK*R

# Algebraically extract camera centre (p163 Zisserman)

C = np.array([det([p2, p3, p4]),
             -det([p1, p3, p4]),
              det([p1, p2, p4]),
             -det([p1, p2, p3])])
# Normalize by last component
C = C/C[3]
# Extract Euclidean vector
C = C[:3]

# Extract principal axis vector as in Zisserman page 161
v = det(M)*m3
# Normalize
v = v/norm(v)
# Translate vector to start at camera centre
v = C+v

# Extract focal lenght in x and y direction
fx = K[0, 0]
fy = K[1, 1]

# Extract off centered optical axis
cx = K[0, 2]
cy = K[1, 2]

# Compute aspect ratio and shifted screen size
aspect_ratio = nx/ny
screen_dimentionx = aspect_ratio * (fx/fy)
screen_dimentiony = 1.0

# Compute field of view
fov = 2.0*math.degrees(math.atan(ny/(2.0*fy)))

print('Film "image" "string filename" ["out.pfm"] "integer xresolution" [{:.0f}] \
        "integer yresolution" [{:.0f}]'.format(nx, ny))
print('SurfaceIntegrator "depth"')
print('LookAt')
print(' '.join(map(str, C)))
print(' '.join(map(str, v)))
print(' '.join(map(str, R[1, :])))
print
print('Camera "perspective" "float fov" [{}]'.format(fov))
print('"float screenwindow"')
print('[ {} {} {} {} ]'.format(
                             -(screen_dimentionx * cx) / (nx/2.0),
                              (screen_dimentionx * (nx-cx)) / (nx/2.0),
                             -(screen_dimentiony * cy) / (ny/2.0),
                              (screen_dimentiony * (ny-cy)) / (ny/2.0)))
print ('WorldBegin')
print ('AttributeBegin')
print ('#Include "fountain-resample.pbrt"')
print ('Include "fountain.pbrt"')
print ('AttributeEnd')
print ('WorldEnd')

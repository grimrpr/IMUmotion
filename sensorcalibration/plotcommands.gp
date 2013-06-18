set view equal xyz
splot "points_ellipsoid.txt" using 1:2:3 with points, \
      "points_retransformed.txt" using 1:2:3 with points, \
      "eigenvectors.txt" using 1:2:3:4:5:6 with vectors head filled lt 2
pause -1

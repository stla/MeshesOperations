library(MeshesOperations, lib.loc = "C:/SL/Rloclib")
library(rgl, lib.loc = "C:/SL/Rloclib")
library(rmarchingcubes)


f <- function(x, y, z, a, cosb, sinb){
    (sqrt((sqrt(x*x + (y*sinb + a*cosb)^2) - 2)^2) - 1)^2 +
      (sqrt((sqrt(z*z + (y*cosb - a*sinb)^2) - 2)^2) - 1)^2
}
a <- 0.6
b <- 0.785
cosb <- cos(b)
sinb <- sin(b)

x <- z <- seq(-3.5, 3.5, len = 150L)
y <- seq(-4.2, 4.2, len = 150L)
g <- expand.grid(X = x, Y = y, Z = z)
voxel <- array(
  with(g, f(X, Y, Z, a, cosb, sinb)),
  dim = c(150L, 150L, 150L)
)

contour_shape <- contour3d(
  griddata = voxel,
  level = 0.1,
  x = x,
  y = y,
  z = z
)

tmesh <- tmesh3d(
  vertices = t(contour_shape[["vertices"]]),
  indices = t(contour_shape[["triangles"]]),
  normals = contour_shape[["normals"]],
  homogeneous = FALSE
)
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
shade3d(tmesh, color = "orangered")

meshes <- MeshesOperations:::connectedComponentsK(
    rmesh0 = list(vertices = contour_shape[["vertices"]], faces = contour_shape[["triangles"]]),
    isTriangle = TRUE,
    triangulate = FALSE,
    clean = TRUE,
    normals = FALSE,
    epsilon = 0
)
library(MeshesOperations)
library(rgl)
library(randomcoloR)
meshes <- convexParts(mesh = NonConvexPolyhedron)
ncp <- length(meshes)
colors <- randomColor(ncp, hue = "random", luminosity = "bright")
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.8)
for(i in seq_len(ncp)){
  shade3d(toRGL(meshes[[i]]), color = colors[i])
}
plotEdges(
  NonConvexPolyhedron[["vertices"]],
  NonConvexPolyhedron[["edges"]]
)


# mesh one: a cube
mesh1 <- cube3d() # (from the rgl package)
# mesh two: another cube
mesh2 <- translate3d( # (from the rgl package)
  cube3d(), 1, 1, 1
)
# compute the union
umesh <- MeshesUnion(list(mesh1, mesh2))
cxp <- convexParts(mesh = umesh)


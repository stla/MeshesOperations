library(MeshesOperations)


mesh <- cyclideMesh(a = 1, c = 0.32, mu = 0.80, nu=350, nv=350)
points <- c(0, 0, 0.32)

distancesToMesh(mesh, points) # should be a - mu = 40

rgl::shade3d(mesh); rgl::points3d(rbind(c(0,0,0)), size=5)


library(MeshesOperations)
mesh <- rgl::cube3d()
points <- rbind(
  c(0, 0, 0),
  c(1, 1, 1)
)
distancesToMesh(mesh, points) # should be 1 and 0

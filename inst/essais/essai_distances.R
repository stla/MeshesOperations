
# cyclide example ####
library(MeshesOperations)
a <- 100; c <- 30; mu <- 80
mesh <- cyclideMesh(a, c, mu, nu = 100L, nv = 100L)
points <- c(c, 0, 0)
distancesToMesh(mesh, points) # should be a - mu = 20

rgl::shade3d(mesh); rgl::points3d(rbind(c(0,0,0)), size=5)


library(MeshesOperations)
mesh <- rgl::cube3d()
points <- rbind(
  c(0, 0, 0),
  c(1, 1, 1)
)
distancesToMesh(mesh, points) # should be 1 and 0

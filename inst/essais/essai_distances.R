library(MeshesOperations)

mesh <- cyclideMesh(a = 70, c = 32, mu = 60, nu=250, nv=250)
mesh <- rgl::cube3d()
points <- c(0, 0, 32)

distancesToMesh(mesh, points) # should be a - mu = 40

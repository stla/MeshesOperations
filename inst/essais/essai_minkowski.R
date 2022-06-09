library(MeshesOperations)
library(SurfaceReconstruction)
library(rgl)

psr <- PoissonReconstruction(SolidMobiusStrip)
shade3d(psr, color="yellow")

ms <- MinkowskiSum(scale3d(cube3d(),1/8,1/8,1/8), psr, normals = TRUE)
rms <- toRGL(ms)
shade3d(rms, color="yellow")
wire3d(rms)

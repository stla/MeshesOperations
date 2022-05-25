library(MeshesOperations, lib.loc = "C:/SL/Rloclib")
library(rgl, lib.loc = "C:/SL/Rloclib")
# plot
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("gainsboro")
#
palette <- trekcolors::trek_pal("klingon")
colors <- colorRampPalette(palette)(6)
invisible(lapply(1:5, function(i){
  shade3d(octahedraCompound[["rglmeshes"]][[i]], color = colors[i])
}))

# animation ####
movie3d(spin3d(axis = c(0, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE,
        startTime = 1/10,
        webshot = FALSE)


command <- "gifski --fps=10 --frames=zzpic*.png -o octahedraCompound.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)


stop("")

# compute the intersection of the 'gmp' meshes
inter <- MeshesIntersection(
  octahedraCompound[["gmpmeshes"]], numbersType = "gmp", clean = TRUE
)
inter2 <- Mesh(
  inter$gmpVertices, inter$faces, numbersType = "gmp", epsilon = 1e-8
)

# plot
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("#363940")
# first the five octahedra with transparency
invisible(lapply(
  octahedraCompound[["rglmeshes"]], shade3d,
  color = "yellow", alpha = 0.1
))
# now the intersection
rglinter <- tmesh3d(
  "vertices"    = t(inter[["vertices"]]),
  "indices"     = t(inter[["faces"]]),
  "homogeneous" = FALSE
)
shade3d(rglinter, color = "whitesmoke")
# and finally the edges
plotEdges(
  inter2[["vertices"]], inter2[["exteriorEdges"]],
  only = inter2[["exteriorVertices"]], color = "darkmagenta",
  verticesAsSpheres = FALSE
)
# this is an icosahedron

# animation ####
movie3d(spin3d(axis = c(0, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE,
        startTime = 1/10,
        webshot = FALSE)


command <- "gifski --fps=10 --frames=zzpic*.png -o octahedraCompound.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)


# we can recognize some exact values in the 'gmpVertices', for instance
#   the value at entry \code{(1, 1)} is \code{1/sqrt(3)}:
asNumeric(1L / gmpVertices[1L, 1L]^2L)


library(microbenchmark)

tetrahedraCompound <- list(

  "double_meshes" = lapply(meshes, function(mesh){
    list("vertices" = mesh[["vertices"]], faces = faces)
  }),

  "gmp_meshes" = lapply(meshes, function(mesh){
    list("vertices" = as.bigq(mesh[["vertices"]]), faces = faces)
  })

)

microbenchmark(
  exactKernel = MeshesIntersection(
    tetrahedraCompound[["double_meshes"]], numbersType = "lazyExact"
  ),
  gmpKernel = MeshesIntersection(
    tetrahedraCompound[["gmp_meshes"]], numbersType = "gmp"
  ),
  times = 5
)

rglmeshes <- tetrahedraCompound[["rglmeshes"]]
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
bg3d("#363940")
colors <- hcl.colors(5, palette = "Spectral")
invisible(lapply(
  1:5, function(i) shade3d(rglmeshes[[i]], color = colors[i])
))

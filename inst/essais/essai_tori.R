library(MeshesOperations, lib.loc = "C:/SL/Rloclib")
library(rgl, lib.loc = "C:/SL/Rloclib")

torusMesh <- function(R, r, S, s, ...){
  vertices <- matrix(NA_real_, nrow = 3L, ncol = S*s)
  Normals <- matrix(NA_real_, nrow = 3L, ncol = S*s)
  indices1 <- indices2 <- matrix(NA_integer_, nrow = 3L, ncol = S*s)
  u_ <- seq(0, 2*pi, length.out = S+1)[-1]
  v_ <- seq(0, 2*pi, length.out = s+1)[-1]
  for(i in 1:S){
    cos_ui <- cos(u_[i])
    sin_ui <- sin(u_[i])
    cx <- R * cos_ui
    cy <- R * sin_ui
    for(j in 1:s){
      rcos_vj <- r*cos(v_[j])
      n <- c(rcos_vj*cos_ui, rcos_vj*sin_ui, r*sin(v_[j]))
      Normals[, (i-1)*s+j] <- -n
      vertices[, (i-1)*s+j] <- c(cx,cy,0) + n
    }
  }
  # triangles
  s <- as.integer(s)
  for(i in 1L:S){
    ip1 <- ifelse(i==S, 1L, i+1L)
    for(j in 1L:s){
      jp1 <- ifelse(j==s, 1L, j+1L)
      indices1[,(i-1)*s+j] <- c((i-1L)*s+j, (i-1L)*s+jp1, (ip1-1L)*s+j)
      indices2[,(i-1)*s+j] <- c((i-1L)*s+jp1, (ip1-1L)*s+jp1, (ip1-1L)*s+j)
    }
  }
  out <- list(
    vb = rbind(vertices, 1),
    it = cbind(indices1, indices2),
    primitivetype = "triangle",
    material = list(...),
    normals = rbind(Normals,1)
  )
  class(out) <- c("mesh3d", "shape3d")
  out
}

tor1 <- torusMesh(3, 1, 16L, 16L)
tor2 <- translate3d(tor1, 0, 0, 5)

shade3d(tor1); shade3d(tor2)

vs1 <- t(tor1$vb[-4, ])
vs2 <- t(tor2$vb[-4, ])
vs <- rbind(vs1, vs2)
fs1 <- t(tor1$it[c(1,3,2), ])
fs2 <- t(tor2$it[c(1,2,3), ]) + max(fs1)
fs <- rbind(fs1, fs2)

mesh <- Mesh(vs, fs, normals = TRUE)
tmesh <- tmesh3d(
  vertices = t(mesh$vertices),
  indices = t(mesh$faces),
  normals = mesh$normals,
  homogeneous = FALSE
)

shade3d(tmesh, col = "red", back = "lines")

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

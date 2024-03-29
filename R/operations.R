#' @title Meshes intersection
#' @description Computes the intersection of the given meshes.
#'
#' @param meshes a list of meshes, each being either a
#'   \strong{rgl} mesh, or as a list with (at least) two fields:
#'   \code{vertices} and \code{faces}; the \code{vertices}
#'   matrix must have the \code{bigq} class if \code{numbersType="gmp"},
#'   otherwise it must be numeric
#' @param clean Boolean, whether to clean the input meshes (merging
#'   duplicated vertices, duplicated faces, removing isolated vertices)
#'   as well as the output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param numbersType the type of the numbers used in C++ for the
#'   computations; must be one of \code{"double"}, \code{"lazyExact"}
#'   (a type provided by CGAL for exact computations), or \code{"gmp"}
#'   (exact computations with rational numbers); of course using
#'   exact computations is slower but more accurate
#'
#' @return A triangle mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if \code{numbersType="gmp"}, and \code{normals} if \code{normals=TRUE}.
#'
#' @importFrom gmp as.bigq asNumeric
#'
#' @export
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#'
#' # mesh one: truncated icosahedron; we triangulate it for plotting
#' mesh1 <- Mesh(
#'   mesh = truncatedIcosahedron,
#'   triangulate = TRUE, normals = FALSE,
#'   numbersType = "lazyExact"
#' )
#'
#' # mesh two: a cube
#' cube <- translate3d( # (from the rgl package)
#'   cube3d(), 2, 0, 0
#' )
#' mesh2 <-
#'   list(vertices = t(cube[["vb"]][-4L, ]), faces = t(cube[["ib"]]))
#'
#' # compute the intersection
#' inter <- MeshesIntersection(list(mesh1, mesh2))
#'
#' # plot
#' rglmesh1 <- toRGL(mesh1)
#' rglinter <- toRGL(inter)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(rglmesh1, color = "yellow", alpha = 0.2)
#' shade3d(cube, color = "cyan", alpha = 0.2)
#' shade3d(rglinter, color = "red")
#' plotEdges(
#'   vertices = inter[["vertices"]], edges = inter[["exteriorEdges"]],
#'   edgesAsTubes = FALSE, lwd = 3, verticesAsSpheres = FALSE
#' )
#'
#' # other example, with 'gmp' rational numbers ####
#' library(MeshesOperations)
#' library(gmp)
#' library(rgl)
#'
#' cube <- cube3d()
#'
#' rglmesh1 <- cube
#' mesh1 <-
#'   list(vertices = t(cube[["vb"]][-4L, ]), faces = t(cube[["ib"]]))
#' mesh1[["vertices"]] <- as.bigq(mesh1[["vertices"]])
#'
#' rotMatrix <- t(cbind( # pi/3 around a great diagonal
#'   as.bigq(c(2, -1, 2), c(3, 3, 3)),
#'   as.bigq(c(2, 2, -1), c(3, 3, 3)),
#'   as.bigq(c(-1, 2, 2), c(3, 3, 3))
#' ))
#' mesh2 <-
#'   list(vertices = t(cube[["vb"]][-4L, ]), faces = t(cube[["ib"]]))
#' mesh2[["vertices"]] <- as.bigq(mesh2[["vertices"]]) %*% rotMatrix
#' rglmesh2 <- rotate3d(cube, pi/3, 1, 1, 1)
#'
#' inter <- MeshesIntersection(list(mesh1, mesh2), numbersType = "gmp")
#' # perfect vertices:
#' inter[["gmpVertices"]]
#' rglinter <- toRGL(inter)
#'
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#' bg3d("#363940")
#' shade3d(rglmesh1, color = "yellow", alpha = 0.2)
#' shade3d(rglmesh2, color = "orange", alpha = 0.2)
#' shade3d(rglinter, color = "hotpink")
#' plotEdges(
#'   inter[["vertices"]], inter[["exteriorEdges"]],
#'   only = inter[["exteriorVertices"]],
#'   color = "firebrick",
#'   tubesRadius = 0.05, spheresRadius = 0.07
#' )
MeshesIntersection <- function(
    meshes, clean = FALSE, normals = FALSE, numbersType = "double"
){
  numbersType <- match.arg(numbersType, c("double", "lazyExact", "gmp"))
  gmp <- numbersType == "gmp"
  stopifnot(is.list(meshes))
  checkMeshes <- lapply(meshes, function(mesh){
    if(inherits(mesh, "mesh3d")){
      vft  <- getVFT(mesh, beforeCheck = TRUE)
      mesh <- vft[["rmesh"]]
    }
    checkMesh(mesh[["vertices"]], mesh[["faces"]], gmp, aslist = TRUE)
  })
  areTriangle <- vapply(checkMeshes, `[[`, logical(1L), "isTriangle")
  triangulate <- !areTriangle
  meshes <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  if(numbersType == "double"){
    inter <- Intersection_K(meshes, clean, normals, triangulate)
  }else if(numbersType == "lazyExact"){
    inter <- Intersection_EK(meshes, clean, normals, triangulate)
  }else{ # numbersType == "gmp"
    inter <- Intersection_Q(meshes, clean, normals, triangulate)
  }
  if(gmp){
    vertices <- as.bigq(t(inter[["vertices"]]))
    inter[["gmpVertices"]] <- vertices
    vertices <- asNumeric(vertices)
  }else{
    vertices <- t(inter[["vertices"]])
  }
  inter[["vertices"]] <- vertices
  edgesDF <- inter[["edges"]]
  inter[["edgesDF"]] <- edgesDF
  inter[["edges"]] <- as.matrix(edgesDF[, c("i1", "i2")])
  exteriorEdges <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
  inter[["exteriorEdges"]] <- exteriorEdges
  inter[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  inter[["faces"]] <- t(inter[["faces"]])
  if(normals){
    inter[["normals"]] <- t(inter[["normals"]])
  }
  attr(inter, "toRGL") <- 3L
  class(inter) <- "cgalMesh"
  inter
}

#' @title Meshes difference
#' @description Computes the difference between two meshes.
#'
#' @param mesh1,mesh2 two meshes, each being either a
#'   \strong{rgl} mesh, or as a list with (at least) two fields:
#'   \code{vertices} and \code{faces}; the \code{vertices}
#'   matrix must have the \code{bigq} class if \code{numbersType="gmp"},
#'   otherwise it must be numeric
#' @param clean Boolean, whether to clean the input mesh (merging duplicated
#'   vertices, duplicated faces, removing isolated vertices) as well as the
#'   output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param numbersType the type of the numbers used in C++ for the
#'   computations; must be one of \code{"double"}, \code{"lazyExact"}
#'   (a type provided by CGAL for exact computations), or \code{"gmp"}
#'   (exact computations with rational numbers); of course using
#'   exact computations is slower but more accurate
#'
#' @return A triangle mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if \code{numbersType="gmp"}, and \code{normals} if \code{normals=TRUE}.
#'
#' @importFrom gmp as.bigq asNumeric
#'
#' @export
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#'
#' # mesh one: a cube
#' cube1 <- cube3d() # (from the rgl package)
#' mesh1 <-
#'   list(vertices = t(cube1[["vb"]][-4L, ]), faces = t(cube1[["ib"]]))
#'
#' # mesh two: another cube
#' cube2 <- translate3d( # (from the rgl package)
#'   cube3d(), 1, 1, 0
#' )
#' mesh2 <-
#'   list(vertices = t(cube2[["vb"]][-4L, ]), faces = t(cube2[["ib"]]))
#'
#' # compute the difference
#' differ <- MeshesDifference(mesh1, mesh2)
#'
#' # plot
#' rgldiffer <- toRGL(differ)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(cube1, color = "yellow", alpha = 0.2)
#' shade3d(cube2, color = "cyan", alpha = 0.2)
#' shade3d(rgldiffer, color = "red")
#' plotEdges(
#'   vertices = differ[["vertices"]], edges = differ[["exteriorEdges"]],
#'   edgesAsTubes = TRUE, verticesAsSpheres = TRUE
#' )
MeshesDifference <- function(
    mesh1, mesh2, clean = FALSE, normals = FALSE, numbersType = "double"
){
  stopifnot(is.list(mesh1))
  stopifnot(is.list(mesh2))
  numbersType <- match.arg(numbersType, c("double", "lazyExact", "gmp"))
  gmp <- numbersType == "gmp"
  if(inherits(mesh1, "mesh3d")){
    vft  <- getVFT(mesh1, beforeCheck = TRUE)
    mesh1 <- vft[["rmesh"]]
  }
  checkMesh1 <-
    checkMesh(mesh1[["vertices"]], mesh1[["faces"]], gmp, aslist = TRUE)
  triangulate1 <- !checkMesh1[["isTriangle"]]
  if(inherits(mesh2, "mesh3d")){
    vft  <- getVFT(mesh2, beforeCheck = TRUE)
    mesh2 <- vft[["rmesh"]]
  }
  checkMesh2 <-
    checkMesh(mesh2[["vertices"]], mesh2[["faces"]], gmp, aslist = TRUE)
  triangulate2 <- !checkMesh2[["isTriangle"]]
  mesh1 <- checkMesh1[c("vertices", "faces")]
  mesh2 <- checkMesh2[c("vertices", "faces")]
  if(numbersType == "double"){
    differ <-
      Difference_K(mesh1, mesh2, clean, normals, triangulate1, triangulate2)
  }else if(numbersType == "lazyExact"){
    differ <-
      Difference_EK(mesh1, mesh2, clean, normals, triangulate1, triangulate2)
  }else{
    differ <-
      Difference_Q(mesh1, mesh2, clean, normals, triangulate1, triangulate2)
  }
  if(gmp){
    vertices <- as.bigq(t(differ[["vertices"]]))
    differ[["gmpVertices"]] <- vertices
    vertices <- asNumeric(vertices)
  }else{
    vertices <- t(differ[["vertices"]])
  }
  differ[["vertices"]] <- vertices
  edgesDF <- differ[["edges"]]
  differ[["edgesDF"]] <- edgesDF
  differ[["edges"]] <- as.matrix(edgesDF[, c("i1", "i2")])
  exteriorEdges <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
  differ[["exteriorEdges"]] <- exteriorEdges
  differ[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  differ[["faces"]] <- t(differ[["faces"]])
  if(normals){
    differ[["normals"]] <- t(differ[["normals"]])
  }
  attr(differ, "toRGL") <- 3L
  class(differ) <- "cgalMesh"
  differ
}

#' @title Meshes union
#' @description Computes the union of the given meshes.
#'
#' @param meshes a list of two or more meshes, each being
#'   either a \strong{rgl} mesh, or as a list with (at least) two fields:
#'   \code{vertices} and \code{faces}; the \code{vertices}
#'   matrix must have the \code{bigq} class if \code{numbersType="gmp"},
#'   otherwise it must be numeric
#' @param clean Boolean, whether to clean the input meshes (merging duplicated
#'   vertices, duplicated faces, removed isolated vertices) as well as the
#'   output mesh
#' @param normals Boolean, whether to return the per-vertex normals of the
#'   output mesh
#' @param numbersType the type of the numbers used in C++ for the
#'   computations; must be one of \code{"double"}, \code{"lazyExact"}
#'   (a type provided by CGAL for exact computations), or \code{"gmp"}
#'   (exact computations with rational numbers); of course using
#'   exact computations is slower but more accurate
#'
#' @return A triangle mesh given as a list with fields \code{vertices},
#'   \code{faces}, \code{edges}, \code{exteriorEdges}, \code{gmpvertices}
#'   if \code{numbersType="gmp"}, and \code{normals} if \code{normals=TRUE}.
#'
#' @importFrom gmp as.bigq asNumeric
#'
#' @export
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#'
#' # mesh one: a cube
#' mesh1 <- cube3d() # (from the rgl package)
#'
#' # mesh two: another cube
#' mesh2 <- translate3d( # (from the rgl package)
#'   cube3d(), 1, 1, 1
#' )
#'
#' # compute the union
#' umesh <- MeshesUnion(list(mesh1, mesh2))
#'
#' # plot
#' rglumesh <- toRGL(umesh)
#' open3d(windowRect = c(50, 50, 562, 562))
#' shade3d(rglumesh, color = "red")
#' plotEdges(
#'   vertices = umesh[["vertices"]], edges = umesh[["exteriorEdges"]],
#'   edgesAsTubes = TRUE, verticesAsSpheres = TRUE
#' )
MeshesUnion <- function(
  meshes, clean = FALSE, normals = FALSE, numbersType = "double"
){
  stopifnot(is.list(meshes))
  stopifnot(length(meshes) >= 2)
  numbersType <- match.arg(numbersType, c("double", "lazyExact", "gmp"))
  gmp <- numbersType == "gmp"
  checkMeshes <- lapply(meshes, function(mesh){
    if(inherits(mesh, "mesh3d")){
      vft  <- getVFT(mesh, beforeCheck = TRUE)
      mesh <- vft[["rmesh"]]
    }
    checkMesh(mesh[["vertices"]], mesh[["faces"]], gmp, aslist = TRUE)
  })
  areTriangle <- vapply(checkMeshes, `[[`, logical(1L), "isTriangle")
  triangulate <- !areTriangle
  meshes <- lapply(checkMeshes, `[`, c("vertices", "faces"))
  if(numbersType == "double"){
    umesh <- Union_K(meshes, clean, normals, triangulate)
  }else if(numbersType == "lazyExact"){
    umesh <- Union_EK(meshes, clean, normals, triangulate)
  }else{
    umesh <- Union_Q(meshes, clean, normals, triangulate)
  }
  if(gmp){
    vertices <- as.bigq(t(umesh[["vertices"]]))
    umesh[["gmpVertices"]] <- vertices
    vertices <- asNumeric(vertices)
  }else{
    vertices <- t(umesh[["vertices"]])
  }
  umesh[["vertices"]] <- vertices
  edgesDF <- umesh[["edges"]]
  umesh[["edgesDF"]] <- edgesDF
  umesh[["edges"]] <- as.matrix(edgesDF[, c("i1", "i2")])
  exteriorEdges <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
  umesh[["exteriorEdges"]] <- exteriorEdges
  umesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
  umesh[["faces"]] <- t(umesh[["faces"]])
  if(normals){
    umesh[["normals"]] <- t(umesh[["normals"]])
  }
  attr(umesh, "toRGL") <- 3L
  class(umesh) <- "cgalMesh"
  umesh
}

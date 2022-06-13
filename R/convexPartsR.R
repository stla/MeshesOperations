#' Title
#'
#' @param vertices x
#' @param faces x
#' @param mesh x
#' @param triangulate x
#'
#' @return x
#' @export
#'
#' @examples
#' #
convexParts <- function(
  vertices, faces, mesh = NULL, triangulate = TRUE
){
  if(!is.null(mesh)){
    if(inherits(mesh, "mesh3d")){
      vft  <- getVFT(mesh, beforeCheck = TRUE)
      mesh <- vft[["rmesh"]]
    }
    vertices <- mesh[["vertices"]]
    faces    <- mesh[["faces"]]
  }
  checkedMesh <- checkMesh(vertices, faces, gmp = FALSE, aslist = TRUE)
  vertices         <- checkedMesh[["vertices"]]
  faces            <- checkedMesh[["faces"]]
  rmesh <- list("vertices" = vertices, "faces" = faces)
  cxparts <- convexDecomposition(rmesh, triangulate)
  ncp <- length(cxparts)
  meshes <- vector("list", ncp)
  for(i in seq_len(ncp)){
    mesh <- cxparts[[i]]
    mesh[["vertices"]] <- t(mesh[["vertices"]])
    edgesDF <- mesh[["edges"]]
    mesh[["edgesDF"]] <- edgesDF
    mesh[["edges"]] <- as.matrix(edgesDF[, c("i1", "i2")])
    exteriorEdges <- as.matrix(subset(edgesDF, exterior)[, c("i1", "i2")])
    mesh[["exteriorEdges"]] <- exteriorEdges
    mesh[["exteriorVertices"]] <- which(table(exteriorEdges) != 2L)
    # if(triangulate){
    #   edges0DF <- mesh[["edges0"]]
    #   mesh[["edges0DF"]] <- edges0DF
    #   mesh[["edges0"]] <- as.matrix(edges0DF[, c("i1", "i2")])
    #   if(normals){
    #     mesh[["normals0"]] <- t(mesh[["normals0"]])
    #   }
    # }
    if(FALSE && (triangulate || homogeneousFaces)){
      mesh[["faces"]] <- do.call(rbind, mesh[["faces"]])
      toRGL <- ifelse(triangulate, 3L, homogeneousFaces)
    }else{
      sizes <- lengths(mesh[["faces"]])
      usizes <- uniqueN(sizes)
      if(usizes == 1L){
        if(sizes[1L] %in% c(3L, 4L)){
          toRGL <- sizes[1L]
        }
        mesh[["faces"]] <- do.call(rbind, mesh[["faces"]])
      }else if(usizes == 2L && all(sizes %in% c(3L, 4L))){
        toRGL <- 34L
      }else{
        toRGL <- FALSE
      }
    }
    attr(mesh, "toRGL") <- toRGL
    class(mesh) <- "cgalMesh"
    meshes[[i]] <- mesh
  }
  meshes
}

#' @title Torus mesh
#' @description Triangle mesh of a torus (in \strong{rgl} format).
#'
#' @param R,r major and minor radii, positive numbers 
#' @param nu,nv numbers of subdivisions, integers (at least 3)
#' @param rgl Boolean, whether to return a \strong{rgl} mesh
#'
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}) if 
#'   \code{rgl=TRUE}, otherwise a \code{cgalMesh} list (vertices, faces, 
#'   and normals).
#' @export
#' 
#' @importFrom rgl tmesh3d
#'
#' @examples
#' library(MeshesOperations)
#' library(rgl)
#' mesh <- torusMesh(R = 3, r = 1)
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "green")
#' wire3d(mesh)
torusMesh <- function(R, r, nu = 50, nv = 30, rgl = TRUE){
  stopifnot(isPositiveNumber(R), isPositiveNumber(r))
  stopifnot(R > r)
  stopifnot(nu >= 3, nv >= 3)
  stopifnot(isBoolean(rgl))
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  nunv <- nu*nv
  vs      <- matrix(NA_real_, nrow = 3L, ncol = nunv)
  normals <- matrix(NA_real_, nrow = nu*nv, ncol = 3L)
  tris1   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  tris2   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  u_ <- seq(0, 2*pi, length.out = nu + 1L)[-1L]
  cosu_ <- cos(u_)
  sinu_ <- sin(u_)
  v_ <- seq(0, 2*pi, length.out = nv + 1L)[-1L]
  cosv_ <- cos(v_)
  sinv_ <- sin(v_)
  Rrcosv_ <- R + r*cosv_
  rsinv_ <- r*sinv_
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  for(i in 1L:(nu-1L)){
    i_nv <- i*nv
    rg <- (i_nv - nv + 1L):i_nv
    cosu_i <- cosu_[i]
    sinu_i <- sinu_[i]
    vs[, rg] <- rbind(
      cosu_i * Rrcosv_,
      sinu_i * Rrcosv_,
      rsinv_
    )
    normals[rg, ] <- cbind(
      cosu_i * cosv_,
      sinu_i * cosv_,
      sinv_
    )
    k1 <- i_nv - nv
    k_ <- k1 + j_
    l_ <- k1 + jp1_
    m_ <- i_nv + j_
    tris1[, k_] <- rbind(m_, l_, k_)
    tris2[, k_] <- rbind(m_, i_nv + jp1_, l_)
  }
  i_nv <- nunv
  rg <- (i_nv - nv + 1L):i_nv
  vs[, rg] <- rbind(
    Rrcosv_,
    0,
    rsinv_
  )
  normals[rg, ] <- cbind(
    cosv_,
    0,
    sinv_
  )
  k1 <- i_nv - nv
  l_ <- k1 + jp1_
  k_ <- k1 + j_
  tris1[, k_] <- rbind(j_, l_, k_)
  tris2[, k_] <- rbind(j_, jp1_, l_)
  if(rgl){
    tmesh3d(
      vertices = vs,
      indices = cbind(tris1, tris2),
      normals = normals, 
      homogeneous = FALSE
    )
  }else{
    out <- list(
      "vertices" = t(vs), "faces" = t(cbind(tris1, tris2)), "normals" = normals
    )
    attr(out, "toRGL") <- 3L
    class(out) <- "cgalMesh"
    out
  }
}

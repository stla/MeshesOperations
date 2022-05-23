# the compound of five octahedra ####
library(rgl)

# all vertices
phi <- (1+sqrt(5))/2
vertices <- rbind(
  c( 0.0 , 0.0 , -phi )
  , c( 0.0 , -phi , 0.0 )
  , c( -phi , 0.0 , 0.0 )
  , c( 0.0 , 0.0 , phi )
  , c( 0.0 , phi , 0.0 )
  , c( phi , 0.0 , 0.0 )
  , c( 0.5 , phi/2 , phi*phi/2 )
  , c( phi/2 , phi*phi/2 , 0.5 )
  , c( phi*phi/2 , 0.5 , phi/2 )
  , c( 0.5 , phi/2 , -phi*phi/2 )
  , c( phi/2 , -phi*phi/2 , 0.5 )
  , c( -phi*phi/2 , 0.5 , phi/2 )
  , c( 0.5 , -phi/2 , phi*phi/2 )
  , c( -phi/2 , phi*phi/2 , 0.5 )
  , c( phi*phi/2 , 0.5 , -phi/2 )
  , c( 0.5 , -phi/2 , -phi*phi/2 )
  , c( -phi/2 , -phi*phi/2 , 0.5 )
  , c( -phi*phi/2 , 0.5 , -phi/2 )
  , c( -0.5 , phi/2 , phi*phi/2 )
  , c( phi/2 , phi*phi/2 , -0.5 )
  , c( phi*phi/2 , -0.5 , phi/2 )
  , c( -0.5 , phi/2 , -phi*phi/2 )
  , c( phi/2 , -phi*phi/2 , -0.5 )
  , c( -phi*phi/2 , -0.5 , phi/2 )
  , c( -0.5 , -phi/2 , phi*phi/2 )
  , c( -phi/2 , phi*phi/2 , -0.5 )
  , c( phi*phi/2 , -0.5 , -phi/2 )
  , c( -0.5 , -phi/2 , -phi*phi/2 )
  , c( -phi/2 , -phi*phi/2 , -0.5 )
  , c( -phi*phi/2 , -0.5 , -phi/2 )
)
vertices1 <- vertices[1+c(0,1,2,3,4,5), ]
vertices2 <- vertices[1+c(6,10,14,27,25,23), ]
vertices3 <- vertices[1+c(7,11,12,28,26,21), ]
vertices4 <- vertices[1+c(8,9,13,29,24,22), ]
vertices5 <- vertices[1+c(15,16,17,18,19,20), ]

# the face indices are common to all octahedra
faces <- rbind(
  c(1L, 2L, 3L),
  c(1L, 5L, 6L),
  c(6L, 2L, 1L),
  c(1L, 3L, 5L),
  c(4L, 3L, 2L),
  c(6L, 4L, 2L),
  c(3L, 4L, 5L),
  c(6L, 5L, 4L)
)

# define the five octahedra meshes
mesh1 <- list(
  "vertices" = vertices1,
  "faces" = faces
)
mesh2 <- list(
  "vertices" = vertices2,
  "faces" = faces
)
mesh3 <- list(
  "vertices" = vertices3,
  "faces" = faces
)
mesh4 <- list(
  "vertices" = vertices4,
  "faces" = faces
)
mesh5 <- list(
  "vertices" = vertices5,
  "faces" = faces
)

# put them in a list to apply the intersection algorithm
meshes <- list(
  mesh1, mesh2, mesh3, mesh4, mesh5
)

octahedraCompound <- list(

  "meshes" = meshes,

  "rglmeshes" = list(
    tmesh3d(
      "vertices"    = t(vertices1),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices2),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices3),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices4),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(vertices5),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    )
  )
)

# MeshesOperations: operations on 3D meshes with R

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/MeshesOperations/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/MeshesOperations/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check-valgrind](https://github.com/stla/MeshesOperations/actions/workflows/R-CMD-check-valgrind.yaml/badge.svg)](https://github.com/stla/MeshesOperations/actions/workflows/R-CMD-check-valgrind.yaml)
<!-- badges: end -->

Using the C++ library **CGAL** to perform operations on 3D Meshes.


## Boolean operations on meshes

#### Intersection

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/Intersection.png)

#### Difference

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/Difference.png)

#### Union

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/Union.png)

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/tetrahedraCompound.gif)


## Shape smoothing

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/HopfTorusSmoothed.gif)

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/StanfordBunnySmoothed.gif)


## Minkowski sum

- Octahedron + sphere:

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/OctahedronPlusSphere.gif)

- Tetrahedron + truncated icosahedron:

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/TetrahedronPlusTruncatedIcosahedron.gif)


## Decomposition in convex parts

![](https://raw.githubusercontent.com/stla/MeshesOperations/master/inst/screenshots/pentagrammicPrism.png)


## License

This package is provided under the GPL-3 license. If you wish to use CGAL for 
commercial purposes, you must obtain a license from the 
[GeometryFactory](https://geometryfactory.com).


## Blog posts

I wrote a 
[blog post](https://laustep.github.io/stlahblog/posts/BooleanOpsOnMeshes.html)
devoted to Boolean operations on meshes (using **RCGAL**, the 
[ancestor](https://laustep.github.io/stlahblog/posts/splittingRCGAL.html) 
of **MehesOperations**), and 
[another one](https://laustep.github.io/stlahblog/posts/MinkowskiSumLeonardo.html) 
showing an example of the Minkowski addition.

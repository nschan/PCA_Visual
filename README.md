Visualizing more than 3 dimensions
================
Niklas Schandry
May / June 2020

# Motivation

I am teaching a course on metagenomic analysis, with the core audience
being molecular biologists. My goal is to explain data analysis that are
commonly use for these multivariate datasets, and how to perform these
using R. I am not intending to teach math, or pure stats, but
data-science. Part of this is finding low-dimensional representations of
a high-dimensional space, basically ordination.

In my experience, two dimensional spaces are easy to grasp, and can be
visualized easily as an xy-plot. Some people are able to visualize a
three dimensional space, but this can already be challenging for many
datapoints.

Literally everyone (myself included) is highly confused when trying to
imagine a four-dimensional space, using axis that point, well, into
three dimensional space. It is not possible using another spacial
dimension, because we cannot imagine more than three orthogonal spacial
dimensions, and trying to do this is our first instinct, as we are
already thinking in a “space-space” (? i am no mathematician)

However, anyone who does any type of analysis does in fact regularly
project additional dimensions onto a two-dimensional space, for example
when we color dots by a categorial or continous variable. We can use the
exact same principle, to add a fourth dimension to spheres in a
three-dimensional space.

To explain, we generate a matrix with independent, and less independent
dimensions, this matrix is centered.

``` r
library(tidyverse)
library(magrittr)
library(patchwork)
```

``` r
set.seed(238954)
Dim1 = rnorm(10,0,25)
Dim2 = rnorm(10,0,9)
Dim3 = rnorm(10,0,15)
Dim4 = Dim3 * (- 2.4) + rnorm(10,0,7)
cent_mat_2 <- cbind(Dim1,Dim2,Dim3,Dim4) %>% as.data.frame() %>% round(1)
```

``` r
cent_mat_2
```

    ##     Dim1  Dim2  Dim3  Dim4
    ## 1  -42.9 -13.8  18.6 -34.0
    ## 2   32.4  -7.1 -23.9  62.2
    ## 3    4.9  -3.3 -38.2  90.4
    ## 4   -2.2   1.5  -2.0  12.1
    ## 5  -14.9   6.1  12.2 -30.2
    ## 6  -28.1   6.2  12.9 -20.2
    ## 7   18.5  -6.9  -7.8  28.4
    ## 8   -1.5  12.2   4.1  -1.1
    ## 9    2.0  -1.6   3.6 -10.5
    ## 10  47.2   9.5  27.6 -69.1

We see that Dim1 has the largest variance of the first three (720.4),
followed by Dim3 (393.7). However, Dim4 widely exceeds the variance of
Dim1 by a factor of \~3.1 (see above)

``` r
cent_mat_2 %>% var %>% round(1)
```

    ##        Dim1   Dim2   Dim3   Dim4
    ## Dim1  720.4   41.4 -127.9  234.4
    ## Dim2   41.4   68.7   61.0 -157.7
    ## Dim3 -127.9   61.0  393.7 -935.2
    ## Dim4  234.4 -157.7 -935.2 2253.1

Initially, we only focus on the first three dimensions.

``` r
p1_2 <- cent_mat_2 %>% 
  ggplot() +
  geom_point(aes(x=Dim1,y=Dim2))

p1_3 <- cent_mat_2 %>% 
  ggplot() +
  geom_point(aes(x=Dim1,y=Dim3))

p2_3  <- cent_mat_2 %>% 
  ggplot() +
  geom_point(aes(x=Dim2,y=Dim3))

(p1_2 + p1_3 + p2_3)  & theme_minimal()
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

This repository contains a blender file (in /blender) in which i have
reproduced the three dimensional arrangement. This scene (in blender)
also contains one very strong light, and a camera, both are centered on
0,0,0. The spheres completely absorb light, and therefore appear as dots
in the camera view. Both, light and camera are fixed in their distance
to 0,0,0. There is an “axis” element (“Empty” in the right panel,
displayed as a black cross) which is at 0,0,0. Both the camera and the
light follow the rotation of this element. This allows us to freely
rotate around the center in three dimensional space, observing the
“projection” in the camera panel on the lower right.

![](imgs/Blender_scene_overview.PNG)

If we try to mimick a principal component analysis, we rotate (click on
the “empty” element and press R and then move around in the panels
showing the coordinate system, pressing R twice allows for a more “free”
rotation, which is cool and confusing) until we maximize distance along
the x axis of this panel, and then “fix” this and rotate to maximize the
distance along the other axis. It is not hard to find a projection that
approximates the output of a PCA pretty well, once we get a feel for how
rotation works in blender. This is me playing a bit in blender:

![](imgs/First_projection.PNG)

This is the PCA result:

``` r
pca_obj_1 <- cent_mat_2[,1:3] %>% pcaMethods::pca(center = F) 
pca_obj_1@scores %>%
  as.data.frame() %>%
  ggplot(aes(x = PC1, y =  PC2)) + 
  geom_point() +
  theme_minimal()
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- --> Now to
understand how we can project a fourth dimension, that is orthogonal to
the other three, we simply map dimension 4 to a color vector. I have
tried to explain this using the plots below, initially we map dim1 and
dim4 in space, then we color the dots by their dim4 value, and then, we
remove the axis that was used for dim4, retaining information about both
dimensions in a one-dimensional display.

``` r
(p1_4 + p1_4_color) + p1_4_color_no_y
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
p1_4_color_no_y
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

I extract these colors.

``` r
plotdat <- ggplot_build(p1_4_color_no_y)
colors <- plotdat$data[[1]] %$% colour
```

``` r
colors
```

    ##  [1] "#404988" "#8CD44D" "#FDE725" "#2B9289" "#3F5089" "#3A608B" "#2CAA82" "#2B7E8D" "#31708D"
    ## [10] "#440154"

In the blender, you will notice that there is a “collection” called
sample colors

\!(imgs/blender\_sample\_color\_collection.PNG)

This is the same samples, but colored according to their value, you can
display this by clicking on the eye-icon to make them visible.

We can now project their position in the fourth dimension as well. We
also know (see above) that dimension four explains much more variance.
We now try to project dimension four (from yellow to purple) along the X
axis in the camera view. Then, we rotate to maximize the distance along
the other axis. Since we take into account variance of the color vector
explicitly (and that of the other vectors implicitly since we “see” the
variance as the distance between points) we are actually able to
intuitively maximize these, the hardest part is gaining intuition for
how the blender rotation works. This is my projection:

\!(imgs/Second\_Projection.PNG)

Again, we can compare this to the PCA result on these four columns.

``` r
pca_obj_1 <- cent_mat_2 %>% pcaMethods::pca(center = F) 
pca_obj_1@scores %>%
  as.data.frame() %>%
  ggplot(aes(x = PC1, y =  PC2)) + 
  geom_point(color=colors)+
  theme_minimal()
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Once we understand this, we can actually imagine even more dimension,
using common mappings e.g. shape, or transparency to add additional
dimensions beyond the third.

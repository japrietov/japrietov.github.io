---
layout: post
title:  "Divide and Conquer!: The Warnock’s Algorithm"
date:   2017-09-27 12:00:00 -0500
categories: jekyll update
---

Computer Graphics attempts to represent objects in the general three-dimensional universe.
Most objects are not transparent and so we are interested in their outer surfaces, which have properties such as shape,
colour and texture which affect the graphical representation.
A wire-frame drawing of a solid object is less realistic because it includes parts of the object which are hidden in reality,
and this generates a need for some form of hidden-line or hidden-surface removal. One famous hidden-surface removal algorithm is the Warnock’s Algorithm,
it was first described in 1969 by John Warnock ([Warnock, 1969](#references)).
It solves the problem of rendering a complicated image by recursive subdivision of a scene until areas are obtained that are trivial to compute([Divide and conquer algorithm](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithm))


## Hidden surface removal (HSR)

Also known as **occlusion culling (OC)** or **visible surface determination (VSD)**, is the process used to determine
which surfaces and parts of surfaces are not visible from a certain viewpoint.We have looked at removing backfaces,
this stage will usually eliminate about half of the surfaces in a scene from further processing.
Once the polygons to be displayed have been clipped, we need to draw them onto the display device (Rasterization).

We also need to remove parts of those surfaces that are partially obscured by other surfaces. These are called _hidden surfaces_.
If hidden surfaces are not dealt with correctly, then we may get distorted images like the right of the image below.
_Rasterization_ and _HSR_ are usually combined into a single phase of the rendering pipeline.


|           No Lines Removed           |             Hidden Lines Removed  |             Hidden Surfaces Removed  |
|---------------------------------------------------------------|:-----------------------------------------------------:|-------------------------------------------------------|
|             ![surrounding](/images/noLinesRemoved.png)          |            ![intersecting](/images/hiddenLinesRemoved.png) |            ![intersecting](/images/hiddenSurfaceRemoved.png) |




Considering the rendering pipeline, the projection, the clipping, and the rasterization steps are handled differently by
the following algorithms ([Computer graphics lecture, 2007](#references)):

- [**Z-buffering**](https://en.wikipedia.org/wiki/Z-buffering) During rasterization the depth/Z value of each pixel (or sample in the case of anti-aliasing, but without
loss of generality the term pixel is used) is checked against an existing depth value.
If the current pixel is behind the pixel in the Z-buffer, the pixel is rejected, otherwise it is shaded and its depth value replaces the one in the Z-buffer.
Z-buffering supports dynamic scenes easily, and is currently implemented efficiently in graphics hardware. This is the current standard.


- [**Painter's algorithm**](https://en.wikipedia.org/wiki/Painter%27s_algorithm) sorts polygons by their barycenter and draws them back to front.
This produces few artifacts when applied to scenes with polygons of similar size forming smooth meshes and back-face culling turned on.
The cost here is the sorting step and the fact that visual artifacts can occur. This algorithm is broken by design for general scenes,
as it cannot handle polygons in various common configurations.


- [**Binary space partitioning (BSP)**](https://en.wikipedia.org/wiki/Binary_space_partitioning) divides a scene along planes corresponding to polygon boundaries.
The subdivision is constructed in such a way as to provide an unambiguous depth ordering from any point in the scene when the BSP tree is traversed.
The disadvantage here is that the BSP tree is created with an expensive pre-process. This means that it is less suitable for scenes consisting of dynamic geometry.
The advantage is that the data is pre-sorted and error free, ready for the previously mentioned algorithms. Note that the BSP is not a solution to HSR, only an aid.


Another hidden surface algorithm is the Warnock Algorithm, which speaks in more detail below.


## The Warnock's Algorithm
he Warnock algorithm is a _hidden surface algorithm_ invented by John Warnock that is typically used in the field of computer graphics ([Warnock, 1969](#references)).
It solves the problem of rendering a complicated image by recursive subdivision of a scene until areas are obtained that are trivial to compute.
In other words, if the scene is simple enough to compute efficiently then it is rendered;
otherwise it is divided into smaller parts which are likewise tested for simplicity or _“simple enough”_.

Some general requirements for _“simple enough”_ is that given an area of interest, classify polygons as

|                   |              |
|:-----------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------:|
|           **Surrounding (Fills the viewport)**          |             **Intersecting (Polygon partially visible)**  |
|             ![surrounding](/images/surrounding1.png )          |            ![intersecting](/images/intersecting1.png ) |
|            **Containing (Polygon completely visible)**            |             **Disjoint (Polygon invisible)**              |
|             ![containing](/images/containing.png )            |            ![disjoint](/images/disjoint1.png )         |




After classify the polygon is due:

- Disjoint polygons can be eliminated.
- Intersecting polygons can be split into disjoint and containing polygons.
- If there is only one polygon, fill area with background, then scan-fill polygon area.
- If there is a singe surrounding polygon, and no intersecting or containing polygons, display the surrounding polygons color.
- There is a surrounding polygon in front of all other polygons, display the surrounding polygons's color in the area.
- A surrounding polygon can be determined to be in front by computing the z coordinates of the surrounding, intersecting, and containing polygons at the four corners of the area.

Then the warnock algorithm approach is a set of steps described below
1. Take a given section of the screen (the entire screen, in the first pass).
2. Check to see if it is _"simple enough"_.
3. If it is, display it.
4. If it isn't, subdivide the screen into four sections and check each of the new sections (Step 1).

In other words, the inputs for Warnock algorithm are detail of polygons and a viewport. The good case is that if the
detail of polygons is very simple then creates the polygons in the viewport.
The continuous step is to divide the viewport into four equally sized quadrants and to recursively identify the algorithm for each quadrant,
with a polygon list changed such that it contains polygons that are detectable in that quadrant.

<p align="center">
  <img src="/images/warnockAlgorithm.png" />
</p>


<p align="center">
  A simple scene after apply the Warnock algorithm
</p>





This is a [divide and conquer algorithm](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithm) with run-time of ***O(np)***, where _n_ is the number of polygons and _p_ is the number of pixels in the viewport.

## References

Warnock, J. (1969). A Hidden Surface Algorithm for Computer Generated Halftone Pictures. Ph.D. Dissertation. The University of Utah, 32.

Computer graphics lecture. (2007). Hidden Surface Determination. Retrieved from https://courses.cs.washington.edu/courses/cse557/07wi/lectures/hidden-surfaces.pdf
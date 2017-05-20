Algorithm
-------
Fabs: the Forward and Backward Stagewise (Fabs) algorithm.


Maintainer
-------
Xingjie Shi  <xingjieshi@njue.edu.cn>


Publication
-------
Shi, X., Huang, J.,  Huang, Y.,  & Ma, S. (2017). A Forward and Backward Stagewise Algorithm for Nonconvex Loss Functions and Convex Penalties. Manuscript. 



Usage
-------
1. Demo
   a) examples.r provides examples of Fabs on calculating solution paths under the ordinary least square loss and smoothed partial rank loss.

2. Main functions:
   a) fabs.r              --- Fabs algorithm for adaptive lasso under the ordinary least square loss and smoothed partial rank loss
   b) fabsBrdige.r   ---  Fabs algorithm for Bridge penalty under the ordinary least square loss and smoothed partial rank loss

3. Functions used in fabs.r and fabsBrdige.r:
   a) standardize.r  --- standardize design matrix
   b) loss.r              --- loss functions
   c) penalty.r         --- penalty functions
   d) derivative.r     --- derivatives of loss functions 
   e) spr.r               --- interface functions for calculate smoothed partial rank loss and it derivatives

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/special_orthogonal_group.R
\docType{class}
\name{SpecialOrthogonalGroup-class}
\alias{SpecialOrthogonalGroup-class}
\alias{SpecialOrthogonalGroup}
\title{The special orthogonal group SO(n),
i.e. the Lie group of rotations in n dimensions.}
\description{
The special orthogonal group SO(n),
i.e. the Lie group of rotations in n dimensions.
}
\section{Methods}{

\describe{
\item{\code{Belongs(point)}}{Evaluate if a point belongs to SO(n).}

\item{\code{Compose(point.1, point.2)}}{Compose two elements of SO(n).}

\item{\code{GroupExpFromIdentity(tangent.vec)}}{Compute the group exponential of the tangent vector at the identity.}

\item{\code{GroupLogFromIdentity(point)}}{Compute the group logarithm of the point at the identity.}

\item{\code{Inverse(point)}}{Compute the group inverse in SO(n).}

\item{\code{JacobianTranslation(point, left.or.right = "left")}}{Compute the jacobian matrix of the differential
of the left/right translations from the identity to point in SO(n).}

\item{\code{MatrixFromRotationVector(rot.vec)}}{Convert rotation vector to rotation matrix.}

\item{\code{RandomUniform(n.samples = 1)}}{Sample in SO(n) with the uniform distribution.}

\item{\code{Regularize(point)}}{In 3D, regularize the norm of the rotation vector,
to be between 0 and pi, following the axis-angle
representation's convention.
If the angle angle is between pi and 2pi,
the function computes its complementary in 2pi and
inverts the direction of the rotation axis.}

\item{\code{RotationVectorFromMatrix(rot.mat)}}{In 3D, convert rotation matrix to rotation vector
      (axis-angle representation).

      Get the angle through the trace of the rotation matrix:
      The eigenvalues are:
      1, cos(angle) + i sin(angle), cos(angle) - i sin(angle)
      so that: trace = 1 + 2 cos(angle), -1 <= trace <= 3

      Get the rotation vector through the formula:
      S_r = angle / ( 2 * sin(angle) ) (R - R^T)

      For the edge case where the angle is close to pi,
      the formulation is derived by going from rotation matrix to unit
      quaternion to axis-angle:
      r = angle * v / |v|, where (w, v) is a unit quaternion.

      In nD, the rotation vector stores the n(n-1)/2 values of the
      skew-symmetric matrix representing the rotation.}

\item{\code{SkewMatrixFromVector(vec)}}{In 3D, compute the skew-symmetric matrix,
known as the cross-product of a vector,
associated to the vector vec.}

\item{\code{VectorFromSkewMatrix(skew.mat)}}{In 3D, compute the vector defining the cross product
      associated to the skew-symmetric matrix skew mat.

      In nD, fill a vector by reading the values
      of the upper triangle of skew_mat.}
}}


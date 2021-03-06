% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AR.data}
\alias{AR.data}
\title{Skeleton data for 6 different actions.}
\format{
A list of 6 different actions:
\describe{
\item{come}{Gesture for telling the robot to come to you. 46 instances of [50x92] matrices.}
\item{five}{Gesture of 'high five'. 45 instances of [50x92] matrices.}
\item{handshake}{Gesture of handshaking with the robot. 45 instances of [50x92] matrices.}
\item{hello}{Gesture for telling hello to the robot. 44 instances of [50x92] matrices.}
\item{ignore}{Ignore the robot, pass by. 46 instances of [50x92] matrices.}
\item{look_at}{Stare at the robot in front of it. 46 instances of [50x92] matrices.}
}
}
\usage{
AR.data
}
\description{
A dataset containing the skeleton data extracted from videos of people performing
six different actions, recorded by a humanoid robot. Each class has several instances with
50 values (25 XY body keypoints extracted using OpenPose \url{https://github.com/CMU-Perceptual-Computing-Lab/openpose})
per frame of the video (92 frames).
}
\examples{
X1 <- AR.data$come
X2 <- AR.data$five
mydbcsp <- new('dbcsp',X1,X2)
}
\keyword{datasets}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROI.R
\name{drawROIOnImage}
\alias{drawROIOnImage}
\title{Select a single region-of-interest (ROI) on optical image}
\usage{
drawROIOnImage(img)
}
\arguments{
\item{img}{An \code{EBImage::Image} object to select ROI on}
}
\value{
Returns the ROI as a logical matrix of the same spatial
dimensions as \code{img}
}
\description{
Select a single region-of-interest (ROI) on optical image
}

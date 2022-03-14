#' Landsat time series of spectral indices.
#'
#' A raster dataset containing annual spectral indices from 1985 to 2020.
#' Spectral indices are:
#' the Moisture Stress Index \insertCite{@MSI; @hunt1989detection}{cdmt};
#' the Tasseled cap Wetness \insertCite{@TCW; @crist1985tm}{cdmt};
#' the Tasseled cap Angle \insertCite{@TCA; @powell2010quantification}{cdmt}.
#' The spatial subset covers a wildfire occurred in the Aosta Valley (Italy) during the winter of 2005.
#' Landsat imagery was acquired by the Thematic Mapper (TM), Enhanced Thematic Mapper + (ETM+) and Operational Land Imager (OLI)
#' sensors. Surface reflectance was derived from Level 1 Tier 1 imagery from Landsat Collection 1 in the USGS archive using the
#' Framework for Operational Radiometric Correction for Environmental monitoring (FORCE) software \insertCite{frantz2019force}{cdmt}.
#' Annual reflectance composites were produced using all the images acquired during the growing season (between June 1 and September 30).
#' The geometric median approach \insertCite{roberts2017high}{cdmt} was employed to generate reflectance composites at the pixel level.
#'
#' @format A \code{PackedSpatRaster} object with 108 layers, three per year.
#'
#' @references
#' \insertRef{crist1985tm}{cdmt}
#'
#' \insertRef{hunt1989detection}{cdmt}
#'
#' \insertRef{powell2010quantification}{cdmt}
"lnd_si"

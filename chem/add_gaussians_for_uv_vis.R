#!/usr/bin/env Rscript

# UV-Vis/Fluorescence
# Note that functions are expecting output from autochem, in the form
# Config,Root,Iteration,Transition Energies (eV),Wavelength (nm),Intensity (au)

library(readr)
library(dplyr)
library(magrittr)

normpdf <- function(x, mu, sigma) {
  u <- (x - mu) / abs(sigma)
  y <- (1 / (sqrt(2 * pi) * abs(sigma))) * exp(-u**2 / 2)
  return(y)
}

smooth_with_gaussians <- function(df, sigma, step) {
  if (missing(sigma)) {
    sigma <- 0.05
  }
  if (missing(step)) {
    step <- 0.01
  }

  h <- 6.62607004e-34
  c <- 299792458
  e <- 1.6021766208e-19

  min_transition_energy <- min(df$`Transition Energies (eV)`) - (5 * sigma)
  max_transition_energy <- max(df$`Transition Energies (eV)`) + (5 * sigma)

  # create new energy scale based on transition energies
  # (makes sense; no transition = no absorbance = no peak)
  new_energies <- seq(min_transition_energy, max_transition_energy, step)
  # convert to wavelengths (eV -> nm)
  new_waves <- (h * c * 1e9 / (new_energies * e))

  # take number of x values, set to 0, then move along the
  # line and fill values of intensity as you go, using the index of the
  # intensity list to find the oscillator strength.

  new_ints <- rep(0, length(new_waves))

  for (i in 1:length(df$`Intensity (au)`)) { # want the index
    new_ints <- new_ints + df$`Intensity (au)`[i] * normpdf(
      new_energies,
      df$`Transition Energies (eV)`[i],
      sigma
    )
  }
  return(list(new_waves = new_waves, new_ints = new_ints))
}

add_gaussians <- function(original_df, sigma, step) {
  new <- smooth_with_gaussians(original_df, sigma, step)
  ret <- data.frame(
    new_waves = new$new_waves,
    new_ints = new$new_ints
  )
  # extend original df to make it the same length as the new one
  rows_to_add <- nrow(ret) - nrow(original_df)
  ret$raw_waves <- c(original_df$`Wavelength (nm)`, rep(NA, rows_to_add))
  ret$raw_ints <- c(original_df$`Intensity (au)`, rep(NA, rows_to_add))
  return(ret)
}

read_csv('uv_vis.csv') %>% group_by(Config) %>% do(add_gaussians(.)) %>% write_csv('uv_vis_spectra.csv')

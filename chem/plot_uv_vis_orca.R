#!/usr/bin/env Rscript

library(tidyverse)

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

  min_transition_energy <- min(df$Transitions) - (5 * sigma)
  max_transition_energy <- max(df$Transitions) + (5 * sigma)

  # create new energy scale based on transition energies
  # (makes sense; no transition = no absorbance = no peak)
  new_energies <- seq(min_transition_energy, max_transition_energy, step)
  # convert to wavelengths (eV -> nm)
  new_waves <- (h * c * 1e9 / (new_energies * e))

  # take number of x values, set to 0, then move along the
  # line and fill values of intensity as you go, using the index of the
  # intensity list to find the oscillator strength.

  new_ints <- rep(0, length(new_waves))

  for (i in 1:length(df$Intensity)) { # want the index
    new_ints <- new_ints + df$Intensity[i] * normpdf(
      new_energies,
      df$Transitions[i],
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
  ret$raw_waves <- c(original_df$Wavelength, rep(NA, rows_to_add))
  ret$raw_ints <- c(original_df$Intensity, rep(NA, rows_to_add))
  return(ret)
}

plot_gaussians <- function(df) {
  return(
    ggplot(df) +
      theme_bw() +
      theme(
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        text = element_text(size = 14),
        axis.title = element_text(size = 12),
        legend.text.align = 0,
        panel.grid = element_blank()
      ) +
      geom_line(aes(new_waves, new_ints), color = "red") +
      geom_segment(aes(
        x = raw_waves, xend = raw_waves,
        y = 0, yend = raw_ints
      ), color = "blue") +
      labs(x = "Wavelength (nm)", y = "Intensity (au)")
  )
}

# Extract data from file
args <- commandArgs(trailingOnly = TRUE)
f <- readLines(args[1])

# loop through file, search for ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
data <- data.frame()
found <- FALSE
for (i in 1:length(f)) {
  if (grepl("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS", f[i])) {
    found <- TRUE
    next
  }
  if (found && f[i] == "") {
    break
  }
  # find states and extract data
  if (found && grepl("^\\s*[0-9]", f[i])) {
    line <- f[i]
    # split the line into columns
    cols <- strsplit(str_trim(line), "\\s+")
    cols <- cols[[1]] # R nests lists (almost like an index and list I think)
    # check if state can be numeric, otherwise we might take the wrong lines
    wavelength <- as.numeric(cols[3])
    intensity <- as.numeric(cols[4]) # fosc
    transition_energy <- as.numeric(cols[2])
    # add the data to the data frame
    data <- rbind(data, list(Wavelength = wavelength, Intensity = intensity, Transitions = transition_energy))
  }
}

inverse_cm_to_eV <- 1 / 8065.6

data %>%
  arrange(Wavelength) %>%
  mutate(Transitions = Transitions * inverse_cm_to_eV) %>%
  add_gaussians() %>%
  plot_gaussians() +
  ggsave("uv_vis.png", dpi = 300, width = 6, height = 5)

# LICOR CO2 peak extraction pipeline - generalized maxima-first version v5 human-review
# -------------------------------------------------------------------------
# Goal:
#   1) Read LICOR .txt output.
#   2) Detect injected-sample CO2 peaks without hard-coding the expected count.
#   3) Extract peak maximum CO2 and instantaneous H2O at that maximum.
#   4) Define peak boundaries from valleys between detected maxima.
#   5) Calculate baseline-corrected CO2 AUC and average H2O across each peak.
#   6) Export polished CSV plus diagnostic files.
#
# Main v5 change:
#   Adds an optional human-review step. The script first creates diagnostic files
#   with provisional peak numbers. After inspecting the PDF, list any false
#   provisional peak numbers in manual_exclude_peak_numbers and rerun. The final
#   CSVs are then regenerated after removing those user-flagged false positives.
#   This keeps the algorithm generalized while giving you a transparent QC step.

# --------------------------- USER SETTINGS ------------------------------
input_file <- "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/PennState/AssistantProfessor/BLUE_Shared_Folder/Methods/LiCOR_Respiration/LiCOR_DataPipeline_Repo/LiCOR_Project/Input/Day1_T0_LICOR_GW_SP.txt"
output_prefix <- "~/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/PennState/AssistantProfessor/BLUE_Shared_Folder/Methods/LiCOR_Respiration/LiCOR_DataPipeline_Repo/LiCOR_Project/Output/Day1_T0_LICOR"

# Use this only for checking a validation file. Set to NULL for unknown future files.
expected_peak_count_for_validation <- 46
validation_tolerance <- 1

# Peak detection parameters. These are intentionally general, not count-specific.
smooth_window_points <- 3        # small median smoothing; 3 preserves narrow injection spikes
baseline_window_points <- 121    # rolling low-quantile baseline window; make odd if possible
baseline_quantile <- 0.20        # lower quantile approximates local baseline under peaks
noise_multiplier <- 3.0          # adaptive threshold = noise_multiplier * MAD noise
min_abs_prominence <- 6.25       # floor for real peaks on detrended CO2; tune down only if true small peaks are missed
min_peak_gap_sec <- 3.0          # minimum allowed spacing between peak maxima
prominence_window_sec <- 25      # local window for calculating peak prominence

# Boundary/AUC settings
edge_fraction <- 0.05            # optional edge crossing threshold relative to peak height
use_edge_refinement <- FALSE     # FALSE uses valleys between maxima; TRUE tightens to edge_fraction crossing
baseline_correct_area <- TRUE    # recommended TRUE for injected peak area

# Shape-pruning settings. These are not count-specific. They remove shoulders that
# pass prominence filtering but do not have a real rise-and-fall shape.
require_two_sided_peak <- TRUE
min_flank_width_sec <- 0.5       # each side of a peak must span at least this much time

### Human review / diagnostics
### Set to TRUE for the first pass. Inspect the diagnostic PDF. Then put false
### provisional peak numbers into manual_exclude_peak_numbers and rerun.
make_diagnostic_pdf <- TRUE
manual_exclude_peak_numbers <- c(22, 24)  # e.g., c(22, 24); use integer labels from diagnostic PDF; NULL if you don't want manual exclusion. 
manual_exclude_peak_times_sec <- NULL     # optional alternative, e.g., c(385.2, 612.7)
time_match_tolerance_sec <- 1.5           # used only for manual_exclude_peak_times_sec

# If TRUE, write an extra provisional table before manual exclusion. This is
# useful because manual_exclude_peak_numbers refer to provisional labels.
write_provisional_outputs <- TRUE

# --------------------------- HELPER FUNCTIONS ---------------------------
make_odd <- function(x) {
  x <- as.integer(round(x))
  if (x < 1) x <- 1
  if (x %% 2 == 0) x <- x + 1
  x
}

rolling_stat <- function(x, k, FUN, ...) {
  k <- make_odd(k)
  n <- length(x)
  half <- floor(k / 2)
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    out[i] <- FUN(x[lo:hi], na.rm = TRUE, ...)
  }
  out
}

rolling_median <- function(x, k) {
  rolling_stat(x, k, stats::median)
}

rolling_quantile <- function(x, k, probs = 0.2) {
  rolling_stat(x, k, stats::quantile, probs = probs, names = FALSE, type = 7)
}

trapz <- function(x, y) {
  if (length(x) < 2) return(0)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

# Make a strictly increasing elapsed time vector even when LICOR repeats whole-second timestamps.
make_elapsed_seconds <- function(date_vec, time_vec) {
  timestamp <- as.POSIXct(paste(date_vec, time_vec), tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  raw_sec <- as.numeric(difftime(timestamp, timestamp[1], units = "secs"))

  elapsed <- raw_sec
  groups <- split(seq_along(raw_sec), raw_sec)
  unique_times <- sort(unique(raw_sec))

  for (g in groups) {
    n_g <- length(g)
    if (n_g > 1) {
      current_time <- raw_sec[g[1]]
      later_times <- unique_times[unique_times > current_time]
      if (length(later_times) > 0) {
        span <- min(later_times) - current_time
      } else {
        positive_steps <- diff(unique_times)
        span <- ifelse(length(positive_steps) > 0, stats::median(positive_steps), 1)
      }
      elapsed[g] <- current_time + seq(0, span * (n_g - 1) / n_g, length.out = n_g)
    }
  }

  # Final guard against any remaining ties/non-increasing values.
  for (i in 2:length(elapsed)) {
    if (!is.finite(elapsed[i]) || elapsed[i] <= elapsed[i - 1]) {
      elapsed[i] <- elapsed[i - 1] + 1e-6
    }
  }
  elapsed
}

find_local_maxima <- function(y) {
  n <- length(y)
  peaks <- integer(0)
  i <- 2
  while (i <= n - 1) {
    if (is.na(y[i])) {
      i <- i + 1
      next
    }

    # Plateau-aware local maximum detection.
    if (y[i] >= y[i - 1] && y[i] >= y[i + 1] && (y[i] > y[i - 1] || y[i] > y[i + 1])) {
      start <- i
      end <- i
      while (end < n && y[end + 1] == y[start]) end <- end + 1
      left_ok <- start == 1 || y[start] > y[start - 1]
      right_ok <- end == n || y[end] > y[end + 1]
      if (left_ok && right_ok) peaks <- c(peaks, round((start + end) / 2))
      i <- end + 1
    } else {
      i <- i + 1
    }
  }
  unique(peaks)
}

calculate_prominence <- function(y, peak_idx, elapsed_sec, window_sec) {
  p <- peak_idx
  t <- elapsed_sec[p]
  lo <- max(1, min(which(elapsed_sec >= t - window_sec)))
  hi <- min(length(y), max(which(elapsed_sec <= t + window_sec)))

  left_min <- if (lo < p) min(y[lo:p], na.rm = TRUE) else y[p]
  right_min <- if (p < hi) min(y[p:hi], na.rm = TRUE) else y[p]

  # Match standard peak-prominence logic: peak minus the higher of the two bases.
  y[p] - max(left_min, right_min)
}

apply_min_gap <- function(candidates, elapsed_sec, min_gap_sec) {
  if (nrow(candidates) == 0) return(candidates)

  # Keep strongest candidates first, then restore chronological order.
  candidates <- candidates[order(-candidates$prominence), , drop = FALSE]
  kept <- candidates[0, , drop = FALSE]

  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i, , drop = FALSE]
    if (nrow(kept) == 0) {
      kept <- rbind(kept, cand)
    } else {
      dt <- abs(elapsed_sec[cand$peak_index] - elapsed_sec[kept$peak_index])
      if (all(dt >= min_gap_sec)) kept <- rbind(kept, cand)
    }
  }

  kept[order(kept$peak_index), , drop = FALSE]
}


parse_manual_numbers <- function(x) {
  if (is.null(x) || length(x) == 0) return(integer(0))
  x <- as.integer(x)
  x[is.finite(x) & !is.na(x) & x > 0]
}

find_valley_between <- function(y, left_peak, right_peak) {
  if (right_peak <= left_peak + 1) return(left_peak)
  segment <- y[left_peak:right_peak]
  left_peak + which.min(segment) - 1
}

find_outer_valley <- function(y, peak, elapsed_sec, direction = c("left", "right"), window_sec = 25) {
  direction <- match.arg(direction)
  t <- elapsed_sec[peak]
  if (direction == "left") {
    idx <- which(elapsed_sec >= t - window_sec & elapsed_sec <= t)
  } else {
    idx <- which(elapsed_sec >= t & elapsed_sec <= t + window_sec)
  }
  if (length(idx) == 0) return(peak)
  idx[which.min(y[idx])]
}

refine_edge <- function(y, baseline_y, peak_idx, boundary_idx, direction = c("left", "right"), edge_fraction = 0.05) {
  direction <- match.arg(direction)
  peak_height <- y[peak_idx] - baseline_y[peak_idx]
  cutoff <- baseline_y[peak_idx] + edge_fraction * peak_height

  if (direction == "left") {
    search_idx <- seq(peak_idx, boundary_idx, by = -1)
  } else {
    search_idx <- seq(peak_idx, boundary_idx, by = 1)
  }

  crossed <- search_idx[y[search_idx] <= cutoff]
  if (length(crossed) == 0) return(boundary_idx)
  crossed[1]
}

# --------------------------- READ DATA ----------------------------------
if (!file.exists(input_file)) {
  stop("Could not find input_file. Set input_file to the LICOR .txt path.")
}

raw <- read.delim(
  input_file,
  skip = 1,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw <- raw[, !grepl("^Unnamed|^$", names(raw)), drop = FALSE]

co2_col <- "CO₂_(µmol_mol⁻¹)"
h2o_col <- "H₂O_(mmol_mol⁻¹)"
date_col <- "System_Date_(Y-M-D)"
time_col <- "System_Time_(h:m:s)"

needed <- c(date_col, time_col, co2_col, h2o_col)
missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

co2 <- as.numeric(raw[[co2_col]])
h2o <- as.numeric(raw[[h2o_col]])
elapsed_sec <- make_elapsed_seconds(raw[[date_col]], raw[[time_col]])

valid <- is.finite(co2) & is.finite(h2o) & is.finite(elapsed_sec)
raw <- raw[valid, , drop = FALSE]
co2 <- co2[valid]
h2o <- h2o[valid]
elapsed_sec <- elapsed_sec[valid]

# --------------------------- DETECT MAXIMA -------------------------------
co2_smooth <- rolling_median(co2, smooth_window_points)
baseline <- rolling_quantile(co2_smooth, baseline_window_points, baseline_quantile)
detrended <- co2_smooth - baseline

noise_mad <- stats::mad(diff(detrended), constant = 1.4826, na.rm = TRUE)
adaptive_min_prominence <- max(min_abs_prominence, noise_multiplier * noise_mad)

all_maxima <- find_local_maxima(detrended)
if (length(all_maxima) == 0) stop("No local maxima found. Check the CO2 column and input file.")

prom <- vapply(
  all_maxima,
  calculate_prominence,
  numeric(1),
  y = detrended,
  elapsed_sec = elapsed_sec,
  window_sec = prominence_window_sec
)

candidate_table <- data.frame(
  peak_index = all_maxima,
  time_sec = elapsed_sec[all_maxima],
  co2_raw = co2[all_maxima],
  co2_smooth = co2_smooth[all_maxima],
  detrended_height = detrended[all_maxima],
  prominence = prom
)

candidates <- candidate_table[candidate_table$prominence >= adaptive_min_prominence, , drop = FALSE]
peaks <- apply_min_gap(candidates, elapsed_sec, min_peak_gap_sec)

if (nrow(peaks) == 0) {
  stop("No peaks passed the adaptive prominence threshold. Try lowering min_abs_prominence or noise_multiplier.")
}

# --------------------------- BOUNDARIES + AUC ----------------------------
peak_idx <- peaks$peak_index
n_peaks <- length(peak_idx)
left_idx <- integer(n_peaks)
right_idx <- integer(n_peaks)

between_valleys <- integer(max(0, n_peaks - 1))
if (n_peaks > 1) {
  for (i in seq_len(n_peaks - 1)) {
    between_valleys[i] <- find_valley_between(detrended, peak_idx[i], peak_idx[i + 1])
  }
}

for (i in seq_len(n_peaks)) {
  left_idx[i] <- if (i == 1) {
    find_outer_valley(detrended, peak_idx[i], elapsed_sec, "left", prominence_window_sec)
  } else {
    between_valleys[i - 1]
  }

  right_idx[i] <- if (i == n_peaks) {
    find_outer_valley(detrended, peak_idx[i], elapsed_sec, "right", prominence_window_sec)
  } else {
    between_valleys[i]
  }

  if (use_edge_refinement) {
    left_idx[i] <- refine_edge(co2_smooth, baseline, peak_idx[i], left_idx[i], "left", edge_fraction)
    right_idx[i] <- refine_edge(co2_smooth, baseline, peak_idx[i], right_idx[i], "right", edge_fraction)
  }

  if (left_idx[i] > peak_idx[i]) left_idx[i] <- peak_idx[i]
  if (right_idx[i] < peak_idx[i]) right_idx[i] <- peak_idx[i]
}

# Prune one-sided shoulder artifacts after boundaries are calculated.
# These are candidates where the valley between adjacent maxima lands exactly on
# the candidate maximum, so the peak has no true left or right flank. This is a
# general shape rule, not an expected-count rule.
left_flank_sec <- elapsed_sec[peak_idx] - elapsed_sec[left_idx]
right_flank_sec <- elapsed_sec[right_idx] - elapsed_sec[peak_idx]
keep_shape <- rep(TRUE, n_peaks)
if (require_two_sided_peak) {
  keep_shape <- left_flank_sec >= min_flank_width_sec & right_flank_sec >= min_flank_width_sec
}

if (any(!keep_shape)) {
  removed <- data.frame(
    original_order = which(!keep_shape),
    peak_time_sec = elapsed_sec[peak_idx[!keep_shape]],
    left_flank_sec = left_flank_sec[!keep_shape],
    right_flank_sec = right_flank_sec[!keep_shape],
    prominence = peaks$prominence[!keep_shape]
  )
  message("Shape pruning removed ", nrow(removed), " one-sided candidate peak(s).")
}

peaks <- peaks[keep_shape, , drop = FALSE]
peak_idx <- peak_idx[keep_shape]
left_idx <- left_idx[keep_shape]
right_idx <- right_idx[keep_shape]
n_peaks <- length(peak_idx)

# Provisional labels are the labels shown on the diagnostic PDF before manual QC.
provisional_peak_number <- seq_len(n_peaks)
peaks$provisional_peak_number <- provisional_peak_number

# Optional human-in-the-loop false-positive removal. This does NOT force a target
# count. It simply removes specific provisional labels or times that a reviewer
# marks as false after inspecting the diagnostic PDF.
manual_numbers <- parse_manual_numbers(manual_exclude_peak_numbers)
manual_keep <- rep(TRUE, n_peaks)

if (write_provisional_outputs) {
  provisional_file <- paste0(output_prefix, "_CO2_provisional_peaks_for_review.csv")
  provisional_table <- data.frame(
    provisional_peak_number = provisional_peak_number,
    left_index = left_idx,
    peak_index = peak_idx,
    right_index = right_idx,
    left_time_sec = elapsed_sec[left_idx],
    peak_time_sec = elapsed_sec[peak_idx],
    right_time_sec = elapsed_sec[right_idx],
    co2_raw = co2[peak_idx],
    h2o_at_peak = h2o[peak_idx],
    prominence = peaks$prominence,
    detrended_height = peaks$detrended_height
  )
  write.csv(provisional_table, provisional_file, row.names = FALSE, fileEncoding = "UTF-8")
  message("Wrote provisional review table: ", provisional_file)
}

if (length(manual_numbers) > 0) {
  invalid_numbers <- setdiff(manual_numbers, provisional_peak_number)
  if (length(invalid_numbers) > 0) {
    warning("Manual exclude number(s) not found among provisional peaks: ",
            paste(invalid_numbers, collapse = ", "))
  }
  manual_keep <- manual_keep & !(provisional_peak_number %in% manual_numbers)
}

if (!is.null(manual_exclude_peak_times_sec) && length(manual_exclude_peak_times_sec) > 0) {
  for (t_ex in manual_exclude_peak_times_sec) {
    nearest <- which.min(abs(elapsed_sec[peak_idx] - t_ex))
    if (length(nearest) == 1 && abs(elapsed_sec[peak_idx[nearest]] - t_ex) <= time_match_tolerance_sec) {
      manual_keep[nearest] <- FALSE
    } else {
      warning("No provisional peak found within ", time_match_tolerance_sec,
              " s of manual exclude time: ", t_ex)
    }
  }
}

if (any(!manual_keep)) {
  removed_manual <- data.frame(
    provisional_peak_number = provisional_peak_number[!manual_keep],
    peak_time_sec = elapsed_sec[peak_idx[!manual_keep]],
    co2_raw = co2[peak_idx[!manual_keep]],
    prominence = peaks$prominence[!manual_keep]
  )
  message("Manual review removed ", nrow(removed_manual), " provisional peak(s): ",
          paste(removed_manual$provisional_peak_number, collapse = ", "))
}

peaks <- peaks[manual_keep, , drop = FALSE]
peak_idx <- peak_idx[manual_keep]
left_idx <- left_idx[manual_keep]
right_idx <- right_idx[manual_keep]
provisional_peak_number <- provisional_peak_number[manual_keep]
n_peaks <- length(peak_idx)

summary_rows <- vector("list", n_peaks)
boundary_rows <- vector("list", n_peaks)

for (i in seq_len(n_peaks)) {
  idx <- left_idx[i]:right_idx[i]
  p <- peak_idx[i]

  # Use raw CO2 maximum within the boundary, not necessarily smoothed/detrended max.
  raw_peak_idx <- idx[which.max(co2[idx])]

  area_y <- if (baseline_correct_area) {
    pmax(co2[idx] - baseline[idx], 0)
  } else {
    co2[idx]
  }

  summary_rows[[i]] <- data.frame(
    sample_number = i,
    `maximum co2 (µmol_mol⁻¹)` = co2[raw_peak_idx],
    `Instantaneous H2O` = h2o[raw_peak_idx],
    `CO2 curve area` = trapz(elapsed_sec[idx], area_y),
    `Average H2O` = mean(h2o[idx], na.rm = TRUE),
    check.names = FALSE
  )

  boundary_rows[[i]] <- data.frame(
    sample_number = i,
    provisional_peak_number = provisional_peak_number[i],
    left_index = left_idx[i],
    peak_index = raw_peak_idx,
    right_index = right_idx[i],
    left_time_sec = elapsed_sec[left_idx[i]],
    peak_time_sec = elapsed_sec[raw_peak_idx],
    right_time_sec = elapsed_sec[right_idx[i]],
    detected_prominence = peaks$prominence[i],
    detected_detrended_height = peaks$detrended_height[i]
  )
}

peak_summary <- do.call(rbind, summary_rows)
peak_boundaries <- do.call(rbind, boundary_rows)

# --------------------------- EXPORTS -------------------------------------
summary_file <- paste0(output_prefix, "_CO2_peak_summary.csv")
boundaries_file <- paste0(output_prefix, "_CO2_peak_summary_boundaries.csv")
candidates_file <- paste0(output_prefix, "_CO2_peak_candidates.csv")

write.csv(peak_summary, summary_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(peak_boundaries, boundaries_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(candidate_table[order(-candidate_table$prominence), ], candidates_file, row.names = FALSE, fileEncoding = "UTF-8")
if (exists("removed_manual")) {
  manual_removed_file <- paste0(output_prefix, "_CO2_manual_removed_peaks.csv")
  write.csv(removed_manual, manual_removed_file, row.names = FALSE, fileEncoding = "UTF-8")
}

message("Detected final peaks after automatic + manual QC: ", n_peaks)
if (length(manual_numbers) > 0 || (!is.null(manual_exclude_peak_times_sec) && length(manual_exclude_peak_times_sec) > 0)) {
  message("Manual exclusions were applied. Provisional labels are retained in the boundaries CSV.")
}
message("Adaptive minimum prominence used: ", round(adaptive_min_prominence, 4))
message("Noise MAD from detrended first differences: ", round(noise_mad, 4))
message("Wrote: ", summary_file)
message("Wrote: ", boundaries_file)
message("Wrote: ", candidates_file)
if (exists("manual_removed_file")) message("Wrote: ", manual_removed_file)

if (!is.null(expected_peak_count_for_validation)) {
  difference <- n_peaks - expected_peak_count_for_validation
  if (abs(difference) > validation_tolerance) {
    warning(
      "Validation warning: detected ", n_peaks,
      " peaks, expected about ", expected_peak_count_for_validation,
      ". This is only a warning; no peaks were added or removed to force this count."
    )
  } else {
    message("Validation check passed: detected count is within tolerance of expected count.")
  }
}

# --------------------------- DIAGNOSTIC PLOT -----------------------------
if (make_diagnostic_pdf) {
  pdf_file <- paste0(output_prefix, "_CO2_peak_diagnostic.pdf")
  grDevices::pdf(pdf_file, width = 12, height = 8)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

  plot(elapsed_sec, co2, type = "l", xlab = "Elapsed time (s)", ylab = co2_col,
       main = "Raw CO2 with detected peak maxima and boundaries")
  points(elapsed_sec[peak_boundaries$peak_index], co2[peak_boundaries$peak_index], pch = 19)
  text(elapsed_sec[peak_boundaries$peak_index], co2[peak_boundaries$peak_index],
       labels = peak_boundaries$sample_number, pos = 3, cex = 0.7)
  abline(v = elapsed_sec[peak_boundaries$left_index], lty = 3)
  abline(v = elapsed_sec[peak_boundaries$right_index], lty = 3)

  plot(elapsed_sec, detrended, type = "l", xlab = "Elapsed time (s)", ylab = "Detrended CO2",
       main = "Detrended CO2 used for maxima-first detection")
  abline(h = adaptive_min_prominence, lty = 2)
  points(elapsed_sec[peaks$peak_index], detrended[peaks$peak_index], pch = 19)
  text(elapsed_sec[peaks$peak_index], detrended[peaks$peak_index],
       labels = seq_len(nrow(peaks)), pos = 3, cex = 0.7)

  grDevices::dev.off()
  message("Wrote: ", pdf_file)
}

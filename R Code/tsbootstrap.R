// R function 
  
  
tsbootstrap <- function(data, ...) 
{# begin-of-function
  UseMethod("tsbootstrap")
}

TSBOOT_TYPES <- c("fixed", "geom")

#' @describeIn tsbootstrap Generate time series replicates using blocks of
#'   rows in the data frame.
#' @export

tsbootstrap.data.frame <- function(data, R = 1L, size = 1L, m = nrow(data), type = "fixed", ...) 
{# begin-of-function
  
  assert_that(is.number(R) && R >= 1)
  assert_that(is.number(size) && size >= 1L && size <= nrow(data))
  assert_that(type %in% TSBOOT_TYPES)
  to_resample_df(tsbootstrap_(nrow(data), R = R, size = size, type = type), data)

}# end-of-function

#' @describeIn tsbootstrap Generate time series replicates using blocks of
#'   groups in the grouped data frame.
#' @export

tsbootstrap.grouped_df <- function(data, R = 1L, size = 1L, type = "fixed", m = dplyr::n_groups(data), ...)
{# begin-of-function
  
  assert_that(is.number(R) && R >= 1)
  assert_that(is.number(size) && size >= 1L && size <= nrow(data))
  assert_that(type %in% TSBOOT_TYPES)
  idx <- group_indices_lst(data)
  f <- function(i) flatten_int(idx[i])  # nolint
  res <- mutate_(tsbootstrap_(length(idx), R = R, size = size, type = type),
                 sample = ~ map(sample, f))
  to_resample_df(res, data)

}# end-of-function

# Adapted from boot:::make.ends
mod <- function(i, n, endcorr) {
  if (endcorr) 1 + (i - 1) %% n
  else i
}

# Fixed Moving Block Bootstrap from boot::tsboot and boot:::ts.array
# Fixed Moving Block Bootstrap from boot::tsboot

.tsboot_mbb <- function(n, m = n, size = 1L, endcorr = TRUE)
{# begin-of-function
  
  endpt <- if (endcorr) n else n - size + 1L
  nn <- ceiling(m / size)
  lens <- c(rep(size, nn - 1L), 1L + (m - 1L) %% size)
  st <- sample.int(endpt, nn, replace = TRUE)
 
  flatten_int(purrr::map2(st, lens, function(s, l)
  {
    if (l >= 1) as.integer(mod(seq(s, s + l - 1L), n, endcorr))
    else integer()
  }))

}# end-of-function

# Fixed Moving Block Bootstrap from boot::tsboot and boot:::ts.array
# adapted from
# use ... to ignore endcorr

.tsboot_geom <- function(n, m = n, size = 1L, ...)
{# begin-of-function
  
  endpt <- n - size + 1L
  # worst case scenario is to draw m. So take m draws from
  # geom and truncate
  lens <- 1L + stats::rgeom(m, 1L / size)
  len_tot <- cumsum(lens)
  # truncate to minimum length >= m
  lens <- lens[seq_len(purrr::detect_index(len_tot, ~ .x >= m))]
  st <- sample.int(endpt, length(lens), replace = TRUE)
  
  f <- function(s, l) 
  {
    if (l >= 1) as.integer(mod(seq.int(s, s + l - 1L), n, TRUE))
    else integer()
  }
  utils::head(flatten_int(purrr::map2(st, lens, f)), m)

}# end-of-function


tsbootstrap_ <- function(n, R = 1L, m = n, size = 1L, type = "fixed", endcorr = TRUE) 
{# begin-of-function
  
  f <- switch(type,
              geom = .tsboot_geom,
              fixed = .tsboot_mbb,
              stop("type = ", type, " is not recognized.", call. = FALSE))
  
  tibble(sample = purrr::rerun(R, f(n, m = m, size = size, endcorr = endcorr)), .id = seq_len(R))

}# end-of-function

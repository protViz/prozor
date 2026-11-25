# Stops unless `value` is a single number > 0 (if `positive`) or >= 0, and a whole number if `whole`; `name` is the
# argument name for the message
.check_number <- function(value, name, positive, whole) {
  scalar <- is.numeric(value) && length(value) == 1 && !is.na(value)
  valid <- scalar && all(value >= 0, value > 0 || !positive, value == round(value) || !whole)
  if (!valid) {
    kind <- paste0(if (positive) "positive" else "non-negative", if (whole) " whole" else "")
    stop("`", name, "` must be a single ", kind, " number, not ", deparse1(value), ".", call. = FALSE)
  }
}

# Stops unless `masses` is a non-empty numeric vector of finite values
.check_masses <- function(masses) {
  if (!is.numeric(masses) || length(masses) == 0) {
    stop(
      "`masses` must be a non-empty numeric vector, not ",
      class(masses)[1],
      " of length ",
      length(masses),
      ".",
      call. = FALSE
    )
  }
  if (!all(is.finite(masses))) {
    stop("`masses` must not contain NA, NaN or infinite values; found ", sum(!is.finite(masses)), ".", call. = FALSE)
  }
}

# The `method` argument: the first of `choices` for the default (all choices), otherwise exactly one of them
.match_method <- function(method, choices) {
  if (identical(method, choices)) {
    return(choices[1])
  }
  if (!(length(method) == 1 && method %in% choices)) {
    stop(
      "`method` must be one of ",
      paste0("\"", choices, "\"", collapse = ", "),
      ", not ",
      deparse1(method),
      ".",
      call. = FALSE
    )
  }
  method
}

# Stops for an assignment to the read-only field `name` of an object of the R6 class `class`
.stop_read_only <- function(name, class) {
  stop("`", name, "` of a ", class, " object is read-only; create a new object with ", class, "$new().", call. = FALSE)
}

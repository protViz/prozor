# Stops unless `value` is a single finite number > 0 (if `positive`) or >= 0, and a whole number if `whole`; `name` is
# the argument name for the message
.check_number <- function(value, name, positive, whole) {
  scalar <- is.numeric(value) && length(value) == 1 && is.finite(value)
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

# Stops unless the data frame or list `x` has the columns `needed`; `arg` is the argument name for the message
.check_columns <- function(x, needed, arg) {
  missing_columns <- setdiff(needed, names(x))
  if (length(missing_columns) > 0) {
    stop(
      "`",
      arg,
      "` must have the columns ",
      .quoted(needed, "`"),
      "; missing: ",
      .quoted(missing_columns, "`"),
      ".",
      call. = FALSE
    )
  }
}

# The argument `method`, named `arg` in the message: the first of `choices` for the default (all choices), otherwise
# exactly one of them
.match_method <- function(method, choices, arg) {
  if (identical(method, choices)) {
    return(choices[1])
  }
  if (!(length(method) == 1 && method %in% choices)) {
    stop(
      "`",
      arg,
      "` must be one of ",
      .quoted(choices, "\""),
      ", not ",
      deparse1(method),
      ".",
      call. = FALSE
    )
  }
  method
}

# Body of the active binding of the read-only field `name` of an object of the R6 class `class`: `current`, or an error
# if a `value` is assigned
.read_only <- function(value, name, class, current) {
  if (!missing(value)) {
    stop(
      "`",
      name,
      "` of a ",
      class,
      " object is read-only; create a new object with ",
      class,
      "$new().",
      call. = FALSE
    )
  }
  current
}

# The elements of `x` in the quote character `quote`, escaped as by encodeString() (NA without quotes), separated by
# commas
.quoted <- function(x, quote) {
  paste(encodeString(as.character(x), quote = quote), collapse = ", ")
}

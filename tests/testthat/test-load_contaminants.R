test_that("load_contaminants() returns each set and no_human removes human entries", {
  # entries kept with no_human = TRUE
  kept <- c(fgcz2024 = 247, fgcz2019 = 340, fgcz2021 = 343, fgcz2022 = 349, universal2024 = 231, special2024 = 16)
  for (set in names(kept)) {
    all_entries <- load_contaminants(set)
    no_human <- load_contaminants(set, no_human = TRUE)
    annot <- vapply(no_human, attr, character(1), "Annot")
    expect_gt(length(all_entries), 0, label = paste("entries of", set))
    expect_lt(length(no_human), length(all_entries), label = paste("entries of", set, "(no_human = TRUE)"))
    expect_length(no_human, kept[[set]])
    expect_false(any(grepl("HUMAN", annot)), label = paste(set, "(no_human = TRUE) has HUMAN entries"))
    # e.g. "Type I hair keratin KA36 - Homo sapiens (Human)." in the fgcz sets
    expect_false(any(grepl("Homo sapiens", annot)), label = paste(set, "(no_human = TRUE) has Homo sapiens entries"))
  }
})

test_that("load_contaminants() loads fgcz2024 by default: special2024, then universal2024", {
  default <- load_contaminants()

  expect_identical(default, load_contaminants("fgcz2024"))
  expect_identical(default, c(load_contaminants("special2024"), load_contaminants("universal2024")))
})

test_that("load_contaminants() rejects unknown sets", {
  expect_error(load_contaminants("fgcz2023"), "`set` must be one of .* not \"fgcz2023\"")
})

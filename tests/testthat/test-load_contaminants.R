test_that("load_contaminants() returns each set and no_human removes human entries", {
  for (set in c("fgcz2019", "fgcz2021", "fgcz2022", "universal2024", "special2024")) {
    all_entries <- load_contaminants(set)
    no_human <- load_contaminants(set, no_human = TRUE)
    annot <- vapply(no_human, attr, character(1), "Annot")
    expect_gt(length(all_entries), 0)
    expect_lt(length(no_human), length(all_entries))
    expect_false(any(grepl("HUMAN", annot)), label = paste(set, "(no_human = TRUE) has HUMAN entries"))
    # e.g. "Type I hair keratin KA36 - Homo sapiens (Human)." in the fgcz sets
    expect_false(any(grepl("Homo sapiens", annot)), label = paste(set, "(no_human = TRUE) has Homo sapiens entries"))
  }
})

test_that("load_contaminants() rejects unknown sets", {
  expect_error(load_contaminants("fgcz2023"), "`set` must be one of .* not \"fgcz2023\"")
})

test_that("all bundled FASTA files are present", {
  # 5 contaminant files in inst/extdata, 2 in inst/p1000_db1_example
  expect_length(bundled_fasta_files(), 7)
})

test_that("every bundled FASTA file is a plain gzipped FASTA of amino-acid sequences", {
  for (file in bundled_fasta_files()) {
    # a tar archive starts with a tar header instead of ">" and makes readPeptideFasta() warn
    expect_match(readLines(file, n = 1), "^>", label = paste(basename(file), "first line"))
    expect_no_warning(fasta <- readPeptideFasta(file))
    seqs <- unlist(fasta)
    expect_gt(length(fasta), 0)
    expect_false(anyDuplicated(names(fasta)) > 0, label = paste(basename(file), "has duplicated names"))
    expect_true(all(nchar(seqs) > 0), label = paste(basename(file), "has only non-empty sequences"))
    expect_true(all(grepl("^[A-Za-z*]+$", seqs)), label = paste(basename(file), "has only amino-acid letters"))
  }
})

test_that("contaminant loaders return sequences and noHuman removes human entries", {
  loaders <- list(
    loadContaminantsFasta2019 = loadContaminantsFasta2019,
    loadContaminantsFasta2021 = loadContaminantsFasta2021,
    loadContaminantsFGCZ2022 = loadContaminantsFGCZ2022,
    load_universal_contaminants_2024 = load_universal_contaminants_2024,
    load_special_prot_2024 = load_special_prot_2024
  )
  for (name in names(loaders)) {
    all_entries <- loaders[[name]]()
    no_human <- loaders[[name]](noHuman = TRUE)
    annot <- vapply(no_human, seqinr::getAnnot, character(1))
    expect_gt(length(all_entries), 0)
    expect_lt(length(no_human), length(all_entries))
    expect_false(any(grepl("HUMAN", annot)), label = paste(name, "(noHuman = TRUE) has HUMAN entries"))
  }
})

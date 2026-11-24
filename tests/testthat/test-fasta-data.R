bundled_fasta_files <- function() {
  c(
    list.files(system.file("extdata", package = "prozor"),
               pattern = "\\.fasta\\.gz$", full.names = TRUE),
    list.files(system.file("p1000_db1_example", package = "prozor"),
               pattern = "\\.fasta\\.gz$", full.names = TRUE)
  )
}

test_that("all bundled FASTA files are present", {
  expect_length(bundled_fasta_files(), 7)
})

test_that("every bundled FASTA file reads as amino-acid sequences", {
  for (file in bundled_fasta_files()) {
    fasta <- readPeptideFasta(file)
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

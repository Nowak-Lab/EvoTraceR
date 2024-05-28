get_pwalign_score_func <- function() {
  bioc_version <- BiocManager::version()
  if (bioc_version >= "3.19") {
    return(pwalign::score)
  } else {
    return(Biostrings::score)
  }
}

get_pwalign_func <- function() {
  bioc_version <- BiocManager::version()
  if (bioc_version >= "3.19") {
    return(pwalign::pairwiseAlignment)
  } else {
    return(Biostrings::pairwiseAlignment)
  }
}

get_pwalign_nedit_func <- function() {
  bioc_version <- BiocManager::version()
  if (bioc_version >= "3.19") {
    return(pwalign::nedit)
  } else {
    return(Biostrings::nedit)
  }
}

get_pwalign_pid_func <- function() {
  bioc_version <- BiocManager::version()
  if (bioc_version >= "3.19") {
    return(pwalign::pid)
  } else {
    return(Biostrings::pid)
  }
}

get_pwalign_subject_func <- function() {
  bioc_version <- BiocManager::version()
  if (bioc_version >= "3.19") {
    return(pwalign::alignedSubject)
  } else {
    return(Biostrings::alignedSubject)
  }
}

get_pwalign_pattern_func <- function() {
  bioc_version <- BiocManager::version()
  if (bioc_version >= "3.19") {
    return(pwalign::alignedPattern)
  } else {
    return(Biostrings::alignedPattern)
  }
}

get_pwalign_nt_sub_matrix_func <- function() {
  bioc_version <- BiocManager::version()
  if (bioc_version >= "3.19") {
    return(pwalign::nucleotideSubstitutionMatrix)
  } else {
    return(Biostrings::nucleotideSubstitutionMatrix)
  }
}
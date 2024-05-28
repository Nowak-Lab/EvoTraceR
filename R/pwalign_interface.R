get_pwalign_score_func <- function() {
  if (packageVersion("Biostrings") > "2.70.0") {
    return(pwalign::score)
  } else {
    return(Biostrings::score)
  }
}

get_pwalign_func <- function() {
  if (packageVersion("Biostrings") > "2.70.0") {
    return(pwalign::pairwiseAlignment)
  } else {
    return(Biostrings::pairwiseAlignment)
  }
}

get_pwalign_nedit_func <- function() {
  if (packageVersion("Biostrings") > "2.70.0") {
    return(pwalign::nedit)
  } else {
    return(Biostrings::nedit)
  }
}

get_pwalign_pid_func <- function() {
  if (packageVersion("Biostrings") > "2.70.0") {
    return(pwalign::pid)
  } else {
    return(Biostrings::pid)
  }
}

get_pwalign_subject_func <- function() {
  if (packageVersion("Biostrings") > "2.70.0") {
    return(pwalign::alignedSubject)
  } else {
    return(Biostrings::alignedSubject)
  }
}

get_pwalign_pattern_func <- function() {
  if (packageVersion("Biostrings") > "2.70.0") {
    return(pwalign::alignedPattern)
  } else {
    return(Biostrings::alignedPattern)
  }
}

get_pwalign_nt_sub_matrix_func <- function() {
  if (packageVersion("Biostrings") > "2.70.0") {
    return(pwalign::nucleotideSubstitutionMatrix)
  } else {
    return(Biostrings::nucleotideSubstitutionMatrix)
  }
}
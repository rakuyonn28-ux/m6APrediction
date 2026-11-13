
#' Encode DNA Sequences as Nucleotide-Position Factors
#'
#' Split equal-length DNA sequences into single nucleotides and return one
#' categorical column per nucleotide position. Each position is a factor with
#' levels `"A"`, `"T"`, `"C"` and `"G"`, so the columns can be used as
#' categorical predictors by a random forest model. This is an internal function
#' for [prediction_multiple()].
#'
#' @param dna_strings Character vector of DNA sequences. Each sequence must
#'   contain the bases A, T, C and G, and all sequences must have the
#'   same length.
#'
#' @return A data frame with one row of sequence and one column per nucleotide position,
#'   named `nt_pos1`, `nt_pos2`, ..., `nt_posN`. Every column is a
#'   factor with levels `"A"`, `"T"`, `"C"` and `"G"`.
#'
#' @seealso [prediction_multiple()], [prediction_single()]
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Status for Multiple RNA Sites
#'
#' Apply a fitted random forest classifier to a feature table and classify
#' each site as `"Positive"` or `"Negative"`. The five-base sequence column
#' `DNA_5mer` is expanded into categorical nucleotide-position predictors with
#' [dna_encoding()] before prediction. A site is classified as `"Positive"`
#' only when its predicted probability is strictly greater than
#' `positive_threshold`; a probability equal to the threshold is
#' `"Negative"`.
#'
#' @param ml_fit A fitted random forest classification model, such as an
#'   object returned by [randomForest::randomForest()]. The model must have
#'   been trained on the same seven features and must return a `"Positive"`
#'   class probability when predicted with `type = "prob"`.
#' @param feature_df A data frame with one row per site and the columns
#'   `gc_content` (numeric), `RNA_type`, `RNA_region`, `exon_length`,
#'   `distance_to_junction`, `evolutionary_conservation` and `DNA_5mer`
#'   (character). `RNA_type` must be one of `"mRNA"`, `"lincRNA"`,
#'   `"lncRNA"` or `"pseudogene"`; `RNA_region` must be one of `"CDS"`,
#'   `"intron"`, `"3'UTR"` or `"5'UTR"`. `exon_length` and
#'   `distance_to_junction` must already be log2-transformed. `DNA_5mer`
#'   holds the five-base DNA context with A in the central position.
#' @param positive_threshold Numeric scalar probability threshold between 0
#'   and 1. Defaults to `0.5`.
#'
#' @return The input `feature_df`, with rows in their original order, plus two
#'   columns:
#'   * `predicted_m6A_prob`: numeric predicted probability of the
#'     `"Positive"` class.
#'   * `predicted_m6A_status`: character vector, either `"Positive"` or
#'     `"Negative"`.
#'
#'   `RNA_type` and `RNA_region` are converted to factors with the training
#'   levels before prediction, and those factor columns are returned.
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#'
#' result <- prediction_multiple(ml_fit, example_df[1:5, ],
#'                               positive_threshold = 0.6)
#' result[, c("predicted_m6A_prob", "predicted_m6A_status")]
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict m6A Status for a Single RNA Site
#'
#' Convenience wrapper around [prediction_multiple()] for one RNA site. It
#' builds a one-row feature table, predicts the m6A probability and status,
#' and returns both values in a named character vector.
#'
#' @param ml_fit A fitted random forest classification model, such as an
#'   object returned by [randomForest::randomForest()]. The model must have
#'   been trained on the same seven features and must return a `"Positive"`
#'   class probability when predicted with `type = "prob"`.
#' @param gc_content Numeric GC content of the sequence, between 0 and 1.
#' @param RNA_type Character or factor giving the RNA type. Must be one of
#'   `"mRNA"`, `"lincRNA"`, `"lncRNA"` or `"pseudogene"`.
#' @param RNA_region Character or factor giving the RNA region. Must be one of
#'   `"CDS"`, `"intron"`, `"3'UTR"` or `"5'UTR"`.
#' @param exon_length Numeric exon length. This value must already be
#'   log2-transformed.
#' @param distance_to_junction Numeric distance to the nearest splice junction.
#'   This value must already be log2-transformed.
#' @param evolutionary_conservation Numeric evolutionary conservation score.
#' @param DNA_5mer Character string giving the five-base DNA context with A in
#'   the central position, for example `"GGACA"`.
#' @param positive_threshold Numeric scalar probability threshold between 0
#'   and 1. Defaults to `0.5`. A probability strictly greater than the
#'   threshold gives `"Positive"`; equality gives `"Negative"`.
#'
#' @return A named character vector of length two:
#'   * `predicted_m6A_prob`: the predicted probability of the `"Positive"`
#'     class, coerced to character by [c()].
#'   * `predicted_m6A_status`: either `"Positive"` or `"Negative"`.
#'
#'   Both elements are character strings, not numeric.
#'
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#'
#' prediction_single(ml_fit,
#'                   gc_content = 0.6368,
#'                   RNA_type = "mRNA",
#'                   RNA_region = "CDS",
#'                   exon_length = 10.7671,
#'                   distance_to_junction = 8.8826,
#'                   evolutionary_conservation = 0.2498,
#'                   DNA_5mer = "GGACA")
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}

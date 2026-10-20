#' Felinos Dataset
#'
#' A categorical dataset containing morphological and behavioral characteristics
#' of 30 feline species. This dataset is commonly used as a benchmark to test 
#' and compare categorical clustering algorithms (e.g., AAC, Hierarchical, DBSCAN).
#'
#' @format A data frame with 30 rows (species) and 14 categorical variables:
#' \describe{
#'   \item{tipopiel}{Skin type (categorical)}
#'   \item{longpoill}{Hair length (categorical)}
#'   \item{retract}{Retractable claws (categorical)}
#'   \item{comport}{Behavior (categorical)}
#'   \item{orielles}{Ears shape/type (categorical)}
#'   \item{larynx}{Larynx structure (categorical)}
#'   \item{tailler}{Size category (categorical)}
#'   \item{poids}{Weight category (categorical)}
#'   \item{longueurs}{Body length category (categorical)}
#'   \item{queue}{Tail type (categorical)}
#'   \item{dents}{Teeth characteristics (categorical)}
#'   \item{typproie}{Prey type preference (categorical)}
#'   \item{arbre}{Tree climbing ability (categorical)}
#'   \item{chasse}{Hunting style (categorical)}
#' }
#' @source Adapted from categorical clustering literature and benchmark datasets.
"Felinos"
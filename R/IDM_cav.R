


#' Cardiac allograft vasculopathy (CAV) data
#'
#' A series of approximately yearly angiographic examinations of heart transplant recipients. Onset of cardiac allograft vasculopathy, a deterioration of the arterial wall, and death are monitored. This is a subset of data from the \code{cav} dataset in \code{R} package \code{msm}.
#'
#' @format
#' A data frame with 2,803 observations of 614 patients and 5 variables. These are grouped by patient number and ordered by years after transplant.
#' \describe{
#'   \item{PTNUM}{Unique number identifying each patient.}
#'   \item{years}{Examination time (years after transplant).}
#'   \item{state}{State of the examination. State 1 represents no CAV, state 2 represents CAV, state 3 represents death.}
#'   \item{dage}{Age of the heart donor (years).}
#'   \item{pdiag}{Primary diagnosis (reason for transplant). IHD = ischaemic heart disease, IDC = idiopathic dilated cardiomyopathy.}
#' }
#' @source Papworth Hospital, U.K.
"IDM_cav"

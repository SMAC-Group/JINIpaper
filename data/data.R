#' Students' Alcohol Consumption Real Data
#'
#' Data collected during the 2005-2006 school year from
#' two public schools in Portugal. The data comes from two sources:
#' school recordings (e.g. grades, number of school absences) and
#' self-reporting questionnaires (e.g. workday and weekend alcohol consumptions,
#' parents' jobs, quality of family relationships, frequency to go out with friends).
#' It was originally collected and first studied in \insertCite{cortez2008using;textual}{JINIpaper}.
#'
#' @format A data frame with 395 observations and 45 variables:
#' \describe{
#' \item{\code{alc}}{student's alcohol consumption (binary: 0 = low, 1 = high)}
#' \item{\code{schoolMS}}{student's school (binary: 0 = Gabriel Pereira, 1 = Mousinho da Silveira)}
#' \item{\code{sexM}}{student's gender (binary: 0 = female, 1 = male)}
#' \item{\code{age}}{standardized student's age}
#' \item{\code{addressU}}{student's home address type (binary: 0 = rural, 1 = urban)}
#' \item{\code{famsizeLE3}}{family size (binary: 0 = greater than 3, 1 = lower or equal to 3)}
#' \item{\code{PstatusT}}{parent's cohabitation status (binary: 0 = living apart, 1 = living together)}
#' \item{\code{Mjobat_home}}{mother's areas of professional activity (binary: 0 = others, 1 = at home)}
#' \item{\code{Mjobhealth}}{mother's areas of professional activity (binary: 0 = others, 1 = health care related)}
#' \item{\code{Mjobservices}}{mother's areas of professional activity (binary: 0 = others, 1 = public services)}
#' \item{\code{Mjobteacher}}{mother's areas of professional activity (binary: 0 = others, 1 = teacher)}
#' \item{\code{Fjobat_home}}{father's areas of professional activity (binary: 0 = others, 1 = at home)}
#' \item{\code{Fjobhealth}}{father's areas of professional activity (binary: 0 = others, 1 = health care related)}
#' \item{\code{Fjobservices}}{father's areas of professional activity (binary: 0 = others, 1 = public services)}
#' \item{\code{Fjobteacher}}{father's areas of professional activity (binary: 0 = others, 1 = teacher)}
#' \item{\code{reasoncourse}}{reason for choosing the school (binary: 0 = others, 1 = course preference)}
#' \item{\code{reasonhome}}{reason for choosing the school (binary: 0 = others, 1 = close to home)}
#' \item{\code{reasonreputation}}{reason for choosing the school (binary: 0 = others, 1 = school reputation)}
#' \item{\code{guardianfather}}{student's guardian (binary: 0 = others, 1 = father)}
#' \item{\code{guardianmother}}{student's guardian (binary: 0 = others, 1 = mother)}
#' \item{\code{traveltime2}}{home to school travel time (binary: 0 = others, 1 = 15 to 30 min)}
#' \item{\code{traveltime3}}{home to school travel time (binary: 0 = others, 1 = more than 30 min)}
#' \item{\code{studytime2}}{weekly study time (binary: 0 = others, 1 = between 2 to 5 hours)}
#' \item{\code{studytime3}}{weekly study time (binary: 0 = others, 1 = between 5 to 10 hours)}
#' \item{\code{studytime4}}{weekly study time (binary: 0 = others, 1 = more than 10 hours)}
#' \item{\code{failures1}}{number of past class failures (binary: 0 = others, 1 = one)}
#' \item{\code{failures2}}{number of past class failures (binary: 0 = others, 1 = two or more)}
#' \item{\code{schoolsupyes}}{extra educational school support (binary: 0 = no, 1 = yes)}
#' \item{\code{famsupyes}}{family educational support (binary: 0 = no, 1 = yes)}
#' \item{\code{paidyes}}{extra paid classes (binary: 0 = no, 1 = yes)}
#' \item{\code{activitiesyes}}{extra-curricular activities (binary: 0 = no, 1 = yes)}
#' \item{\code{nurseryyes}}{attended nursery school (binary: 0 = no, 1 = yes)}
#' \item{\code{higheryes}}{willing to take higher education (binary: 0 = no, 1 = yes)}
#' \item{\code{internetyes}}{Internet access at home (binary: 0 = no, 1 = yes)}
#' \item{\code{romanticyes}}{involved in a romantic relationship (binary: 0 = no, 1 = yes)}
#' \item{\code{famrel}}{standardized quality of family relationship (original data from 1 to 5)}
#' \item{\code{freetime2}}{free time after school (binary: 0 = others, 1 = low)}
#' \item{\code{freetime3}}{free time after school (binary: 0 = others, 1 = medium)}
#' \item{\code{freetime4}}{free time after school (binary: 0 = others, 1 = high)}
#' \item{\code{freetime5}}{free time after school (binary: 0 = others, 1 = very high)}
#' \item{\code{goout}}{standardized going out with friends variable (original data from 1 to 5)}
#' \item{\code{absences}}{standardized number of school absences (original data from 0 to 75)}
#' \item{\code{G1}}{standardized first period grade (original data from 3 to 19)}
#' \item{\code{G2}}{standardized second period grade (original data from 0 to 19)}
#' \item{\code{G3}}{standardized final period grade (original data from 0 to 20)}
#' }
#' @references
#' \insertAllCited{}
#' @source \url{https://pcortez.dsi.uminho.pt/}
#' @importFrom Rdpack reprompt
"student"

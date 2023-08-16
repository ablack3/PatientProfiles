


#' Add Charlson Comorbidity Index to a cohort table
#'
#' @param x A Generated Cohort Set
#'
#' @return A GeneratedCohort with a charlson_score column
#' @export
add_charlson_score <- function(x) {
  require(CDMConnector)

  checkmate::assert_class(x, "GeneratedCohortSet")
  checkmate::assert_class(x, "tbl")
  cdm <- attr(x, "cdm_reference")
  checkmate::assert_class(cdm, "cdm_reference")
  con <- attr(cdm, "dbcon")
  checkmate::assert_true(DBI::dbIsValid(con))
  assert_write_schema(cdm)
  write_schema <- attr(cdm, "write_schema")

  assert_tables(cdm, c("concept", "concept_ancestor", "condition_era"))

  charlson_scoring <- tibble::tribble(
    ~diag_category_id, ~diag_category_name, ~weight, ~concept_id,
     1, 'Myocardial infarction',       1, 4329847,
     2, 'Congestive heart failure',    1, 316139,
     3, 'Peripheral vascular disease', 1, 321052,
     4, 'Cerebrovascular disease',     1, c(381591, 434056),
     5, 'Dementia',                    1, 4182210,
     6, 'Chronic pulmonary disease',   1, 4063381,
     7, 'Rheumatologic disease',       1, c(257628, 134442, 80800, 80809, 256197, 255348),
     8, 'Peptic ulcer disease',        1, 4247120,
     9, 'Mild liver disease',          1, c(4064161, 4212540),
    10, 'Diabetes (mild to moderate)', 1, 201820,
    11, 'Diabetes with chronic complications', 2, c(443767, 442793),
    12, 'Hemoplegia or paralegia', 2, c(192606, 374022),
    13, 'Renal disease', 2, 4030518,
    14, 'Any malignancy', 2, 443392,
    15, 'Moderate to severe liver disease', 3, c(4245975, 4029488, 192680, 24966),
    16, 'Metastatic solid tumor', 6, 432851,
    17, 'AIDS', 6, 439727) %>%
    tidyr::unnest(concept_id) %>%
    dplyr::mutate(concept_id = as.integer(concept_id))

  tempname <- paste0("temp", floor(as.numeric(Sys.time())*10) %% 1e6, "_charlson")


  DBI::dbWriteTable(con,
                    inSchema(write_schema, tempname, dbms = dbms(con)),
                    charlson_scoring)

  charlson <- dplyr::tbl(con, inSchema(write_schema, tempname, dbms = dbms(con))) %>%
    dplyr::inner_join(cdm$concept_ancestor, by = c("concept_id" = "ancestor_concept_id")) %>%
    dplyr::select(diag_category_name, concept_id = "descendant_concept_id", "weight")

  conditions <- cdm$condition_era %>%
    dplyr::select(subject_id = "person_id",
                  concept_id = "condition_concept_id",
                  start_date = "condition_era_start_date")

  x %>%
    dplyr::inner_join(conditions, by = "subject_id") %>%
    dplyr::inner_join(charlson, by = "concept_id") %>%
    dplyr::filter(.data$start_date <= .data$cohort_start_date) %>%
    dplyr::select("cohort_definition_id",
                    "cohort_start_date",
                    "subject_id",
                    "diag_category_name",
                    "weight") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$cohort_definition_id,
                    .data$cohort_start_date,
                    .data$subject_id) %>%
    dplyr::summarise(charlson_score = sum(weight, na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(x, ., by = c("cohort_definition_id",
                                  "cohort_start_date",
                                  "subject_id")) %>%
    dplyr::mutate(charlson_score = coalesce(charlson_score, 0L)) %>%
    computeQuery()
}


library(CDMConnector) # using v1.1 (dev version)
library(dplyr)

example_datasets()

con <- DBI::dbConnect(duckdb::duckdb(), eunomia_dir("synthea-lung_cancer-10k"))
cdm <- cdm_from_con(con, cdm_schema = "main", write_schema = "main")

# look at common conditions to find lung cancer codes
cdm %>%
  cdm_flatten(domain = "condition") %>%
  count(observation_concept_id, observation_concept_name, sort = TRUE)

cdm <- generate_concept_cohort_set(cdm,
                                   concept_set = list(lung_cancer = c(4310703,4115276,3739564,4110591)),
                                   name = "cohort",
                                   overwrite = TRUE)
cdm$cohort

cdm$cohort <- add_charlson_score(cdm$cohort)

cdm$cohort %>%
  select(cohort_definition_id, subject_id, charlson_score)

cdm$cohort %>%
  pull(charlson_score) %>%
  summary

DBI::dbDisconnect(con, shutdown = TRUE)

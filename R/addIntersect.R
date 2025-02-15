# Copyright 2023 DARWIN EU (C)
#
# This file is part of PatientProfiles
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' It creates columns to indicate overlap information between two tables
#'
#' @param x Table with individuals in the cdm
#' @param cdm Object that contains a cdm reference. Use CDMConnector to obtain a
#' cdm reference.
#' @param tableName name of the cohort that we want to check for overlap
#' @param filterVariable the variable that we are going to use to filter (e.g.
#' cohort_definition_id)
#' @param filterId the value of filterVariable that we are interested in, it can
#' be a vector
#' @param idName the name of each filterId, must have same length than
#' filterId
#' @param value value of interest to add: it can be count, flag, date or time
#' @param window window to consider events of
#' @param indexDate Variable in x that contains the date to compute the
#' intersection.
#' @param targetStartDate date of reference in cohort table, either for start
#' (in overlap) or on its own (for incidence)
#' @param targetEndDate date of reference in cohort table, either for end
#' (overlap) or NULL (if incidence)
#' @param order last or first date to use for date/time calculations
#' @param nameStyle naming of the added column or columns, should include
#' required parameters
#' @param tablePrefix The stem for the permanent tables that will
#' be created. If NULL, temporary tables will be used throughout.
#'
#' @return table with added columns with overlap information
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(PatientProfiles)
#'
#' cdm <- mockPatientProfiles()
#' result <- cdm$cohort1 %>%
#'   addIntersect(
#'     cdm = cdm, tableName = "cohort2", value = "date"
#'   ) %>%
#'   dplyr::collect()
#' }
#'
addIntersect <- function(x,
                         cdm,
                         tableName,
                         value,
                         filterVariable = NULL,
                         filterId = NULL,
                         idName = NULL,
                         window = list(c(0, Inf)), # list
                         indexDate = "cohort_start_date",
                         targetStartDate = getStartName(tableName),
                         targetEndDate = getEndName(tableName),
                         order = "first",
                         nameStyle = "{value}_{id_name}_{window_name}",
                         tablePrefix = NULL) {
  # initial checks
  personVariable <- checkX(x)
  checkmate::assertCharacter(tableName, len = 1, any.missing = FALSE)
  checkCdm(cdm, tableName)
  personVariableTable <- checkX(cdm[[tableName]])
  extraValue <- checkValue(value, cdm[[tableName]], tableName)
  filterTbl <- checkFilter(filterVariable, filterId, idName, cdm[[tableName]])
  windowTbl <- checkWindow(window)
  checkVariableInX(indexDate, x)
  checkVariableInX(targetStartDate, cdm[[tableName]], FALSE, "targetStartDate")
  checkVariableInX(targetEndDate, cdm[[tableName]], TRUE, "targetEndDate")
  checkmate::assertChoice(order, c("first", "last"))
  checkNameStyle(nameStyle, filterTbl, windowTbl, value)
  checkmate::assertCharacter(tablePrefix, len = 1, null.ok = TRUE)
  if (!is.null(idName)) {
    idName <- checkSnakeCase(idName)
  }

  startTibble <- x
  originalColnames <- colnames(x)

  # define overlapTable that contains the events of interest
  overlapTable <- cdm[[tableName]]
  if (!is.null(filterTbl)) {
    overlapTable <- overlapTable %>%
      dplyr::filter(.data[[filterVariable]] %in% .env$filterId)
  } else {
    filterVariable <- "id"
    filterTbl <- dplyr::tibble("id" = 1, "id_name" = "all")
    overlapTable <- dplyr::mutate(overlapTable, "id" = 1)
  }
  if (is.null(targetEndDate)) {
    overlapTable <- overlapTable %>%
      dplyr::select(
        !!personVariable := dplyr::all_of(personVariableTable),
        "id" = dplyr::all_of(filterVariable),
        "overlap_start_date" = dplyr::all_of(targetStartDate),
        "overlap_end_date" = dplyr::all_of(targetStartDate),
        dplyr::all_of(extraValue)
      )
  } else {
    overlapTable <- overlapTable %>%
      dplyr::select(
        !!personVariable := dplyr::all_of(personVariableTable),
        "id" = dplyr::all_of(filterVariable),
        "overlap_start_date" = dplyr::all_of(targetStartDate),
        "overlap_end_date" = dplyr::all_of(targetEndDate),
        dplyr::all_of(extraValue)
      )
  }

  result <- x %>%
    dplyr::select(
      dplyr::all_of(personVariable),
      "index_date" = dplyr::all_of(indexDate)
    ) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(overlapTable, by = personVariable)

  if (is.null(tablePrefix)) {
    result <- CDMConnector::computeQuery(result)
  } else {
    result <- CDMConnector::computeQuery(
      result, paste0(tablePrefix, "_join"), FALSE, attr(cdm, "write_schema"), TRUE
    )
  }

  resultCountFlag <- NULL
  resultDateTimeOther <- NULL
  # Start loop for different windows

  for (i in c(1:nrow(windowTbl))) {
    resultW <- result
    if (!is.infinite(windowTbl$upper[i])) {
      resultW <- resultW %>%
        dplyr::mutate(indicator = dplyr::if_else(.data$index_date >= as.Date(!!CDMConnector::dateadd(
          date = "overlap_start_date", number = -windowTbl$upper[i]
        )), 1, 0))
    } else {
      resultW <- resultW %>% dplyr::mutate(indicator = 1)
    }

    if (!is.infinite(windowTbl$lower[i])) {
      resultW <- resultW %>%
        dplyr::mutate(indicator = dplyr::if_else(.data$index_date > as.Date(!!CDMConnector::dateadd(
          date = "overlap_end_date", number = -windowTbl$lower[i]
        )), 0, .data$indicator))
    }
    if (is.null(tablePrefix)) {
      resultW <- CDMConnector::computeQuery(resultW)
    } else {
      resultW <- CDMConnector::computeQuery(
        resultW, paste0(tablePrefix, "_window"), FALSE, attr(cdm, "write_schema"), TRUE
      )
    }
    # add count or flag
    if ("count" %in% value | "flag" %in% value) {
      resultCF <- resultW %>%
        dplyr::group_by(.data[[personVariable]], .data$index_date, .data$id) %>%
        dplyr::summarise(count = sum(.data$indicator, na.rm = TRUE), .groups = "drop") %>%
        dplyr::left_join(filterTbl, by = "id", copy = TRUE) %>%
        dplyr::select(-"id") %>%
        dplyr::mutate("window_name" = !!tolower(windowTbl$window_name[i]))
      if ("flag" %in% value) {
        resultCF <- resultCF %>% dplyr::mutate(flag = dplyr::if_else(.data$count > 0, 1, 0))
      }
      if (!("count" %in% value)) {
        resultCF <- dplyr::select(resultCF, -"count")
      }
      if (is.null(tablePrefix)) {
        resultCF <- CDMConnector::computeQuery(resultCF)
      } else {
        resultCF <- CDMConnector::computeQuery(
          resultCF, paste0(tablePrefix, "_count_flag_", i), FALSE,
          attr(cdm, "write_schema"), TRUE
        )
      }
      if (i == 1) {
        resultCountFlag <- resultCF
      } else {
        resultCountFlag <- dplyr::union_all(resultCountFlag, resultCF)
      }
    }
    # add date, time or other
    if (length(value[!(value %in% c("count", "flag"))]) > 0) {
      resultDTO <- resultW %>%
        dplyr::filter(.data$indicator == 1) %>%
        dplyr::group_by(.data[[personVariable]], .data$index_date, .data$id)
      if (order == "first") {
        resultDTO <- resultDTO %>%
          dplyr::summarise(
            date = min(.data$overlap_start_date, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        resultDTO <- resultDTO %>%
          dplyr::summarise(
            date = max(.data$overlap_start_date, na.rm = TRUE),
            .groups = "drop"
          )
      }
      resultDTO <- resultDTO %>%
        dplyr::right_join(
          resultW %>%
            dplyr::select(dplyr::all_of(c(personVariable, "index_date", "id"))) %>%
            dplyr::distinct(),
          by = c(personVariable, "index_date", "id")
        )
      if ("days" %in% value) {
        resultDTO <- resultDTO %>%
          dplyr::mutate(
            days = !!CDMConnector::datediff("index_date", "date", interval = "day")
          )
      }
      if (length(extraValue) > 0) {
        resultDTO <- resultDTO %>%
          dplyr::left_join(
            resultW %>%
              dplyr::select(
                dplyr::all_of(personVariable), "index_date", "id",
                "date" = "overlap_start_date", dplyr::all_of(extraValue)
              ) %>%
              dplyr::inner_join(
                resultDTO %>%
                  dplyr::select(dplyr::all_of(
                    c(personVariable, "index_date", "id", "date")
                  )),
                by = c(personVariable, "index_date", "id", "date")
              ) %>%
              dplyr::group_by(.data[[personVariable]], .data$index_date, .data$id) %>%
              dplyr::summarise(
                dplyr::across(
                  dplyr::all_of(extraValue), ~ str_flatten(.x, collapse = "; ")
                ),
                .groups = "drop"
              ),
            by = c(personVariable, "index_date", "id")
          )
      }
      resultDTO <- resultDTO %>%
        dplyr::left_join(filterTbl, by = "id", copy = TRUE) %>%
        dplyr::select(-"id") %>%
        dplyr::mutate("window_name" = !!tolower(windowTbl$window_name[i]))
      if (!("date" %in% value)) {
        resultDTO <- dplyr::select(resultDTO, -"date")
      }
      if (is.null(tablePrefix)) {
        resultDTO <- CDMConnector::computeQuery(resultDTO)
      } else {
        resultDTO <- CDMConnector::computeQuery(
          resultDTO, paste0(tablePrefix, "_date_time_", i), FALSE,
          attr(cdm, "write_schema"), TRUE
        )
      }
      if (i == 1) {
        resultDateTimeOther <- resultDTO
      } else {
        resultDateTimeOther <- dplyr::union_all(resultDateTimeOther, resultDTO)
      }
    }
  }

  if (any(c("flag", "count") %in% value)) {
    resultCountFlag <- resultCountFlag %>%
      tidyr::pivot_longer(
        dplyr::any_of(c("count", "flag")),
        names_to = "value",
        values_to = "values"
      ) %>%
      tidyr::pivot_wider(
        names_from = c("value", "id_name", "window_name"),
        values_from = "values",
        names_glue = nameStyle,
        values_fill = 0
      ) %>%
      dplyr::rename(!!indexDate := "index_date") %>%
      dplyr::rename_all(tolower)

    namesToEliminate <- intersect(names(x), names(resultCountFlag))
    namesToEliminate <- namesToEliminate[
      !(namesToEliminate %in% c(personVariable, indexDate))
    ]
    x <- x %>%
      dplyr::select(-dplyr::all_of(namesToEliminate)) %>%
      dplyr::left_join(
        resultCountFlag,
        by = c(personVariable, indexDate)
      )
    currentColnames <- colnames(x)
    x <- x %>%
      dplyr::mutate(dplyr::across(
        dplyr::all_of(currentColnames[!(currentColnames %in% originalColnames)]),
        ~ dplyr::if_else(is.na(.x), 0, .x)
      ))
    if (is.null(tablePrefix)) {
      x <- CDMConnector::computeQuery(x)
    } else {
      x <- CDMConnector::computeQuery(
        x, paste0(tablePrefix, "_intersect"), FALSE, attr(cdm, "write_schema"),
        TRUE
      )
    }
  }

  if (length(value[!(value %in% c("count", "flag"))]) > 0) {
    values <- value[!(value %in% c("count", "flag"))]
    for (val in values) {
      resultDateTimeOtherX <- resultDateTimeOther %>%
        dplyr::select(
          "subject_id", "index_date", dplyr::all_of(val), "id_name",
          "window_name"
        ) %>%
        tidyr::pivot_longer(
          dplyr::all_of(val),
          names_to = "value",
          values_to = "values"
        ) %>%
        tidyr::pivot_wider(
          names_from = c("value", "id_name", "window_name"),
          values_from = "values",
          names_glue = nameStyle
        ) %>%
        dplyr::rename(!!indexDate := "index_date") %>%
        dplyr::rename_all(tolower)

      namesToEliminate <- intersect(names(x), names(resultDateTimeOtherX))
      namesToEliminate <- namesToEliminate[
        !(namesToEliminate %in% c(personVariable, indexDate))
      ]

      x <- x %>%
        dplyr::select(-dplyr::all_of(namesToEliminate)) %>%
        dplyr::left_join(resultDateTimeOtherX,
          by = c(personVariable, indexDate)
        )
    }

    if (is.null(tablePrefix)) {
      x <- CDMConnector::computeQuery(x)
    } else {
      x <- CDMConnector::computeQuery(
        x, paste0(tablePrefix, "_intersect"), FALSE, attr(cdm, "write_schema"),
        TRUE
      )
    }
  }

  colnames <- expand.grid(value = value, id_name = filterTbl$id_name, window_name = windowTbl$window_name) %>%
    dplyr::mutate(column = glue::glue(nameStyle, value = .data$value, id_name = .data$id_name, window_name = .data$window_name)) %>%
    dplyr::mutate(val = ifelse(value %in% c("flag", "count"), 0,
      ifelse(value %in% "date", as.Date(NA),
        ifelse(value %in% "days", as.numeric(NA), as.character(NA))
      )
    )) %>%
    dplyr::select(.data$column, .data$val) %>%
    dplyr::anti_join(dplyr::tibble(column = colnames(x)), by = "column")

  if (colnames %>% dplyr::tally() %>% dplyr::pull() != 0) {
    x <- x %>%
      dplyr::cross_join(
        colnames %>%
          tidyr::pivot_wider(
            names_from = .data$column,
            values_from = .data$val
          ),
        copy = TRUE
      )
  }

  if (is.null(tablePrefix)) {
    x <- CDMConnector::computeQuery(x)
  } else {
    x <- CDMConnector::computeQuery(
      x, paste0(tablePrefix, "_intersect"), FALSE, attr(cdm, "write_schema"),
      TRUE
    )
  }

  # put back the initial attributes to the output tibble
  x <- x %>% addAttributes(startTibble)

  return(x)
}

#' Get the name of the start date column for a certain table in the cdm
#'
#' @param tableName Name of the table
#'
#' @return Name of the start date column in that table
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(PatientProfiles)
#' getStartName("condition_occurrence")
#' }
#'
getStartName <- function(tableName) {
  if (tableName %in% namesTable$table_name) {
    return(namesTable$start_date_name[namesTable$table_name == tableName])
  } else {
    return("cohort_start_date")
  }
}

#' Get the name of the end date column for a certain table in the cdm
#'
#' @param tableName Name of the table
#'
#' @return Name of the end date column in that table
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(PatientProfiles)
#' getEndName("condition_occurrence")
#' }
#'
getEndName <- function(tableName) {
  if (tableName %in% namesTable$table_name) {
    return(namesTable$end_date_name[namesTable$table_name == tableName])
  } else {
    return("cohort_end_date")
  }
}

#' Get the name of the concept_id column for a certain table in the cdm
#'
#' @param tableName Name of the table
#'
#' @return Name of the concept_id column in that table
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(PatientProfiles)
#' getConceptName("condition_occurrence")
#' }
#'
getConceptName <- function(tableName) {
  if (tableName %in% namesTable$table_name) {
    return(namesTable$concept_id_name[namesTable$table_name == tableName])
  } else {
    return("cohort_definition_id")
  }
}

library(readxl)
library(dplyr)
library(doParallel)
library(tibble)
library(tidyr)
library(stringr)

## Parallel computing
registerDoParallel(cores = detectCores())

## read sub-directories for class names
directory <- "Data"
classes <- gsub(file.path(directory, "/"), "", list.dirs(directory)[c(-1)], fixed=TRUE)

## check file extension
sample_file_ext <- unique(unlist(
  lapply(classes, function(x) 
    tools::file_ext(list.files(file.path(directory, "/",x,"/"))))))

## read names of files -- xlsx
file_name_list <- lapply(classes, function(x) 
  list.files(path = file.path(directory,x), pattern = "*.xlsx", full.names = TRUE))

sample_names <- as.data.frame(unlist(lapply(classes, function(x) 
  tools::file_path_sans_ext(list.files(path = file.path(directory,x), 
                                       pattern = "*.xlsx", full.names = FALSE)))))

## Extract m/z and intensity into list of lists
peak_list <- lapply(file_name_list, function(x) foreach(i = 1:length(x)) %do% 
                      {read_excel(x[[i]],col_names = TRUE, skip = 7)}) 
## Target m/z values
pos_mz_values <- c(194.079, 197.0709, 211.0866, 216.1383, 232.1332, 242.1539, 244.1696, 
                   258.1852, 260.1645, 270.1852, 272.2009, 274.1802, 276.1594, 286.1802, 
                   288.1958, 298.2165, 300.2322, 304.1907, 312.2322, 314.2115, 320.183, 
                   326.2478)

neg_mz_values <- c(475.291, 503.322, 529.33845, 531.354, 649.3801, 675.3961, 677.412, 703.428, 
                   705.4437, 714.5079, 716.5236, 728.5236, 745.5025, 747.5182, 759.5182,
                   761.5338, 773.5338, 774.4846)

tolerance <- 0.005

## Filter peaks based on sample name aka Pos or Neg
filter_peaks_mode <- function(spectrum, file_name) {
  ## Determine mode from file name
  mode_indicator <- if(grepl("Pos", file_name, ignore.case = TRUE)) "Pos" else "Neg"
  ## Select m/z values based on mode
  mz_values <- if(mode_indicator == "Pos") pos_mz_values else neg_mz_values
  ## Filter
  filtered_peaks <- filter_mz_values(spectrum, mz_values, tolerance)
  return(filtered_peaks)
}

## Filter m/z values within tolerance
filter_mz_values <- function(spectrum, mz_peaks, tolerance) {
  ## Identify all "Mass" columns
  mass_columns <- grep("Mass", names(spectrum), value = TRUE)
  if(length(mass_columns) == 0 || is.null(mz_peaks) || length(mz_peaks) == 0) {
    return(data.frame())
  }
  ## Initialize an empty data frame to accumulate matches
  all_matched_peaks <- data.frame()
  
  for(mass_col in mass_columns) {
    ## Apply filtering for each target mz value
    for(target_mz in mz_peaks) {
      valid_mass_values <- !is.na(spectrum[[mass_col]])
      matched <- valid_mass_values & abs(spectrum[[mass_col]] - target_mz) <= tolerance
      
      if(any(matched, na.rm = TRUE)) {
        ## Extract matched peaks
        matched_peaks <- spectrum[matched, c(mass_col, grep("Intensity", names(spectrum), value = TRUE))]
        ## Add a column for the target m/z value for grouping
        matched_peaks$TargetMZ <- target_mz
        
        all_matched_peaks <- dplyr::bind_rows(all_matched_peaks, matched_peaks)
      }
    }
  }
  
  if(nrow(all_matched_peaks) == 0) {
    return(data.frame())
  }
  
  ## Identify all intensity columns
  intensity_cols <- grep("Intensity", names(all_matched_peaks), value = TRUE)
  
  ## Group by TargetMZ and conditionally sum intensities for each intensity column
  aggregated_peaks <- all_matched_peaks %>%
    dplyr::group_by(TargetMZ) %>%
    dplyr::summarise(across(all_of(intensity_cols), ~if(all(is.na(.x))) {NA_real_} else {sum(.x, na.rm = TRUE)}), .groups = 'drop')
  
  return(aggregated_peaks)
}

## Flatten sample_names to match peak_list structure
flat_sample_names <- unlist(sample_names)

## Flatten peak_list to match flat_sample_names structure
flat_peak_list <- unlist(peak_list, recursive = FALSE)

## Match file names to corresponding tibbles
names(flat_peak_list) <- flat_sample_names

## Filter each file to peak list using sample_names
matched_mz <- mapply(function(spectrum, file_name) {
  filter_peaks_mode(spectrum, file_name)
}, flat_peak_list, flat_sample_names, SIMPLIFY = FALSE)

## Save/load matched_mz
#save(matched_mz, file="matched_mz.RData")
#load("matched_mz.RData")

## Function to add missing m/z rows to a single tibble
add_missing_mz_rows_to_tibble <- function(tibble, target_mz_values) {
  ## Create a data frame of all possible target m/z values
  all_mz <- tibble(TargetMZ = target_mz_values)
  
  ## Join the all_mz frame with the existing tibble to find missing values
  ## This uses a left join to ensure all target m/z values are represented
  combined <- all_mz %>%
    left_join(tibble, by = "TargetMZ")
  
  ## Replace NA values in intensity columns with 0
  combined[is.na(combined)] <- 0
  
  return(combined)
}

## Infer mode by sample name
infer_mode_from_sample_name <- function(sample_name) {
  if(grepl("Pos", sample_name, ignore.case = TRUE)) {
    return("Pos")
  } else if(grepl("Neg", sample_name, ignore.case = TRUE)) {
    return("Neg")
  } else {
    return(NA) # Undefined mode
  }
}

## Fill in missing m/z rows based on inferred mode:
matched_mz <- lapply(names(matched_mz), function(name) {
  tibble <- matched_mz[[name]]
  mode <- infer_mode_from_sample_name(name) # Adjust this line to use your mode detection logic
  
  target_mz_values <- if (mode == "Pos") pos_mz_values else neg_mz_values
  
  add_missing_mz_rows_to_tibble(tibble, target_mz_values)
})

## Ensure matched_mz list is named by sample names
names(matched_mz) <- flat_sample_names

## Normalize
medianNonZero <- function(x) {
  median(x[x > 0], na.rm = TRUE)
}
normalize_intensities <- function(tibble) {
  intensity_cols <- grep("Intensity", names(tibble), value = TRUE)
  
  for(col_name in intensity_cols) {
    median_intensity <- median(tibble[[col_name]][tibble[[col_name]] > 0], na.rm = TRUE)
    
    # Avoid division by zero if median_intensity is zero
    if(median_intensity > 0) {
      scaling_factor <- 10000 / median_intensity
      tibble[[col_name]] <- tibble[[col_name]] * scaling_factor
    }
  }
  
  return(tibble)
}
matched_mz_norm <- lapply(matched_mz, normalize_intensities)

## Group tibbles by sample name to combine
base_names <- gsub("_(Neg|Pos)$", "", names(matched_mz_norm))
grouped_names <- split(names(matched_mz_norm), base_names)

## Function to combine Pos and Neg tibbles
combine_pos_neg_tibbles <- function(tibble_names, matched_mz_norm) {
  pos_tibble <- matched_mz_norm[[tibble_names[grepl("_Pos$", tibble_names)]]]
  neg_tibble <- matched_mz_norm[[tibble_names[grepl("_Neg$", tibble_names)]]]
  ## Remove Pos/Neg from column names
  names(pos_tibble) <- gsub("_Pos$", "", names(pos_tibble))
  names(neg_tibble) <- gsub("_Neg$", "", names(neg_tibble))
  ## Merge columns
  combined_tibble <- bind_rows(pos_tibble, neg_tibble)
  return(combined_tibble)
}

## Combine tibbles for each base name
combined_matched_mz <- lapply(grouped_names, function(tibble_names) combine_pos_neg_tibbles(tibble_names, matched_mz_norm))

## Save/load matched_mz
#save(matched_mz, matched_mz_norm, combined_matched_mz, file="matched_mz.RData")
#load("matched_mz.RData")

## Statistical analysis
## Add group (PA vs RPA) and condition (Control, 25, 50, 100) columns to tibbles
preprocess_tibble <- function(tibble, name) {
  group <- ifelse(grepl("RPA_", name), "RPA", "PA")
  condition <- str_extract(name, "\\d+|Control")
  tibble <- tibble %>%
    mutate(Group = group, Condition = condition)
  return(tibble)
}

combined_matched_mz <- lapply(names(combined_matched_mz), function(name) {
  tibble <- combined_matched_mz[[name]]
  preprocess_tibble(tibble, name)
})

names(combined_matched_mz) <- unique(base_names)

## Combine into a single matrix
combined_matrix <- bind_rows(combined_matched_mz)

## Reshape matrix to long format
long_combined_matrix <- combined_matrix %>%
  ## Exclude unnamed column number 1
  pivot_longer(
    cols = starts_with("Intensity"),
    names_to = "Sample",
    names_prefix = "Intensity...",
    values_to = "Intensity"
  )

## Prep for Kruskal-Wallis Test (these data are not normally distributed and tried log and Box-Cox transforms)
long_combined_matrix$Group <- as.factor(long_combined_matrix$Group)
long_combined_matrix$Condition <- as.factor(long_combined_matrix$Condition)
unique_mz <- unique(long_combined_matrix$TargetMZ)

## Perform Kruskal-Wallis Test
kruskalwallis_results <- list()
for (mz in unique_mz) {
  data_subset <- filter(long_combined_matrix, TargetMZ == mz)
  ## Combine Group and Condition into a single factor
  interaction_groupandcondition <- interaction(data_subset$Group, data_subset$Condition)
  kruskalwallis_result <- kruskal.test(data_subset$Intensity ~ interaction_groupandcondition)
  kruskalwallis_results[[as.character(mz)]] <- kruskalwallis_result
}
kruskalwallis_results
## Significant peaks: 216.1383, 304.1907, 312.2322, 475.291, 529.33845, 531.354,
##                    649.3801, 675.3961, 677.412, 703.428, 705.4437, 745.5025, 773.5338

## Dunn's test to determine which groups are significant for each m/z peak
library(FSA)
kruskal_dunn_results <- list()
for (mz in unique_mz) {
  ## Extract data for each m/z peak
  data_subset <- filter(long_combined_matrix, TargetMZ == mz)
  ## Check if the Kruskal-Wallis test was significant (0.05 threshold)
  if (kruskalwallis_results[[as.character(mz)]]$p.value < 0.05) {
    ## Create interaction term for Group and Condition
    data_subset$GroupCondition <- interaction(data_subset$Group, data_subset$Condition)
    ## Perform Dunn's test
    dunn_test_result <- dunnTest(x = data_subset$Intensity,
                                 g = data_subset$GroupCondition,
                                 method = "holm")
    ## Store result
    dunn_results[[as.character(mz)]] <- dunn_test_result
  }
}
dunn_results

## Export Kruskal-Wallis and Dunn Test Results
library(openxlsx)
results_excel <- createWorkbook("Dunn Test Results")
for (mz in names(dunn_results)){
  sheetName <- addWorksheet(results_excel, mz)
  sheetName
  dunn_test_result <- dunn_results[[mz]]
  dunn_df <- data.frame(
    Comparison = dunn_test_result$res$Comparison,
    Z = dunn_test_result$res$Z,
    P.adj = dunn_test_result$res$P.adj
  )
  writeData(results_excel, sheetName, dunn_df)
}
saveWorkbook(results_excel, file = "pseudomonas_dunn_results.xlsx", overwrite = TRUE)


#' Imports a collated dataset into the R session as an NxtSE object
#'
#' Creates a \linkS4class{NxtSE} object from the data (that was collated using 
#' [collateData]). This object is used
#' for downstream differential analysis of IR and alternative splicing events
#' using [ASE-methods], data generation for visualization of scatter plots and
#' heatmaps via [make_plot_data] methods, and coverage visualisation using 
#' [plotCoverage]
#'
#' @details
#' `makeSE` retrieves the data collated by
#' [collateData], and initialises a \linkS4class{NxtSE} object. It references
#' the required on-disk assay data using DelayedArrays, thereby utilising
#' 'on-disk' memory to conserve memory usage.
#'
#' For extremely large datasets, loading the entire data into memory may consume
#' too much memory. In such cases, make a subset of the \linkS4class{NxtSE}
#' object (e.g. subset by samples) before loading the data into memory (RAM) 
#' using [realize_NxtSE]. Alternatively supply a data frame to the `colData`
#' parameter of the `makeSE()` function. Only samples listed in the first column
#' of the `colData` data frame will be imported into the \linkS4class{NxtSE} 
#' object.
#' 
#' It should be noted that downstream applications of SpliceWiz, including
#' [ASE-methods], [plotCoverage], are much faster if the \linkS4class{NxtSE}
#' is realized. It is recommended to realize the \linkS4class{NxtSE} object
#' before extensive usage.
#'
#' If COV files assigned via [collateData] have been moved relative to the
#' `collate_path`, the created \linkS4class{NxtSE} object will not be linked to
#' any COV files and [plotCoverage] cannot be used. To reassign these
#' files, a vector of file paths corresponding to all the COV files of the data
#' set can be assigned using `covfile(se) <- vector_of_cov_files`. See
#' the example below for details.
#'
#' If `RemoveOverlapping = TRUE`, `makeSE` will remove introns that overlap
#' other introns with higher junction read counts in the dataset. This means
#' that SpliceWiz will assess a set of non-overlapping introns which belong
#' to likely major isoforms, ensuring that overlapping
#' IR events are not 'double-counted'.
#'
#' NB: Since version 1.3.4, SpliceWiz has improved the algorithm of generating
#' the set of non-overlapping introns (prior versions appear to generate
#' sets of introns that still overlap). To use the prior algorithm for
#' compatibility with prior analysis, set `RemoveOverlapping = FALSE`.
#'
#' @param collate_path (Required) The output path of [collateData] pointing
#'   to the collated data
#' @param colData (Optional) A data frame containing the sample annotation
#'   information. The first column must contain the sample names.
#'   Omit `colData` to generate a NxtSE object of the whole dataset without
#'   any assigned annotations.
#'   Alternatively, if the names of only a subset of samples are given, then
#'   `makeSE()` will construct the NxtSE object based only on the samples given.
#'   The colData can be set later using [colData]
#' @param RemoveOverlapping (default = `TRUE`) Whether to filter out overlapping
#'   IR events belonging to minor isoforms. See details.
#' @param realize (default = `FALSE`) Whether to load all assay data into
#'   memory. See details
#' @param verbose (default = `TRUE`) Whether loading messages are displayed
#'
#' @return A \linkS4class{NxtSE} object containing the compiled data in
#' DelayedArrays (or as matrices if `realize = TRUE`), pointing to the assay 
#' data contained in the given `collate_path`
#'
#' @examples
#'
#' # The following code can be used to reproduce the NxtSE object
#' # that can be fetched with SpliceWiz_example_NxtSE()
#'
#' buildRef(
#'     reference_path = file.path(tempdir(), "Reference"),
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#'
#' bams <- SpliceWiz_example_bams()
#' processBAM(bams$path, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "SpliceWiz_Output")
#' )
#'
#' expr <- findSpliceWizOutput(file.path(tempdir(), "SpliceWiz_Output"))
#' collateData(expr,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "Collated_output")
#' )
#'
#' se <- makeSE(collate_path = file.path(tempdir(), "Collated_output"))
#'
#' # "Realize" NxtSE object to load all H5 assays into memory:
#'
#' se <- realize_NxtSE(se)
#'
#' # If COV files have been removed since the last call to collateData()
#' # reassign them to the NxtSE object, for example:
#'
#' covfile_path <- system.file("extdata", package = "swCore")
#' covfile_df <- findSamples(covfile_path, ".cov")
#'
#' covfile(se) <- covfile_df$path
#'
#' @md
#' @export
makeSE <- function(
        collate_path, colData, RemoveOverlapping = TRUE,
        realize = FALSE, verbose = TRUE
) {
    # How many steps
    N <- 7

    # Includes iterative filtering for IR events with highest mean PSI
        # To annotate IR events of major isoforms

    fullExperiment <- FALSE
    if(missing(colData)) fullExperiment <- TRUE

    colData <- .makeSE_validate_args(collate_path, colData)
    colData <- .makeSE_colData_clean(colData)

    collate_path <- normalizePath(collate_path)

    se <- .makeSE_load_NxtSE(collate_path, N, verbose)
    # Subset, if required
    if(all(colnames(se) %in% colData$sample)) {
        fullExperiment <- TRUE
    } else {
        se <- se[, colData$sample]
    }
    if (ncol(colData) > 1) {
        colData_use <- colData[, -1, drop = FALSE]
        rownames(colData_use) <- colData$sample
        colData(se) <- as(colData_use, "DataFrame")
    }

    if (realize == TRUE) {
        dash_progress("Realizing NxtSE object...", N)
        if(verbose) .log("...realizing NxtSE object", "message")
        se <- realize_NxtSE(se)
    }
    
    filtered_rowData_file <- file.path(collate_path, "filteredIntrons.Rds")
    if (RemoveOverlapping == TRUE) {
        dash_progress("Removing overlapping introns...", N)
        if(verbose) .log("...removing overlapping introns...", "message")

        if(
            metadata(se)$BuildVersion >= "1.3.4" && 
            fullExperiment & file.exists(filtered_rowData_file)
        ) {
            # Take a shortcut
            tmpFiltered <- readRDS(filtered_rowData_file)
            se <- se[tmpFiltered,]
        } else {
            se <- .makeSE_iterate_IR_new(se, verbose)
        } 
    } else if(RemoveOverlapping == "legacy") {
        dash_progress("Removing overlapping introns...", N)
        if(verbose) .log("...removing overlapping introns...", "message")

        if(
            metadata(se)$BuildVersion < "1.3.4" &
            metadata(se)$BuildVersion >= "1.1.6" &            
            fullExperiment & file.exists(filtered_rowData_file)
        ) {
            # Take a shortcut
            tmpFiltered <- readRDS(filtered_rowData_file)
            se <- se[tmpFiltered,]
        } else {
            se <- .makeSE_iterate_IR(se, verbose)
        }
    }

    if(verbose) .log("NxtSE loaded", "message")
    return(se)
}

################################################################################

# Checks:
# - whether the given path contains a valid collateData() output
# - whether
.makeSE_validate_args <- function(collate_path, colData) {
    item.todo <- c("rowEvent", "Included", "Excluded", "Depth", "Coverage",
        "minDepth", "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")

    if (!file.exists(file.path(collate_path, "colData.Rds"))) {
        .log(paste("In makeSE():",
            file.path(collate_path, "colData.Rds"),
            "was not found"))
    }
    colData.Rds <- readRDS(file.path(collate_path, "colData.Rds"))
    if (!("df.anno" %in% names(colData.Rds))) {
        .log(paste("In makeSE():",
            file.path(collate_path, "colData.Rds"),
            "must contain df.anno containing annotations"))
    }
    if (missing(colData)) {
        colData <- colData.Rds$df.anno
    } else {
        colData <- as.data.frame(colData)
        if (!("sample" %in% colnames(colData)))
            colnames(colData)[1] <- "sample"

        # Check all names unique
        if(!identical(colData$sample, unique(colData$sample)))
            .log("All sample names in colData must be unique.")

        # Make sample the 1st column if it is not
        whichColIsSample <- which(colnames(colData) == "sample")
        colData <- cbind(
            colData[, whichColIsSample, drop = FALSE],
            colData[, -whichColIsSample, drop = FALSE]
        )

        if (!all(colData$sample %in% colData.Rds$df.anno$sample)) {
            .log(paste("In makeSE():",
                "some samples in colData were not found in given path"),
                "message")
            colData <- colData[colData$sample %in% colData.Rds$df.anno$sample, ]
        }
    }
    return(colData)
}

# Converts charactor vectors to factors, removes columns with all NA's
.makeSE_colData_clean <- function(colData) {
    remove_na <- NULL
    if (ncol(colData) > 1) {
        for (i in seq(2, ncol(colData))) {
            if (is(colData[, i], "character")) {
                colData[, i] <- factor(unlist(colData[, i]))
            } else if (is(colData[, i], "logical")) {
                colData[, i] <- factor(unlist(
                    ifelse(colData[, i], "TRUE", "FALSE")))
            } else if (all(is.na(unlist(colData[, i])))) {
                remove_na <- append(remove_na, i)
            }
        }
    }
    if (!is.null(remove_na)) {
        colData <- colData[, -remove_na]
    }
    return(colData)
}

# Loads a NxtSE RDS
.makeSE_load_NxtSE <- function(collate_path, N, verbose = TRUE) {
    dash_progress("Loading NxtSE object from file...", N)
    if(verbose) .log("Loading NxtSE object from file...", "message")

    filepath <- file.path(collate_path, "seed.Rds")
    se <- readRDS(filepath)
    se@metadata[["sourcePath"]] <- collate_path

    if(
            !("BuildVersion" %in% names(se@metadata)) ||
            metadata(se)$BuildVersion < collateData_version
    ) {
        .log(paste(
            "NB: This NxtSE was generated with an older version of SpliceWiz",
            "(version",
            ifelse(
                "BuildVersion" %in% names(se@metadata),
                metadata(se)$BuildVersion,
                "< 0.99.4"
            ), ")"
        ), "warning")
    }

    # Add reference
    dash_progress("Loading reference (NxtSE)...", N)
    if(verbose) .log("...loading reference (NxtSE)", "message")
    metadata(se)$ref <- readRDS(file.path(collate_path, "cov_data.Rds"))

    dash_progress("Linking assays...", N)
    if(verbose) .log("...linking assays", "message")
    se@assays <- .makeSE_adjust_paths(se)
    
    metadata_opts <- c("Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")
    junc_opts <- c("junc_PSI", "junc_counts", "junc_counts_uns")

    for(mopt in metadata_opts) {
        if(mopt %in% names(se@metadata)) {
            se@metadata[[mopt]] <- .makeSE_expand_assay_path(
                se@metadata[[mopt]], collate_path)
        }
    }
    for(jopt in junc_opts) {
        if(jopt %in% names(se@metadata)) {
            se@metadata[[jopt]] <- .makeSE_expand_assay_path(
                se@metadata[[jopt]], collate_path)
        }
    }

    # Locate relative paths of COV files, or have all-empty if not all are found
    dash_progress("Linking COV files...", N)
    if(verbose) .log("...linking COV files", "message")
    suppressWarnings({
        covfiles <- normalizePath(
            file.path(collate_path, se@metadata[["cov_file"]])
        )    
    })
    
    display_cov_missing_message <- FALSE
    if (!all(file.exists(covfiles))) {
        se@metadata[["cov_file"]] <- rep("", ncol(se))
        display_cov_missing_message <- TRUE
    }
    names(se@metadata[["cov_file"]]) <- colnames(se)

    # Import rowData
    dash_progress("Loading rowData...", N)
    if(verbose) .log("...loading rowData", "message")
    rowData <- readRDS(file.path(collate_path, "rowEvent.Rds"))
    rowData(se) <- rowData
    metadata(se)[["row_gr"]] <- coord2GR(rowData$EventRegion)

    # Encapsulate as NxtSE object
    se <- as(se, "NxtSE")

    if(display_cov_missing_message) {
        .log(paste(
            "Note: Some coverage files were not set or not found.\n\n",
            "To set coverage files, use `covfile(se) <- filenames`.",
            "Setting COV files is only supported if every COV file for",
            "every sample in the NxtSE object is valid.\n"
        ), "message")
    }

    return(se)
}

.makeSE_adjust_paths <- function(se) {
    assays <- se@assays
    path <- se@metadata[["sourcePath"]]
    
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        a <- .makeSE_expand_assay_path(getListElement(assays, i), path)
        assays <- setListElement(assays, i, a)
    }
    return(assays)
}

.makeSE_expand_assay_path <- function(assay, path) {
    DelayedArray::modify_seeds(assay,
        function(x) {
            x@filepath <- file.path(path, x@filepath)
            x
        }
    )
}

# Iterates through introns; removes overlapping minor introns
.makeSE_iterate_IR <- function(se, verbose = TRUE) {

    junc_PSI <- junc_PSI(se)
    
    se.IR <- se[rowData(se)$EventType == "IR", , drop = FALSE]

    row_to_junc <- match(rowData(se.IR)$EventRegion, rownames(junc_PSI))
    names(row_to_junc) <- rowData(se.IR)$EventRegion
    row_to_junc <- row_to_junc[!is.na(row_to_junc)]
    
    se.coords <- names(row_to_junc)
    se.coords.gr <- junc_gr(se)[row_to_junc]
    names(se.coords.gr) <- names(row_to_junc)
    
    mcols(se.coords.gr) <- data.frame(
        means = .makeSE_progress_rowMeans(
            junc_PSI[names(se.coords.gr),,drop = FALSE])
    )
    
    if (length(se.coords.gr) > 0) {
        if(verbose) .log(paste("Iterating through IR events to determine introns",
            "of main isoforms"), type = "message")
        include <- .makeSE_iterate_IR_select_events_fast(se.coords.gr)
        se.coords.final <- se.coords.gr[include]
        se.coords.excluded <- se.coords.gr[!include]

        # Iteration to find events not overlapping with se.IR.final
        include <- .makeSE_iterate_IR_retrieve_excluded_introns(
            se.coords.final, se.coords.excluded)
        iteration <- 0
        while (length(include) > 0 & length(se.coords.final) > 0) {
            iteration <- iteration + 1
            if(verbose) .log(paste("Iteration", iteration), type = "message")
            dash_progress(paste("Iteration", iteration), 8)
            se.coords.excluded <- se.coords.excluded[include]

            include <- .makeSE_iterate_IR_select_events_fast(se.coords.excluded)

            if (length(include) > 0 && !all(include)) {
                se.coords.final <- c(se.coords.final,
                    se.coords.excluded[include])
                se.coords.excluded <-
                    se.coords.excluded[!include]
                include <- .makeSE_iterate_IR_retrieve_excluded_introns(
                    se.coords.final, se.coords.excluded)
            } else if (length(include) > 0) {
                se.coords.final <- c(se.coords.final,
                    se.coords.excluded)
                include <- c()
            } else {
                include <- c()
            }
        }

        se <- se[c(
            which(rowData(se.IR)$EventRegion %in% names(se.coords.final)),
            which(rowData(se)$EventType != "IR")
        ), ]
        rm(se.coords.final, se.coords.excluded, include)
    }
    
    rm(se.coords, se.coords.gr)
    gc()
    return(se)
}

.makeSE_progress_rowMeans <- function(junc_PSI) {
    pb <- progress::progress_bar$new(
        format = " calculating junction means [:bar] :percent eta: :eta",
        total = ncol(junc_PSI), clear = FALSE, width= 100)
    junc_means <- rep(0, nrow(junc_PSI))
    for(i in seq_len(ncol(junc_PSI))) {
        junc_means <- junc_means + junc_PSI[,i]
        pb$tick()
    }
    junc_means <- junc_means / ncol(junc_PSI)
    junc_means
}

# Selects introns of major isoforms
.makeSE_iterate_IR_select_events <- function(se.coords.gr, junc_PSI) {
    if(length(se.coords.gr) == 0) return(logical(0))
    if(length(se.coords.gr) == 1) return(TRUE)
    
    gr <- se.coords.gr
    gr.reduced <- reduce(gr)

    OL <- findOverlaps(gr, gr.reduced)
    junc_PSI.group <- as.data.table(
        junc_PSI[names(se.coords.gr), , drop = FALSE])
    junc_PSI.group$means <- rowMeans(junc_PSI.group)
    junc_PSI.group$group <- to(OL)
    junc_PSI.group[, c("max_means") := max(get("means")),
        by = "group"]
        
    res <- junc_PSI.group$means == junc_PSI.group$max_means
    
    rm(junc_PSI.group)
    gc()
    return(res)
}

# Selects introns of major isoforms
.makeSE_iterate_IR_select_events_fast <- function(gr) {
    if(length(gr) == 0) return(logical(0))
    if(length(gr) == 1) return(TRUE)

    OL <- findOverlaps(gr, reduce(gr))
    
    junc_PSI.group <- data.table(
        means = mcols(gr)$means,
        group = to(OL)
    )
    junc_PSI.group[, c("max_means") := max(get("means")), by = "group"]
        
    res <- junc_PSI.group$means == junc_PSI.group$max_means
    
    rm(junc_PSI.group)
    gc()
    return(res)
}

# Find excluded introns that does not overlap with given selection of introns
.makeSE_iterate_IR_retrieve_excluded_introns <- function(
        se.coords.final, se.coords.excluded) {
        
    if (length(se.coords.excluded) > 0) {
        if(length(se.coords.final) == 0) {
            return(rep(TRUE, length(se.coords.excluded)))
        }
        final.gr <- se.coords.final
        excluded.gr <- se.coords.excluded

        OL <- findOverlaps(excluded.gr, final.gr)
        include <- which(!(
            seq_len(length(excluded.gr))) %in% sort(unique(from(OL))))
    } else {
        include <- c()
    }
    return(include)
}

# NEW METHOD
# Iterates through introns; removes overlapping minor introns
.makeSE_iterate_IR_new <- function(se, verbose = TRUE) {

    # junc_PSI <- junc_PSI(se)
    # junc_counts <- junc_counts(se)
    junc_counts <- junc_counts(se)
    
    se.IR <- se[rowData(se)$EventType == "IR", , drop = FALSE]
    if(nrow(se.IR) == 0) return(se)

    # Match IR events to junction counts in se
    row_to_junc <- match(rowData(se.IR)$EventRegion, rownames(junc_counts))
    names(row_to_junc) <- rowData(se.IR)$EventRegion

    # Remove IR events not found in junction counts
    se.IR <- se.IR[!is.na(row_to_junc),]
    row_to_junc <- row_to_junc[!is.na(row_to_junc)]
    
    # Quick way to get GRanges of IR events in se
    se.coords.gr <- junc_gr(se)[row_to_junc]
    names(se.coords.gr) <- names(row_to_junc)
    
    if (length(se.coords.gr) > 0) {
        # Overlap between IR events and repertoire of junction counts
        OL <- findOverlaps(
            se.coords.gr, # coord2GR(EventRegion), 
            junc_gr(se) # coord2GR(rownames(junc_counts))
        )

        # Hash table of overlaps
        OL_DT <- data.table(
          from = from(OL), to = to(OL)
        )

        # Create unique list of junction counts (to minimize rowMeans calcs)
        uniq_to <- sort(unique(to(OL)))
        matMeans <- rowMeans(junc_counts[uniq_to, ])

        # Match junction means to hash table
        OL_DT$to_idx <- match(OL_DT$to, uniq_to)
        OL_DT$countMeans <- matMeans[OL_DT$to_idx]

        # Reduce hash table by max junction count and its coordinates, per IR
        OL_DT <- OL_DT[OL_DT[, .I[which.max(get("countMeans"))], by="from"]$V1]
        OL_DT$coord <- rownames(junc_counts)[OL_DT$to]

        # Reduce IR events to whether max junction coord matches IR event
        coord_maxSplice <- character(length(se.coords.gr))
        coord_maxSplice[OL_DT$from] <- OL_DT$coord
        se.IR <- se.IR[coord_maxSplice == rowData(se.IR)$EventRegion,]
        
        se <- se[c(
            which(rowData(se)$EventName %in% rowData(se.IR)$EventName),
            which(rowData(se)$EventType != "IR")
        ), ]
        
        rm(OL_DT, coord_maxSplice)
    }
    
    rm(se.IR)
    gc()
    return(se)
}

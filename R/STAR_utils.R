# Creates a temporary FASTA file from locally-stored TwoBit
.STAR_get_FASTA <- function(reference_path) {
    genome.fa <- file.path(reference_path, "resource", "genome.fa")
    if (!file.exists(genome.fa)) {
        genome.fa <- paste0(genome.fa, ".temp")
        genome.2bit <- file.path(reference_path, "resource", "genome.2bit")
        if (!file.exists(genome.2bit)) {
            .log(paste(genome.2bit, "not found"))
        }
        .log("Extracting temp genome FASTA from TwoBit file", "message")
        tmp <- rtracklayer::import(TwoBitFile(genome.2bit))
        rtracklayer::export(
            tmp,
            genome.fa, "fasta"
        )
        rm(tmp)
        gc()
    }
    return(genome.fa)
}

# Creates a temporary unzipped GTF for STAR
.STAR_get_GTF <- function(reference_path) {
    transcripts.gtf <- file.path(reference_path, "resource", "transcripts.gtf")
    if (!file.exists(transcripts.gtf)) {
        if (!file.exists(paste0(transcripts.gtf, ".gz"))) {
            .log(paste(paste0(transcripts.gtf, ".gz"), "not found"))
        }
        .log("Extracting temp Annotation GTF from GZ file", "message")
        R.utils::gunzip(paste0(transcripts.gtf, ".gz"), remove = FALSE,
            overwrite = TRUE)
        file.rename(transcripts.gtf, paste0(transcripts.gtf, ".temp"))
        transcripts.gtf <- paste0(transcripts.gtf, ".temp")
    }
    return(transcripts.gtf)
}

.STAR_clean_temp_FASTA_GTF <- function(reference_path) {
    .log("Cleaning temp genome / gene annotation files", "message")
    genome.fa <- file.path(reference_path, "resource", "genome.fa.temp")
    transcripts.gtf <- file.path(reference_path,
        "resource", "transcripts.gtf.temp")
    if (file.exists(genome.fa)) file.remove(genome.fa)
    if (file.exists(transcripts.gtf)) file.remove(transcripts.gtf)
}

.validate_STAR_fastq_samples <- function(fastq_1, fastq_2) {
    if (!is_valid(fastq_2)) {
        # assume single
        if (!all(file.exists(fastq_1))) {
            .log("Some fastq files were not found")
        }
    } else {
        if (length(fastq_2) != length(fastq_1)) {
            .log(paste("There must be equal numbers of",
                "forward and reverse fastq samples"))
        }
        if (!all(file.exists(fastq_1)) || !all(file.exists(fastq_2))) {
            .log("Some fastq files were not found")
        }
    }
    paired <- (length(fastq_1) == length(fastq_2))
    gzipped <- all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if (!gzipped &&
        (
            any(grepl(paste0("\\", ".gz", "$"), fastq_1)) ||
            (paired && any(grepl(paste0("\\", ".gz", "$"), fastq_2)))
        )
    ) {
        .log(paste("A mixture of gzipped and uncompressed",
            "fastq files found.", "You must supply either all",
            "gzipped or all uncompressed fastq files"))
    }
}

.validate_STAR_trim_sequence <- function(sequence) {
    if (length(sequence) != 1) {
        .log("Multiple adaptor sequences are not supported")
    }
    tryCatch({
        ACGT_sum <- sum(Biostrings::letterFrequency(
            Biostrings::DNAString(sequence),
            letters = "AGCT", OR = 0))
    }, error = function(e) ACGT_sum <- 0)
    if (nchar(sequence) != ACGT_sum) {
        .log("Adaptor sequence can only contain A, C, G or T")
    }
}

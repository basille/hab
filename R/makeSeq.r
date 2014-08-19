## makeSeq
##
##' Define sequences of a trajectory based on the chronology of steps,
##' allowing for gaps.
##'
##' @title Define sequences of a trajectory
##' @param date A \code{POSIXt} object, chronologically ordered (by
##' individal if \code{id} is provided).
##' @param id A vector giving the id of each step (must be the same length
##' as \code{date}). If not provided, the function assumes all steps come
##' from a single individual.
##' @param gap The maximum time interval between two successives steps
##' before starting a new sequence.
##' @param units The unit of \code{gap} (default is \code{hour}).
##' @return For \code{makeSeq}, a numeric vector, with the sequence number
##' for each \code{date}, in the form of an integer series for each
##' individual.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## Dummy data:
##' (steps <- data.frame(date = Sys.time() + c(1:10, 12, 14, 20:25) *
##'     3600, id = c(rep("toto", 4), rep("tata", 9), rep("toto",
##'     5))))
##'
##' ## 1-hour gap, all steps considered from a single individual:
##' makeSeq(steps$date)
##' ## 2-hour gap, all steps considered from a single individual:
##' makeSeq(steps$date, gap = 2)
##' ## 1-hour gap, individual id considered:
##' steps$seq <- makeSeq(steps$date, steps$id)
##' steps
##'
##' ## Note that `seq` uses a series starting from 1 for each individual. Uses
##' ## something like this to have unique seq IDs among individuals:
##' steps$seqid <- paste(steps$id, steps$seq, sep = "-")
##' steps
##'
##' ## Summary of the sequence, individual id considered:
##' summarySeq(steps$seq, steps$id)
makeSeq <- function(date, id = NULL, gap = 1, units = c("hours", "mins", "secs", "days", "weeks"))
{
    ## Check that the date is as POSIXt
    if (!inherits(date, "POSIXt"))
        stop("`date` should be of class `POSIXct`")
    ## Match `units`
    units <- match.arg(units)
    ## If `id` is not provided, `date` is considered as coming from one
    ## single individual
    if (is.null(id)) {
        ## Compute the time intervals between locs, in unit `units`
        ddiff <- difftime(date[-1], date[-length(date)], units = units)
        ## Check that the `date`` is correctly ordered
        if (any(ddiff < 0))
            stop("The `date` must be chronologically ordered.")
        ## The sequence starts with 1, and is the cumulative sum of the
        ## intervals greater than `gap`
        seq <- cumsum(c(1, ddiff > gap))
        ## Intervals between each sequence are given by `ddiff` greater
        ## than `gap`, and is stored in the attribute `inter`
        attr(seq, "inter") <- ddiff[ddiff > gap]
    }
    ## `date` from several individuals
    else {
        ## Check the length of `id`
        if (length(date) != length(id))
            stop("`id` must have the same length as `date`.")
        ## Prepare the global sequence filled with NAs
        seq <- rep(NA, length(date))
        ## Prepare the global interval as a list
        seqinter <- list()
        ## For each `id`, runs the same as above
        for (i in (unique(id))) {
            ## Compute the time intervals between locs, in unit `units`
            ddiff <- difftime(date[id == i][-1], date[id ==
                i][-length(date[id == i])], units = units)
            ## Check that the `date`` is correctly ordered
            if (any(ddiff < 0))
                stop("The `date` must be chronologically ordered by individual.")
            ## The sequence starts with 1, and is the cumulative sum of the
            ## intervals greater than `gap`
            seq[id == i] <- cumsum(c(1, ddiff > gap))
            ## Intervals between each sequence are given by `ddiff` greater
            ## than `gap`
            seqinter[[i]] <- ddiff[ddiff > gap]
        }
        ## Store the intervals in the list
        attr(seq, "inter") <- seqinter
    }
    ## The gap is stored in the attribute `gap`
    attr(seq, "gap") <- gap
    ## The unit for gap is stored in the attribute `gap-unit`
    attr(seq, "gap-unit") <- units
    return(seq)
}

## summarySeq
##
##' @rdname makeSeq
##' @param seq An vector giving a sequence in a trajectory, as given by
##' \code{makeSeq}, or as a series of integers.
##' @return For \code{summarySeq}, an object of class \code{summarySeq}.
##' @export
summarySeq <- function(seq, id = NULL)
{
    ## If `id` is not provided, `seq` is considered as coming from one
    ## single individual
    if (is.null(id))
        ## `summ` is a list giving the number of sequences (max seq #), the
        ## length of each sequence and the intervals between sequences
        summ <- list(nr = max(seq), length = unclass(table(seq)),
            inter = attr(seq, "inter"))
    ## If `id` is provided, seq` from several individuals
    else {
        ## `summ` is the same as above, but constrained by `id`
        summ <- lapply(unique(id), function(i) list(nr = max(seq[id ==
            i]), length = unclass(table(seq[id == i])), inter = attr(seq,
            "inter")[[i]]))
        ## Names of `summ` given by `id`
        names(summ) <- unique(id)
    }
    ## `units` is stored in an attribute
    attr(summ, "units") <- attr(seq, "gap-unit")
    ## `id` is stored in an attribute (NULL if not provided)
    attr(summ, "id") <- id
    ## `summ` is of class `summarySeq`
    class(summ) <- "summarySeq"
    return(summ)
}

## print.summarySeq
##
##' @rdname makeSeq
##' @param x An object of class \code{summarySeq}.
##' @export
print.summarySeq <- function(x, ...)
{
    ## Check that `x` is of class `x`
    if (!inherits(x, "summarySeq"))
        stop("x should be of class `summarySeq`")
    ## If no `id` was provided, just one individual
    if (is.null(attr(x, "id"))) {
        ## Print the number of sequences
        cat("Number of sequences:", x$nr, "\n")
        ## Print a summary of sequence lengths
        cat("\nLength of sequences:\n")
        print(summary(x$length))
        ## Print a summary of sequence intervals
        cat("\nTime intervals between sequences (",
            attr(x, "units"), "):\n", sep = "")
        print(summary(x$inter))
    }
    ## If `id` was provided, several individuals
    else {
        ## Print the total number of sequences
        cat("Total number of sequences:", sum(unlist(lapply(x,
            function(i) i$nr))), "\n")
        ## Print the number of sequences per individual (mean, sd)
        cat("\nNumber of sequences per individual (mean = ",
            round(mean(unlist(lapply(x, function(i) i$nr))), 2),
            ", sd = ", round(sd(unlist(lapply(x, function(i) i$nr))), 2),
            "):\n", sep = "")
        print(unlist(lapply(x, function(i) i$nr)))
        ## Print a summary of sequence lengths (all indidividuals combined)
        cat("\nLength of sequences (all individuals):\n")
        ## Print a summary of sequence intervals (all indidividuals
        ## combined)
        print(summary(unlist(lapply(x, function(i) i$length))))
        cat("\nTime intervals between sequences (all individuals, ",
            attr(x, "units"), "):\n", sep = "")
        print(summary(unlist(lapply(x, function(i) i$inter))))
    }
}

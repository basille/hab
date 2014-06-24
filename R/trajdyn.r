## trajdyn
##
##' Modified version of \code{\link[adehabitatLT]{trajdyn}}, which
##' allows to 1) increment by \code{by} points/steps, 2) only display
##' \code{k} previous points/steps, 3) show the current step as a
##' vector, 4) to modify point and line display, 5) show the current
##' burst/loc/(infolocs) and 6) to update interactively a new or
##' existing variable in \code{infolocs}.
##'
##' \code{y} selects the number of previous points/steps to increment
##' at each step. It defaults to \code{1}, i.e. an increment of 1
##' point/step. Choosing anything different than a positive number
##' sets it back to \code{1}.
##'
##' \code{k} selects the number of previous points/steps to
##' display. It defaults to \code{Inf}, i.e. all
##' points/steps. Choosing anything different than a positive number
##' sets it back to \code{Inf}.
##'
##' \code{v} shows (or removes) the current step as a vector. Note
##' that the vector connects the current location to the next
##' available location (even if there are NAs in the data set).
##'
##' \code{s} shows the current buffer and localisation numbers,
##' together with the associated infolocs (if it exists). If a
##' \code{nvar} is requested, only this variable is shown.
##'
##' The argument \code{nvar} allows to work on a given variable: if a
##' variable of that name already exists in \code{infolocs(x)}, the
##' values of the variable are retrieved from there; otherwise, a
##' variable filled with NAs is used. If a \code{nvar} is requested,
##' \code{d} "deletes" the value and resets it to \code{NA}; every
##' other letter not in use in any option is saved in \code{nvar}. The
##' letters available are: "c", "e", "f", "h", "j", "m", "t", "u",
##' "w", "x".
##'
##' On a QWERTY keyboard:
##'
##' we t u
##'   f hj
##'  xcm
##'
##' On a AZERTY keyboard:
##'
##'  e t u
##'   f hj m
##' wxc
##'
##' It is possible to use point and line parameters globally for every
##' trajectory displayed. In this case, \code{ppar} and \code{lpar} need
##' just be a list of graphical parameters, such as \code{list(pch = 21,
##' col = "black", bg = "white")}. It is also possible to use parameters
##' for single steps, using as graphical parameter a list of vectors of
##' length equal to each trajectory. Such information can be based on
##' \code{infolocs}, see Example.
##'
##' @title Interactive Display of Objects of Class \code{ltraj}
##' @param na.rm Logical, whether to remove missing locations.
##' @param addvec Numeric, whether to hihglight the current location (1,
##' default), the current step (2) or nothing (0).
##' @param by The number of previous points/steps to increment at each
##' step. Default is an increment of 1 point/step.
##' @param only The number of previous points/steps to
##' display. Default is \code{Inf}, i.e. all points/steps.
##' @param ppar A list of arguments that allows the user to modify point
##' display, using any argument available to \code{points}. Default is
##' \code{list(pch = 16)}. See Details.
##' @param lpar A list of arguments that allows the user to modify line
##' display, using any argument available to \code{lines}. Default is
##' \code{list(lwd = 2)}. See Details.
##' @param nvar A character string giving the name of a variable.
##' @return If a \code{nvar} is provided, return the original ltraj
##' with updated values in \code{infolocs(nvar)}.
##' @author Modified by Mathieu Basille
##' \email{basille@@ase-research.org}
##' @export
##' @examples
##' \dontrun{
##' data(puechcirc)
##' ##'
##' ## Use of `by` and `only` to select the previous k points/steps:
##' trajdyn(puechcirc, by = 10, only = 20)
##' ##'
##' ## Use of `ppar` and `lpar` globally:
##' trajdyn(puechcirc, ppar = list(col = "red"), lpar = list(col = "blue"))
##' ##'
##' ## Create some random `infolocs`:
##' info <- list(data.frame(col = sample(c("red", "grey"),
##'          80, rep = TRUE), stringsAsFactors = FALSE),
##'      data.frame(col = sample(c("blue", "darkred"),
##'          69, rep = TRUE), stringsAsFactors = FALSE),
##'      data.frame(col = sample(c("darkgreen", "purple"),
##'          66, rep = TRUE), stringsAsFactors = FALSE))
##' ## Watch the row names:
##' info <- mapply(function(x, y) {
##'     row.names(x) <- row.names(y)
##'     return(x)
##' }, info, puechcirc, SIMPLIFY = FALSE)
##' infolocs(puechcirc) <- info
##' ##'
##' ## Use the infolocs to color points and steps:
##' trajdyn(puechcirc, by = 1, only = 20, ppar = list(pch = 19,
##'     col = infolocs(puechcirc, "col", simplify = TRUE)),
##'     lpar = list(col = infolocs(puechcirc, "col", simplify = TRUE)))
##' ##'
##' ## The same without removing the missing locations:
##' trajdyn(puechcirc, by = 1, only = 20, ppar = list(pch = 19,
##'     col = infolocs(puechcirc, "col", simplify = TRUE)),
##'     lpar = list(col = infolocs(puechcirc, "col", simplify = TRUE)),
##'     na.rm = FALSE)
##' ##'
##' ## Use of `nvar` to dynamically fill in new data:
##' (newtraj <- trajdyn(puechcirc, nvar = "Var"))
##' }
trajdyn <- function (x, burst = attr(x[[1]], "burst"), na.rm = TRUE, hscale = 1, vscale = 1, addvec = 1,
    by = 1, only = Inf, recycle = TRUE, ppar = list(pch = 16), lpar = list(lwd = 2),
    nvar = NULL, display = c("guess", "windows", "tk"), ...)
{
    if (!inherits(x, "ltraj"))
        stop("x should be of class 'ltraj'")
    e1 <- new.env(parent = baseenv())
    ## Are there any lists in point/line parameters?
    plist <- sapply(ppar, is.list)
    llist <- sapply(lpar, is.list)
    ## Check the length of list parameters
    if (any(plist)) {
        for (k in (1:length(ppar))[plist])
            if (!isTRUE(all.equal(unlist(lapply(ppar[[k]], length)), unlist(lapply(x, nrow)), check.attributes = FALSE)))
                stop("Point parameters for individual locations must have the same length as the corresponding burst")
    }
    if (any(llist)) {
        for (k in (1:length(lpar))[llist])
            if (!isTRUE(all.equal(unlist(lapply(lpar[[k]], length)), unlist(lapply(x, nrow)), check.attributes = FALSE)))
                stop("Line parameters for individual steps must have the same length as the corresponding burst")
    }
    ## End of modification
    ## Remove missing values, and recompute trajectory parameters
    ## typeII <- attr(x, "typeII")
    ## x <- lapply(x, function(i) {
    ##     jj <- i[!is.na(i$x), ]
    ##     attr(jj, "id") <- attr(i, "id")
    ##     attr(jj, "burst") <- attr(i, "burst")
    ##     return(jj)
    ## })
    ## class(x) <- c("ltraj", "list")
    ## attr(x, "typeII") <- typeII
    ## attr(x, "regular") <- is.regular(x)
    if (na.rm) {
        ## Remove NAs from individual point/line parameters
        nas <- lapply(x, function(i) !is.na(i$x))
        names(nas) <- id(x)
        ## Only if the list of parameter is of length > 0
        if (length(ppar) > 0 & any(plist))
            for (k in (1:length(ppar))[plist])
                ppar[[k]] <- mapply(function(x, y) {
                  x[y]
                }, ppar[[k]], nas, SIMPLIFY = FALSE)
        if (length(lpar) > 0 & any(llist))
            for (k in (1:length(lpar))[llist])
                lpar[[k]] <- mapply(function(x, y) {
                  x[y]
                }, lpar[[k]], nas, SIMPLIFY = FALSE)
        x <- na.omit(x)
    }
    ## Store point/line parameters in e1 for first burst
    ppark <- ppar
    if (any(plist)) {
        for (k in (1:length(ppark))[plist])
            ppark[k] <- ppark[[k]][burst]
    }
    assign("ppark", ppark, envir = e1)
    lpark <- lpar
    if (any(llist)) {
        for (k in (1:length(lpark))[llist])
            lpark[k] <- lpark[[k]][burst]
    }
    assign("lpark", lpark, envir = e1)
    ## End of modification
    u <- x
    assign("x", x[burst = burst], envir = e1)
    assign("v", x[burst = burst], envir = e1)
    ## Prepare 'info': get the whole infolocs of the ltraj
    info <- infolocs(x)
    ## If a nvar is requested
    if (!is.null(nvar)) {
        ## Save the ltraj for backup and prepare a list in e1
        ltr_bkp <- x
        ## If info is NULL, set 'status_info' to "noinfo" (defaults to
        ## 'nvar_exists')
        status_info <- "nvar_exists"
        if (is.null(info))
            status_info <- "noinfo"
        ## Check that 'nvar' is a character
        if (!inherits(nvar, "character"))
            stop("nvar should be a character string")
        ## Get the column nvar in infolocs
        info <- infolocs(x, nvar)
        ## If info is NULL, create an empty info list and set
        ## 'status_info' to "novar"
        if (is.null(info)) {
            status_info <- "novar"
            info <- lapply(x, function(le) setNames(data.frame(rep(NA,
                nrow(le)), row.names = row.names(le)), nvar))
        }
    }
    ## Store info in e1
    assign("info", info, envir = e1)
    ## End of modification
    assign("ajouli", FALSE, envir = e1)
    assign("ajoupo", FALSE, envir = e1)
    assign("ajoubu", FALSE, envir = e1)
    assign("addpoints", TRUE, envir = e1)
    assign("addlines", TRUE, envir = e1)
    ## Prepare the step vector status (0: none, 1: point, 2: vector)
    assign("addvec", addvec, envir = e1)
    ## End of modification
    assign("lim", TRUE, envir = e1)
    assign("buadd", burst, envir = e1)
    assign("K", 1, envir = e1)
    assign("N", nrow(get("x", envir = e1)[[1]]), envir = e1)
    assign("cusr", rep(0 + NA, 4), envir = e1)
    assign("cplt", rep(0 + NA, 4), envir = e1)
    ## Assign 'by' in e1
    assign("by", by, envir = e1)
    ## Assign 'only' in e1
    assign("only", only, envir = e1)
    opt <- options(warn = -1)
    on.exit(options(opt))
    dsp <- substring(match.arg(display), 1, 1)
    if (dsp == "g")
        dsp <- switch(.Platform$OS.type, windows = "w", "t")
    if (dsp == "t" && !require(tkrplot))
        stop("'tkrplot' package needed\n")
    if (dsp == "t")
        assign("hoho", 1, envir = e1)
    replot <- function() {
        opar <- par(mar = c(0, 0, 0, 0), bg = "white")
        ## Retrieve point/line parameters and keep only last 'only'
        pparky <- get("ppark", envir = e1)
        if (any(plist)) {
            for (k in (1:length(pparky))[plist])
                pparky[[k]] <- pparky[[k]][tail(1:get("K",
                  envir = e1), ifelse(is.null(get("only", envir = e1)),
                  get("K", envir = e1), get("only", envir = e1)))]
        }
        lparky <- get("lpark", envir = e1)
        if (any(llist)) {
            for (k in (1:length(lparky))[llist])
                lparky[[k]] <- lparky[[k]][tail(1:get("K",
                  envir = e1), ifelse(is.null(get("only", envir = e1)),
                  get("K", envir = e1), get("only", envir = e1)))]
        }
        ## End of modification
        tmptmp <- get("x", envir = e1)
        attr(tmptmp[[1]], "id") <- " "
        assign("x", tmptmp, envir = e1)
        if (get("lim", envir = e1)) {
            ## Allows for NAs
            assign("xlim", range(get("x", envir = e1)[[1]]$x, na.rm = TRUE),
                envir = e1)
            ## Allows for NAs
            assign("ylim", range(get("x", envir = e1)[[1]]$y, na.rm = TRUE),
                envir = e1)
        }
        plot(get("x", envir = e1), id = attr(get("x", envir = e1)[[1]],
            "id"), addlines = FALSE, addp = FALSE, final = FALSE,
            xlim = get("xlim", envir = e1), ylim = get("ylim",
                envir = e1), ...)
        assign("cusr", par("usr"), envir = e1)
        assign("cplt", par("plt"), envir = e1)
        scatterutil.sub(as.character(get("x", envir = e1)[[1]]$date[get("K",
            envir = e1)]), 1, "topleft")
        if (get("ajoubu", envir = e1)) {
            lapply(u[burst = get("buadd", envir = e1)], function(zz) {
                if (get("addpoints", envir = e1))
                  points(zz[, c("x", "y")], pch = 16, col = "grey")
                if (get("addlines", envir = e1))
                  lines(zz[, c("x", "y")], pch = 16, col = "grey")
            })
        }
        ## addlines before addpoints
        if (get("addlines", envir = e1)) {
            ## if (get("K", envir = e1) > 1) {
                ## Plot only the last 'only' steps, and allows for line
                ## modification:
                ## lines(get("x", envir = e1)[[1]][1:get("K", envir = e1),
                ##   c("x", "y")], lwd = 2)
                bla <- get("x", envir = e1)[[1]][tail(1:get("K",
                  envir = e1), ifelse(is.null(get("only", envir = e1)),
                  get("K", envir = e1), get("only", envir = e1))),
                  c("x", "y", "dx", "dy")]
                do.call(segments, c(list(x0 = bla$x, y0 = bla$y, x1 =
                  bla$x + bla$dx, y1 = bla$y + bla$dy), lparky))
                ## do.call(lines, c(get("x", envir = e1)[[1]][tail(1:get("K",
                ##   envir = e1), ifelse(is.null(get("only", envir = e1)),
                ##   get("K", envir = e1), get("only", envir = e1))),
                ##   c("x", "y")], lparky))
            }
        if (get("addpoints", envir = e1))
            ## Plot only the last 'only' points, and allows for point
            ## modification:
            ## points(get("x", envir = e1)[[1]][1:get("K", envir = e1),
            ##     c("x", "y")], pch = 16)
            do.call(points, c(get("x", envir = e1)[[1]][tail(1:get("K",
                envir = e1), ifelse(is.null(get("only", envir = e1)),
                get("K", envir = e1), get("only", envir = e1))),
                c("x", "y")], pparky))
        if (get("ajouli", envir = e1))
            lines(c(get("a1", envir = e1)[1], get("a2", envir = e1)[1]),
                c(get("a1", envir = e1)[2], get("a2", envir = e1)[2]),
                lwd = 2, col = "red")
        if (get("ajoupo", envir = e1))
            points(get("a5", envir = e1)[1], get("a5", envir = e1)[2],
                pch = 16, col = "red", cex = 1.7)
        iti <- unlist(get("x", envir = e1)[[1]][get("K", envir = e1),
            c("x", "y")])
        ## Add the step point or vector
        ## points(iti[1], iti[2], col = "blue", pch = 16, cex = 1.4)
        vec <- get("addvec", envir = e1)
        if (vec == 1) {
            current <- get("x", envir = e1)[[1]][get("K", envir = e1), ]
            points(current$x, current$y, col = "blue", pch = 16, cex = 1.4)
        }
        if (vec == 2) {
            current <- get("x", envir = e1)[[1]][get("K", envir = e1), ]
            points(current$x, current$y, col = "blue", pch = 16, cex = 1.4)
            arrows(current$x, current$y, current$x + current$dx, current$y + current$dy, lwd = 3, length = .1, col = "blue")
        }
        ## End of modification
        par(opar)
    }
    ## Remove the final \n
    help.txt <- paste("\n-------- to obtain this help, type 'h' ------------------",
        "n/p            -- Next/Previous relocation", "a              -- show All relocations",
        "y              -- browse the trajectory bY n steps",
        "g              -- Go to...", "0-9            -- show a given part of the path",
        "k              -- display only the K previous steps",
        "s              -- Show burst/relocation/(infolocs)",
        "b              -- change Burst", "i              -- add/remove other bursts on the graph",
        "z/o            -- Zoom in/Out", "Left-Click     -- measure the distance between two points",
        "Right-Click    -- identify a relocation",
        "r/l/v          -- add or remove points/Lines/Vector",
        "q              -- Quit", "---------------------------------------------------------\n", sep = "\n")
    ## If a nvar is requested, complete the help text
    if (!is.null(nvar))
        help.txt <- paste0(help.txt,
            "d              -- Delete nvar for this relocation\n",
            "other letters  -- Set nvar to the letter value\n",
            "---------------------------------------------------------\n\n")
    assign("D", 0, envir = e1)
    assign("a1", 0, envir = e1)
    assign("a2", 0, envir = e1)
    if (dsp == "t") {
        tt <- tcltk::tktoplevel()
        tcltk::tkwm.title(tt, "Exploration of Animal Movements")
        img <- tkrplot::tkrplot(tt, replot, hscale = hscale,
            vscale = vscale)
        txt <- tcltk::tktext(tt, bg = "white", font = "courier 10")
        scr <- tcltk::tkscrollbar(tt, repeatinterval = 5, command = function(...) tcltk::tkyview(txt,
            ...))
        tcltk::tkconfigure(txt, yscrollcommand = function(...) tcltk::tkset(scr,
            ...))
        tcltk::tkpack(img, side = "top")
        tcltk::tkpack(txt, side = "left", fill = "both", expand = TRUE)
        tcltk::tkpack(scr, side = "right", fill = "y")
        iw <- as.numeric(tcltk::tcl("image", "width", tcltk::tkcget(img,
            "-image")))
        ih <- as.numeric(tcltk::tcl("image", "height", tcltk::tkcget(img,
            "-image")))
    }
    showz <- function() switch(dsp, w = replot(), t = {
        tkrplot::tkrreplot(img)
    })
    type <- function(s) switch(dsp, w = cat(s), t = {
        tcltk::tkinsert(txt, "end", s)
        tcltk::tksee(txt, "end")
    })
    type(help.txt)
    cc <- function(x, y) {
        if (dsp == "t") {
            x <- (as.double(x) - 1)/iw
            y <- 1 - (as.double(y) - 1)/ih
        }
        px <- (x - get("cplt", envir = e1)[1])/(get("cplt", envir = e1)[2] -
            get("cplt", envir = e1)[1])
        py <- (y - get("cplt", envir = e1)[3])/(get("cplt", envir = e1)[4] -
            get("cplt", envir = e1)[3])
        ux <- px * (get("cusr", envir = e1)[2] - get("cusr",
            envir = e1)[1]) + get("cusr", envir = e1)[1]
        uy <- py * (get("cusr", envir = e1)[4] - get("cusr",
            envir = e1)[3]) + get("cusr", envir = e1)[3]
        c(ux, uy)
    }
    mm.w <- function(buttons, x, y) {
        if (buttons == 0) {
            i <- get("D", envir = e1)
            if (i == 0) {
                assign("a1", cc(x, y), envir = e1)
                assign("D", 1, envir = e1)
            }
            if (i == 1) {
                assign("a2", cc(x, y), envir = e1)
                assign("D", 0, envir = e1)
                di <- sqrt(sum((get("a2", envir = e1) - get("a1",
                  envir = e1))^2))
                cat(paste("distance:", round(di, 6), "\n"))
                lines(c(get("a1", envir = e1)[1], get("a2", envir = e1)[1]),
                  c(get("a1", envir = e1)[2], get("a2", envir = e1)[2]),
                  lwd = 2, col = "red")
            }
            return()
        }
        if (buttons == 2) {
            w <- get("v", envir = e1)[[1]][1:get("K", envir = e1),
                ]
            assign("a3", cc(x, y), envir = e1)
            di <- sqrt((w$x - get("a3", envir = e1)[1])^2 + (w$y -
                get("a3", envir = e1)[2])^2)
            print(w[which.min(di), ])
            cat("\n")
            points(w[which.min(di), c("x", "y")], pch = 16, col = "red",
                cex = 1.7)
            return()
        }
    }
    mm.t <- function(x, y) {
        i <- get("D", envir = e1)
        if (i == 0) {
            assign("a1", cc(x, y), envir = e1)
            assign("D", 1, envir = e1)
        }
        if (i == 1) {
            assign("a2", cc(x, y), envir = e1)
            assign("D", 0, envir = e1)
            di <- sqrt(sum((get("a2", envir = e1) - get("a1",
                envir = e1))^2))
            type(paste("distance:", di, "\n"))
            assign("ajouli", TRUE, envir = e1)
            showz()
            assign("ajouli", FALSE, envir = e1)
        }
        return()
    }
    mm.t2 <- function(x, y) {
        w <- get("v", envir = e1)[[1]][1:get("K", envir = e1),
            ]
        assign("a3", cc(x, y), envir = e1)
        di <- sqrt((w$x - get("a3", envir = e1)[1])^2 + (w$y -
            get("a3", envir = e1)[2])^2)
        assign("a5", unlist(w[which.min(di), c("x", "y")]), envir = e1)
        assign("ajoupo", TRUE, envir = e1)
        showz()
        assign("ajoupo", FALSE, envir = e1)
        tmp <- w[which.min(di), ]
        se <- unlist(lapply((max(nchar(names(tmp)) + nchar(sapply(tmp,
            as.character)) + 1) - nchar(names(tmp)) - nchar(sapply(tmp,
            as.character))), function(zz) paste(rep(" ", zz),
            collapse = "")))
        so <- unlist(lapply(1:length(tmp), function(i) paste(paste(names(tmp)[i],
            as.character(tmp[1, i]), sep = se[i]), "\n")))
        type(paste("Relocation", row.names(w)[which.min(di)],
            ":\n"))
        sapply(so, type)
        type("\n")
        return()
    }
    mm.mouse <- function(buttons, x, y) {
        assign("a8", cc(x, y), envir = e1)
        return()
    }
    mm.mouset <- function(x, y) {
        assign("a8", cc(x, y), envir = e1)
        return()
    }
    kb <- function(A) {
        key <- tolower(A)
        if (key == "q") {
            if (dsp == "t")
                tcltk::tkdestroy(tt)
            return("OK - Finished")
        }
        if (key %in% c(0:9)) {
            if (key > 0)
                assign("K", round(seq(1, get("N", envir = e1),
                  length = 11))[as.numeric(key) + 1], envir = e1)
            if (key == 0)
                assign("K", 1, envir = e1)
            showz()
        }
        if (key == "z") {
            assign("tmppx", (get("cusr", envir = e1)[1:2] - get("cusr",
                envir = e1)[1])/2, envir = e1)
            assign("xlim", c((get("a8", envir = e1)[1] - (get("tmppx",
                envir = e1)[2] - get("tmppx", envir = e1)[1])/2),
                (get("a8", envir = e1)[1] + (get("tmppx", envir = e1)[2] -
                  get("tmppx", envir = e1)[1])/2)), envir = e1)
            assign("tmppy", (get("cusr", envir = e1)[3:4] - get("cusr",
                envir = e1)[3])/2, envir = e1)
            assign("ylim", c((get("a8", envir = e1)[2] - (get("tmppy",
                envir = e1)[2] - get("tmppy", envir = e1)[1])/2),
                (get("a8", envir = e1)[2] + (get("tmppy", envir = e1)[2] -
                  get("tmppy", envir = e1)[1])/2)), envir = e1)
            assign("lim", FALSE, envir = e1)
            showz()
        }
        if (key == "o") {
            assign("lim", TRUE, envir = e1)
            showz()
        }
        if (key == "n") {
            if (get("K", envir = e1) <= get("N", envir = e1))
                ## Browse by increments of 'by'
                ## assign("K", get("K", envir = e1) + 1, envir = e1)
                assign("K", get("K", envir = e1) + get("by",
                  envir = e1), envir = e1)
            if (get("K", envir = e1) > get("N", envir = e1)) {
                if (recycle)
                  assign("K", 1, envir = e1)
                if (!recycle) {
                  assign("K", get("N", envir = e1), envir = e1)
                  cat("End of burst !\n")
                }
            }
            showz()
        }
        if (key == "l") {
            assign("addlines", !get("addlines", envir = e1),
                envir = e1)
            showz()
        }
        if (key == "g") {
            if (dsp == "w") {
                recom <- TRUE
                while (recom) {
                  rr <- readline("Enter a relocation number: ")
                  recom <- FALSE
                  if (!(rr %in% row.names(get("x", envir = e1)[[1]]))) {
                    cat("invalid number\n")
                    recom <- TRUE
                  }
                }
                assign("K", which(row.names(get("x", envir = e1)[[1]]) ==
                  as.numeric(rr)), envir = e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(row.names(get("x", envir = e1)[[1]])[1])
                tu <- tcltk::tktoplevel(tt, width = 500, height = 50)
                tcltk::tkwm.title(tu, "Enter a relocation number")
                tcltk::tkwm.resizable(tu, 0, 0)
                en <- tcltk::tkentry(tu, textvariable = lv, width = 50)
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    rr <- tcltk::tclvalue(lv)
                    if (!(rr %in% row.names(get("x", envir = e1)[[1]]))) {
                      tcltk::tkmessageBox(message = "invalid number",
                        type = "ok")
                    }
                    else {
                      assign("K", which(row.names(get("x", envir = e1)[[1]]) ==
                        as.numeric(rr)), envir = e1)
                      showz()
                      tcltk::tkdestroy(tu)
                    }
                  })
                tcltk::tkpack(en, side = "top", fill = "both")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
            }
        }
        if (key == "r") {
            assign("addpoints", !get("addpoints", envir = e1),
                envir = e1)
            showz()
        }
        ## Display the step vector
        if (key == "v") {
            assign("addvec", ifelse(get("addvec", envir = e1) == 2, 0, get("addvec", envir = e1) + 1), envir = e1)
            showz()
        }
        if (key == "b") {
            assign("K", 1, envir = e1)
            if (dsp == "w") {
                assign("hoho", select.list(unlist(lapply(u, function(y) attr(y,
                  "burst")))), envir = e1)
                type(paste("Choice of the burst:", get("hoho",
                  envir = e1), "\n\n"))
                assign("x", u[burst = get("hoho", envir = e1)],
                  envir = e1)
                assign("v", u[burst = get("hoho", envir = e1)],
                  envir = e1)
                assign("N", nrow(get("x", envir = e1)[[1]]),
                  envir = e1)
                ## Retrieve point/line parameters for this burst and store
                ## it in e1
                ppark <- ppar
                if (any(plist)) {
                    for (k in (1:length(ppark))[plist])
                        ppark[k] <- ppark[[k]][get("hoho", envir = e1)]
                }
                assign("ppark", ppark, envir = e1)
                lpark <- lpar
                if (any(llist)) {
                    for (k in (1:length(lpark))[llist])
                        lpark[k] <- lpark[[k]][get("hoho", envir = e1)]
                }
                assign("lpark", lpark, envir = e1)
                ## End of modification
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(unlist(lapply(u, function(y) attr(y,
                  "burst"))))
                bubu <- unlist(lapply(u, function(y) attr(y,
                  "burst")))
                tu <- tcltk::tktoplevel(tt)
                tcltk::tkwm.title(tu, "Choose a burst of relocations")
                tcltk::tkwm.resizable(tu, 0, 0)
                tfr <- tcltk::tkframe(tu)
                tli <- tcltk::tklistbox(tfr, bg = "white", font = "courier 12",
                  listvariable = lv)
                scr2 <- tcltk::tkscrollbar(tfr, repeatinterval = 5,
                  command = function(...) tcltk::tkyview(tli,
                    ...))
                tcltk::tkconfigure(tli, yscrollcommand = function(...) tcltk::tkset(scr2,
                  ...))
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    assign("hoho", ifelse(nchar(tcltk::tclvalue(tcltk::tkcurselection(tli))) ==
                      0, 1, as.numeric(tcltk::tclvalue(tcltk::tkcurselection(tli))) +
                      1), envir = e1)
                    type(paste("Choice of the burst:", bubu[get("hoho",
                      envir = e1)], "\n\n"))
                    tcltk::tkdestroy(tu)
                  })
                tcltk::tkpack(tli, side = "left", fill = "both",
                  expand = TRUE)
                tcltk::tkpack(scr2, side = "right", fill = "y")
                tcltk::tkpack(tfr, side = "right", fill = "y")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
                assign("x", u[burst = bubu[get("hoho", envir = e1)]],
                  envir = e1)
                assign("v", u[burst = bubu[get("hoho", envir = e1)]],
                  envir = e1)
                assign("N", nrow(get("x", envir = e1)[[1]]),
                  envir = e1)
                ## Retrieve point/line parameters for this burst and store
                ## it in e1
                ppark <- ppar
                if (any(plist)) {
                    for (k in (1:length(ppark))[plist])
                        ppark[k] <- ppark[[k]][get("hoho", envir = e1)]
                }
                assign("ppark", ppark, envir = e1)
                lpark <- lpar
                if (any(llist)) {
                    for (k in (1:length(lpark))[llist])
                        lpark[k] <- lpark[[k]][get("hoho", envir = e1)]
                }
                assign("lpark", lpark, envir = e1)
                ## End of modification
                showz()
            }
        }
        if (key == "i") {
            if (get("ajoubu", envir = e1)) {
                assign("ajoubu", FALSE, envir = e1)
                showz()
            }
            else {
                if (dsp == "w") {
                  assign("buadd", select.list(unlist(lapply(u,
                    function(y) attr(y, "burst"))), multiple = TRUE),
                    envir = e1)
                  if (length(get("buadd", envir = e1) > 0)) {
                    type(paste("show bursts:", paste(get("buadd",
                      envir = e1), collapse = " "), "\n\n"))
                    assign("ajoubu", TRUE, envir = e1)
                    showz()
                  }
                }
                if (dsp == "t") {
                  lv <- tcltk::tclVar(unlist(lapply(u, function(y) attr(y,
                    "burst"))))
                  bubu <- unlist(lapply(u, function(y) attr(y,
                    "burst")))
                  tu <- tcltk::tktoplevel(tt)
                  tcltk::tkwm.title(tu, "Choose one or several bursts")
                  tcltk::tkwm.resizable(tu, 0, 0)
                  tfr <- tcltk::tkframe(tu)
                  tli <- tcltk::tklistbox(tfr, bg = "white",
                    font = "courier 12", listvariable = lv, selectmode = "multiple")
                  scr2 <- tcltk::tkscrollbar(tfr, repeatinterval = 5,
                    command = function(...) tcltk::tkyview(tli,
                      ...))
                  tcltk::tkconfigure(tli, yscrollcommand = function(...) tcltk::tkset(scr2,
                    ...))
                  submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                    command = function() {
                      argg <- ifelse(nchar(tcltk::tclvalue(tcltk::tkcurselection(tli))) ==
                        0, 1, 0)
                      if (argg == 0) {
                        assign("ajoubu", TRUE, envir = e1)
                        assign("buadd", bubu[as.numeric(unlist(strsplit(tcltk::tclvalue(tcltk::tkcurselection(tli)),
                          " "))) + 1], envir = e1)
                        type(paste("show bursts:", paste(get("buadd",
                          envir = e1), collapse = " "), "\n\n"))
                        showz()
                        tcltk::tkdestroy(tu)
                      }
                    })
                  tcltk::tkpack(tli, side = "left", fill = "both",
                    expand = TRUE)
                  tcltk::tkpack(scr2, side = "right", fill = "y")
                  tcltk::tkpack(tfr, side = "right", fill = "y")
                  tcltk::tkpack(submit.but, side = "bottom")
                  tcltk::tkwait.window(tu)
                  assign("x", u[burst = bubu[get("hoho", envir = e1)]],
                    envir = e1)
                  assign("v", u[burst = bubu[get("hoho", envir = e1)]],
                    envir = e1)
                  assign("N", nrow(get("x", envir = e1)[[1]]),
                    envir = e1)
                  showz()
                }
            }
        }
        if (key == "p") {
            if (get("K", envir = e1) > 1)
                ## Browse by increments of 'by'
                ## assign("K", get("K", envir = e1) - 1, envir = e1)
                assign("K", get("K", envir = e1) - get("by",
                  envir = e1), envir = e1)
            ## Condition becomes K <= 1
            ## if (get("K", envir = e1) == 1) {
            if (get("K", envir = e1) <= 1) {
                if (recycle)
                  assign("K", get("N", envir = e1), envir = e1)
                if (!recycle) {
                  assign("K", 1, envir = e1)
                  cat("Beginning of burst!\n")
                }
            }
            showz()
        }
        if (key == "a") {
            assign("K", get("N", envir = e1), envir = e1)
            showz()
        }
        if (key == "h")
            type(help.txt)
        ## If 'y', change 'by'
        if (key == "y") {
            if (dsp == "w") {
                rr <- as.numeric(readline("Increment by n steps: "))
                assign("by", ifelse(is.na(rr) | rr < 0, 1, rr),
                  envir = e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(ifelse(is.null(get("by",
                  envir = e1)), "NULL", as.character(get("by",
                  envir = e1))))
                tu <- tcltk::tktoplevel(tt, width = 500, height = 50)
                tcltk::tkwm.title(tu, "Increment by n steps")
                tcltk::tkwm.resizable(tu, 0, 0)
                en <- tcltk::tkentry(tu, textvariable = lv, width = 50)
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    rr <- as.numeric(tcltk::tclvalue(lv))
                    assign("by", ifelse(is.na(rr) | rr < 0, 1,
                      rr), envir = e1)
                    showz()
                    tcltk::tkdestroy(tu)
                  })
                tcltk::tkpack(en, side = "top", fill = "both")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
            }
        }
        ## If 'k', change 'only'
        if (key == "k") {
            if (dsp == "w") {
                rr <- as.numeric(readline("Display only k steps: "))
                assign("only", ifelse(is.na(rr) | rr < 0, Inf,
                    rr), envir = e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(ifelse(is.null(get("only",
                  envir = e1)), "NULL", as.character(get("only",
                  envir = e1))))
                tu <- tcltk::tktoplevel(tt, width = 500, height = 50)
                tcltk::tkwm.title(tu, "Display only k steps")
                tcltk::tkwm.resizable(tu, 0, 0)
                en <- tcltk::tkentry(tu, textvariable = lv, width = 50)
                submit.but <- tcltk::tkbutton(tu, text = "    OK     ",
                  command = function() {
                    rr <- as.numeric(tcltk::tclvalue(lv))
                    assign("only", ifelse(is.na(rr) | rr < 0,
                      Inf, rr), envir = e1)
                    showz()
                    tcltk::tkdestroy(tu)
                  })
                tcltk::tkpack(en, side = "top", fill = "both")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
            }
        }
        ## If 's', show burst/loc/infolocs
        ## Function 'showloc', which will be reused later...
        showloc <- function() {
            bubu <- unlist(lapply(u, function(y) attr(y, "burst")))
            text <- paste0("Burst: ", bubu[get("hoho", envir = e1)], "; Loc: ",
                row.names(get("x", envir = e1)[[1]])[get("K", envir = e1)])
            ## If there is no infolocs, just show burst/loc
            if (is.null(info))
                type(paste0(text, "\n"))
            ## Else show info
            else {
                ## If no nvar is requested, all variables shown
                if (is.null(nvar)) {
                  infk <- get("info", envir = e1)[[get("hoho",
                    envir = e1)]][which(row.names(get("info",
                    envir = e1)[[get("hoho", envir = e1)]]) ==
                    row.names(get("x", envir = e1)[[1]])[get("K", envir = e1)]),
                      , drop = FALSE]
                  infk <- apply(infk, 2, function(x) ifelse(is.numeric(x),
                    round(x, 3), x))
                  inftext <- paste(names(infk), as.character(infk), sep = ": ",
                    collapse = "; ")
                }
                ## If a nvar is requested, only show nvar
                else {
                  infv <- get("info", envir = e1)[[get("hoho",
                    envir = e1)]][which(row.names(get("info",
                    envir = e1)[[get("hoho", envir = e1)]]) ==
                    row.names(get("x", envir = e1)[[1]])[get("K", envir = e1)]),
                    nvar]
                  infv <- ifelse(is.numeric(infv), round(infv, 3), infv)
                  inftext <- paste0(nvar, ": ", infv)
                }
                type(paste0(text, "; ", inftext, "\n"))
            }
        }
        if (key == "s") {
            showloc()
        }
        ## If a nvar is requested, all letters not attributed can be
        ## used to save in the new variable
        if (!is.null(nvar) & key %in% c("c", "d", "e", "f", "h",
            "j", "m", "t", "u", "w", "x")) {
            ## Get back the list prepared
            infonew <- get("info", envir = e1)
            ## If D, delete the attribute and set it back to NA (note
            ## that the nvar take initial NAs into account)
            if (key == "d")
                infonew[[get("hoho", envir = e1)]][which(row.names(infonew[[get("hoho",
                  envir = e1)]]) == row.names(get("x", envir = e1)[[1]])[get("K",
                  envir = e1)]), nvar] <- NA
            ## Otherwise, set it to the value of the key
            else infonew[[get("hoho", envir = e1)]][which(row.names(infonew[[get("hoho",
                envir = e1)]]) == row.names(get("x", envir = e1)[[1]])[get("K",
                envir = e1)]), nvar] <- key
            ## Update lvar in e1
            assign("info", infonew, envir = e1)
            ## Displays the burst/loc/nvar
            showloc()
        }
        return()
    }
    showz()
    toto <- switch(dsp, w = getGraphicsEvent("", onKeybd = kb,
        onMouseDown = mm.w, onMouseMove = mm.mouse), t = {
        tcltk::tkbind(tt, "<Key>", kb)
        tcltk::tkbind(img, "<Button-1>", mm.t)
        tcltk::tkbind(img, "<Motion>", mm.mouset)
        tcltk::tkbind(img, "<Button-3>", mm.t2)
        tcltk::tkwait.window(tt)
    })
    ## If a nvar is requested, return the ltraj with modified infolocs
    if (!is.null(nvar)) {
        ## If 'noinfo', use 'info' directly
        if (status_info == "noinfo")
            infolocs(ltr_bkp) <- get("info", envir = e1)
        ## If 'novar', simply merge 'info' with the existing
        ## 'infolocs'
        else if (status_info == "novar")
            infolocs(ltr_bkp) <- mapply(cbind, infolocs(ltr_bkp),
                get("info", envir = e1), SIMPLIFY = FALSE)
        ## If 'nvar' already exists in 'infolocs', replace the column
        ## 'nvar' by the new one in 'info'
        else infolocs(ltr_bkp) <- mapply(function(x, y) {
            x[[nvar]] <- y[[nvar]]
            return(x)
        }, infolocs(ltr_bkp), get("info", envir = e1), SIMPLIFY = FALSE)
        return(ltr_bkp)
    }
}

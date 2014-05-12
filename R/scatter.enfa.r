## scatter.enfa
##
scatter.enfa <- function(x, xax = 1, yax = 2, pts = FALSE,
    nc = TRUE, grid = TRUE, percent = 95, clabel = 1, side = c("top",
        "bottom", "axes", "none"), posieig = c("none", "top",
        "bottom"), Adensity, Udensity, Aangle, Uangle, Aborder,
    Uborder, Acol, Ucol, Alty, Ulty, Abg, Ubg, Ainch, Uinch,
    ...)
{
    side <- match.arg(side)
    if (!inherits(x, "enfa"))
        stop("Object of class 'enfa' expected")
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    x1 <- x$li[, xax]
    x1 <- c(x1 - diff(range(x1)/50), x1 + diff(range(x1))/50)
    xlim <- range(x1)
    y1 <- x$li[, yax]
    y1 <- c(y1 - diff(range(y1)/50), y1 + diff(range(y1))/50)
    ylim <- range(y1)
    pmar <- t(x$mar * x$cw) %*% as.matrix(x$co[, c(xax, yax)])
    scatterutil.base(dfxy = x$li[, c(xax, yax)], xax = 1, yax = 2,
        xlim = xlim, ylim = ylim, grid = grid, addaxes = FALSE,
        cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "",
        csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL,
        area = NULL, add.plot = FALSE)
    if (pts) {
        if (missing(Acol))
            Acol <- gray(0.8)
        if (missing(Ucol))
            Ucol <- "black"
        if (missing(Abg))
            Abg <- gray(0.8)
        if (missing(Ubg))
            Ubg <- "black"
        if (missing(Ainch))
            Ainch <- 0.03
        if (missing(Uinch))
            Uinch <- Ainch * max(x$pr)
        symbols(x$li[, c(xax, yax)], circles = rep(1, length(x$pr)),
            fg = Acol, bg = Abg, inches = Ainch, add = TRUE)
        symbols(x$li[x$pr > 0, c(xax, yax)], circles = x$pr[x$pr >
            0], fg = Ucol, bg = Ubg, inches = Uinch, add = TRUE)
        abline(v = 0)
        abline(h = 0)
        if (nc)
            symbols(pmar, circles = 1, fg = "black", bg = "white",
                inches = Ainch * 2, add = TRUE)
    }
    else {
        if (missing(Adensity))
            Adensity <- NULL
        if (missing(Udensity))
            Udensity <- NULL
        if (missing(Aangle))
            Aangle <- 45
        if (missing(Uangle))
            Uangle <- 45
        if (missing(Aborder))
            Aborder <- NULL
        if (missing(Uborder))
            Uborder <- NULL
        if (missing(Acol))
            Acol <- gray(0.95)
        if (missing(Ucol))
            Ucol <- gray(0.6)
        if (missing(Alty))
            Alty <- NULL
        if (missing(Ulty))
            Ulty <- NULL
        pcff <- function(xy) {
            mo <- apply(xy, 2, mean)
            dis <- apply(xy, 1, function(x) sum((x - mo)^2))
            xy <- xy[dis < quantile(dis, percent/100), ]
            return(xy[chull(xy[, 1], xy[, 2]), ])
        }
        mcpA <- pcff(x$li[, c(xax, yax)])
        mcpU <- pcff(x$li[rep(1:length(x$pr), x$pr), c(xax, yax)])
        polygon(mcpA, density = Adensity, angle = Aangle, border = Aborder,
            col = Acol, lty = Alty)
        polygon(mcpU, density = Udensity, angle = Uangle, border = Uborder,
            col = Ucol, lty = Ulty)
        abline(v = 0)
        abline(h = 0)
        if (nc)
            points(pmar, pch = 21, bg = "white", cex = 1.5)
    }
    dfarr <- x$co[, c(xax, yax)]
    born <- par("usr")
    k1 <- min(dfarr[, 1])/born[1]
    k2 <- max(dfarr[, 1])/born[2]
    k3 <- min(dfarr[, 2])/born[3]
    k4 <- max(dfarr[, 2])/born[4]
    k <- c(k1, k2, k3, k4)
    dfarr <- 0.75 * dfarr/max(k)
    s.arrow(dfarr, clabel = clabel, addaxes = FALSE, add.plot = TRUE)
    if (side != "none") {
        if (xax == 1)
            xleg <- "mar"
        else xleg <- paste("sp", xax - 1)
        if (yax == 1)
            yleg <- "mar"
        else yleg <- paste("sp", yax - 1)
        xl <- par("usr")[1]
        xr <- par("usr")[2]
        yd <- par("usr")[3]
        yu <- par("usr")[4]
        if (side == "top") {
            tra <- paste(" xax =", xleg, "\n yax =", yleg)
            wd <- strwidth(tra, cex = 1)
            ht <- strheight(tra, cex = 1) * 1.5
            rect(xl, yu - ht, xl + wd, yu, col = "white", border = 0)
            text(xl + wd/2, yu - ht/2, tra, cex = 1)
        }
        if (side == "bottom") {
            tra <- paste(" xax =", xleg, "\n yax =", yleg)
            wd <- strwidth(tra, cex = 1)
            ht <- strheight(tra, cex = 1) * 1.5
            rect(xl, yd + ht, xl + wd, yd, col = "white", border = 0)
            text(xl + wd/2, yd + ht/2, tra, cex = 1)
        }
        if (side == "axes") {
            trax <- paste(" xax =", xleg)
            wdx <- strwidth(trax, cex = 1)
            htx <- strheight(trax, cex = 1) * 1.5
            rect(xr - 1.05 * wdx, 0 + 0.05 * htx, xr - 0.05 *
                wdx, 0 + 1.05 * htx, col = "white", border = 0)
            text(xr - 1.05 * wdx/2, 0 + htx/2, trax, cex = 1)
            tray <- paste(" yax =", yleg)
            wdy <- strwidth(tray, cex = 1)
            hty <- strheight(tray, cex = 1) * 1.5
            rect(0 + 0.05 * wdy, yu - 1.05 * hty, 0 + 1.05 *
                wdy, yu - 0.05 * hty, col = "white", border = 0)
            text(0 + wdy/2, yu - 1.05 * hty/2, tray, cex = 1)
        }
    }
    add.scatter.eig(x$s, x$nf, xax = xax - 1, yax = yax - 1,
        posi = posieig, sub = "Sp. eigenvalues")
    box()
}

## library(adehabitatHS)
## data(lynxjura)
## map <- lynxjura$map
## locs <- lynxjura$locs
## locs <- locs[slot(locs, "data")[,2]!="D",]
## slot(map,"data")[,4] <- sqrt(slot(map,"data")[,4])
## tab <- slot(map, "data")
## pr <- slot(count.points(locs, map), "data")[,1]
## pc <- dudi.pca(tab, scannf = FALSE)
## (enfa1 <- enfa(pc, pr, scannf = FALSE, nf = 3))
##
## scatter(enfa1)
## scatter(enfa1, grid = FALSE, side = "axes", posieig = "bottom")
## scatter(enfa1, grid = FALSE, side = "axes", posieig = "top", xax = 2, yax = 3)

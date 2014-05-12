## rdsteps
##
### Mathieu Basille, basille@ase-research.org
### Last modified: 2012-11-19
### 'rdSteps' computes random steps. It takes a data frame as input,
### mandatory fields are 'id', 'x/ystart', 'rel./abs.angle', and
### 'dist'.
###   'distMax' allows to use a maximum distance in the empirical
###     distributions of distances used for random steps.
###   'simult' allows to draw simultaneously step lenght and turning
###     angle from the same step (instead of two independent
###     distributions)
###   'other' allows to draw random steps from distributions of all
###     other individuals.
### If they don't exist, two new columns 'xend' and 'yend' are created
### with the end of the random step coordinates computed based on
### 'dx/dy'.
### A new column 'case' is created (or overwritten) which takes the
### value 1 for observed steps and 0 for random ones.
### The 'protect' argument is to copy the value of the set of
### 'protected' variables from the observed step to the random
### ones. For example, use 'protect = c("Area", "Sex")' to copy the
### value of Area and Sex.
##
## Working version without 'simult'
rdSteps <- function(df, emp = df, nr = 10, distMax = Inf, other = TRUE,
    id = "id", xstart = "x", ystart = "y", date = "date", dx = "dx",
    dy = "dy", dist = "dist", dt = "dt", abs.angle = "abs.angle",
    rel.angle = "rel.angle", xend = "xend", yend = "yend", case = "case",
    protect = NULL, reproducible = TRUE)
{
    if (!exists(xend, df))
        df[, xend] <- df[, xstart] + df[, dx]
    if (!exists(yend, df))
        df[, yend] <- df[, ystart] + df[, dy]
    df[, case] <- 1
    angles <- na.omit(emp[, rel.angle])
    idA <- emp[!is.na(emp[, rel.angle]), id]
    dists <- na.omit(emp[emp[, dist] <= distMax, dist])
    idD <- emp[!is.na(emp[, dist]) & emp[, dist] <= distMax,
        id]
    rdStepsId <- function(ldfk) {
        idk <- as.character(ldfk[1, id])
        if (other)
            return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
                ], anglesk = angles[idA != idk], distsk = dists[idD !=
                idk], i = i))))
        else return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
            ], anglesk = angles, distsk = dists, i = i))))
    }
    rdStep <- function(pt, anglesk = anglesk, distsk = distsk,
        i = i) {
        if (is.na(pt[, xstart]) | is.na(pt[, rel.angle]))
            return()
        else {
            if (reproducible)
                set.seed(i)
            rhord <- sample(anglesk, nr, replace = TRUE)
            alphard <- pt[, abs.angle] - pt[, rel.angle] + rhord
            if (reproducible)
                set.seed(i)
            distrd <- sample(distsk, nr, replace = TRUE)
            rd <- pt
            rd[1, ] <- NA
            rd[, id] <- pt[1, id]
            rd <- rd[rep(1, nr), ]
            rd[, xstart] <- pt[1, xstart]
            rd[, ystart] <- pt[1, ystart]
            rd[, date] <- pt[1, date]
            rd[, dx] <- cos(alphard * 180/pi) * distrd
            rd[, dy] <- sin(alphard * 180/pi) * distrd
            rd[, dist] <- distrd
            rd[, dt] <- pt[1, dt]
            rd[, abs.angle] <- alphard
            rd[, rel.angle] <- rhord
            rd[, xend] <- rd[, xstart] + rd[, dx]
            rd[, yend] <- rd[, ystart] + rd[, dy]
            rd[, case] <- 0
            if (!is.null(protect))
                rd[, protect] <- pt[1, protect]
            return(rbind(pt, rd))
        }
    }
    ldf <- split(df, f = df[, id])
    return(do.call("rbind", lapply(ldf, rdStepsId)))
}

rdSteps <- function(df, emp = df, nr = 10, distMax = Inf, simult = FALSE,
    other = TRUE, id = "id", xstart = "x", ystart = "y", date = "date",
    dx = "dx", dy = "dy", dist = "dist", dt = "dt", abs.angle = "abs.angle",
    rel.angle = "rel.angle", xend = "xend", yend = "yend", case = "case",
    protect = NULL, reproducible = TRUE) {
    if (!exists(xend, df))
        df[, xend] <- df[, xstart] + df[, dx]
    if (!exists(yend, df))
        df[, yend] <- df[, ystart] + df[, dy]
    df[, id] <- factor(df[, id])
    emp[, id] <- factor(emp[, id])
    df[, case] <- 1
    if (simult) {
        tmp <- na.omit(emp[, c(dist, rel.angle, id)])
        dists <- tmp[, dist]
        angles <- tmp[, rel.angle]
        idA <- idD <- tmp[, id]
    }
    else {
        angles <- na.omit(emp[, rel.angle])
        idA <- emp[!is.na(emp[, rel.angle]), id]
        dists <- na.omit(emp[emp[, dist] <= distMax, dist])
        idD <- emp[!is.na(emp[, dist]) & emp[, dist] <= distMax,
            id]
    }
    rdStepsId <- function(ldfk) {
        idk <- as.character(ldfk[1, id])
        if (other)
            return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
                ], anglesk = angles[idA != idk], distsk = dists[idD !=
                idk], i = i))))
        else return(do.call("rbind", lapply(1:nrow(ldfk), function(i) rdStep(ldfk[i,
            ], anglesk = angles, distsk = dists, i = i))))
    }
    rdStep <- function(pt, anglesk = anglesk, distsk = distsk,
        i = i) {
        if (is.na(pt[, xstart]) | is.na(pt[, rel.angle]))
            return()
        else {
            if (simult) {
                if (reproducible)
                  set.seed(i)
                spl <- sample(1:length(anglesk), nr, replace = TRUE)
                rhord <- anglesk[spl]
                alphard <- pt[, abs.angle] - pt[, rel.angle] +
                  rhord
                distrd <- distsk[spl]
            }
            else {
                if (reproducible)
                  set.seed(i)
                rhord <- sample(anglesk, nr, replace = TRUE)
                alphard <- pt[, abs.angle] - pt[, rel.angle] +
                  rhord
                if (reproducible)
                  set.seed(i)
                distrd <- sample(distsk, nr, replace = TRUE)
            }
            rd <- pt
            rd[1, ] <- NA
            rd[, id] <- pt[1, id]
            rd <- rd[rep(1, nr), ]
            rd[, xstart] <- pt[1, xstart]
            rd[, ystart] <- pt[1, ystart]
            rd[, date] <- pt[1, date]
            rd[, dx] <- cos(alphard * 180/pi) * distrd
            rd[, dy] <- sin(alphard * 180/pi) * distrd
            rd[, dist] <- distrd
            rd[, dt] <- pt[1, dt]
            rd[, abs.angle] <- alphard
            rd[, rel.angle] <- rhord
            rd[, xend] <- rd[, xstart] + rd[, dx]
            rd[, yend] <- rd[, ystart] + rd[, dy]
            rd[, case] <- 0
            if (!is.null(protect))
                rd[, protect] <- pt[1, protect]
            return(rbind(pt, rd))
        }
    }
    ldf <- split(df, f = df[, id])
    return(do.call("rbind", lapply(ldf, rdStepsId)))
}

## data(puechcirc)
## puechcirc ## class ltraj
## uu <- ld(puechcirc)
## bli <- rdSteps(uu, simul = TRUE)
##
## load("~/Travail/Data/Québec/Côte-Nord/Loup/Outputs/LoupLocsCN.RData")
## head(LoupLocsCN)
## bla <- LoupLocsCN[1:100, ]
## bla$Id <- factor(bla$Id)
## bli <- rdSteps(bla, id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Sector", "Pack", "Id", "Sex"))
## blu <- rdSteps(bla, emp = LoupLocsCN, id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Sector", "Pack", "Id", "Sex"))
## rdSteps(bla[1, ], id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend")
## rdSteps(bla[1, ], emp = LoupLocsCN, id = "Id", xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Schedule", "Burst"))
## rdSteps(bla[2, ], id = "Id", other = FALSE, xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend")
## rdSteps(bla[2, ], emp = LoupLocsCN, id = "Id", xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Schedule", "Burst"))
## debug(rdSteps)
## debug(rdStepsId)
## debug(rdStep)
## bli <- rdSteps(LoupLocsCN, id = "Id", xstart = "X", ystart = "Y", date = "Date", dx = "dX", dy = "dY", dist = "Dist", dt = "dt", abs.angle = "Abs.angle", rel.angle = "Rel.angle", xend = "Xend", yend = "Yend", protect = c("Sector", "Pack", "Id", "Sex"))

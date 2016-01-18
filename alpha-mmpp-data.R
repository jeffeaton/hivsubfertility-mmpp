################################################################################
####                                                                        ####
####  This script creates individual episode data for MMPP model from       ####
####  processed ALPHA data.                                                 ####
####                                                                        ####
####  Created by Jeff Eaton on 16 January 2016.                             ####
####                                                                        ####
################################################################################


setwd("~/Documents/Research/hiv-subfertility/mmpp/")

load("alpha-fertility-data.RData")

library(survival)
library(plyr)


############################
####  Prepare episodes  ####
############################

## Restrict residences

res <- subset(res, sex == "Female")


## Add HIV data

lastneg <- aggregate(testdate ~ site+idno, subset(hiv, !result), max)
firstpos <- aggregate(testdate ~ site+idno, subset(hiv, result), min)
names(lastneg)[3] <- "lastneg"
names(firstpos)[3] <- "firstpos"

res <- merge(res, firstpos, all.x=TRUE)
res <- merge(res, lastneg, all.x=TRUE)

res$lastneg[is.na(res$lastneg)] <- -Inf
res$firstpos[is.na(res$firstpos)] <- Inf


## records with no HIV data 
with(res, table(site, firstpos==Inf & lastneg==-Inf))

## inconsistent HIV data
with(res, table(site, lastneg >= firstpos))

res <- subset(res, lastneg < firstpos)

## HIV positive before age 15
nrow(subset(res, firstpos <= dob + 15))
res <- subset(res, firstpos > dob + 15)

res$res_id <- seq_len(nrow(res))

## Keep only residences starting before 2005
res <- subset(res, entry < 2005)



## Link fertility events to episodes
fert$fert_id <- seq_len(nrow(fert))
fert <- merge(fert, res, all.x=TRUE)
with(fert, table(deliverydate >= entry & deliverydate < exit))
fert <- subset(fert, deliverydate >= entry & deliverydate < exit)
sum(duplicated(fert$fert_id))  # confirm no duplicated fertility events


e.expose <- with(res,
                 data.frame(site = site,
                            id   = res_id,
                            dob  = dob,
                            time = dob,
                            age  = 0,
                            type = "expose",
                            hivst = "HIV-",
                            fertelig = FALSE))

e.entry <- with(res,
                data.frame(site = site,
                           id   = res_id,
                           dob  = dob,
                           time = entry,
                           age  = entry - dob,
                           type = "entry",
                           hivst = ifelse(entry < lastneg, "HIV-",
                                   ifelse(entry >= firstpos, "HIV+",
                                          "Unknown")),
                           fertelig = FALSE))

e.exit <- with(res,
               data.frame(site = site,
                          id   = res_id,
                          dob  = dob,
                          time = exit,
                          age  = exit - dob,
                          type = ifelse(death, "death", "exit"),
                          hivst = ifelse(exit < lastneg, "HIV-",
                                  ifelse(exit > firstpos, "HIV+",
                                         "Unknown")),
                          fertelig = TRUE))

e.lastneg <- with(subset(res, lastneg > -Inf & lastneg < exit),
                  data.frame(site = site,
                             id   = res_id,
                             dob  = dob,
                             time = lastneg,
                             age  = lastneg - dob,
                             type = "lastneg",
                             hivst = "HIV-",
                             fertelig = ifelse(lastneg > entry, TRUE, FALSE)))

e.firstpos <- with(subset(res, firstpos < exit),
                   data.frame(site = site,
                              id   = res_id,
                              dob  = dob,
                              time = firstpos,
                              age  = firstpos - dob,
                              type = "firstpos",
                              hivst = "Unknown",
                              fertelig = ifelse(firstpos > entry, TRUE, FALSE)))

e.fert <- with(fert,
               data.frame(site = site,
                          id   = res_id,
                          dob  = dob,
                          time = deliverydate,
                          age  = deliverydate - dob,
                          type = "birth",
                          hivst = ifelse(deliverydate < lastneg, "HIV-",
                                  ifelse(deliverydate > firstpos, "HIV+",
                                         "Unknown")),
                          fertelig = TRUE))

events <- rbind(e.expose, e.entry, e.exit, e.firstpos, e.lastneg, e.fert)
events <- events[order(events$site, events$id, events$time),]


system.time(epis <- split(events, events[c("site", "id")], drop=TRUE))
system.time(epis <- lapply(epis, function(x){val <- x[-1,]; val$dur <- x$time[-1] - x$time[-nrow(x)]; val}))
system.time(epis <- rbind.fill(epis))
epis$start <- epis$time - epis$dur
epis$dur <- NULL

## Keep only time <2005
epis$type <- as.character(epis$type)
e <- survSplit(epis, 2005, "time", "type", "start", episode="period")
e <- subset(e, period == 0)


## Cut into age groups

e$age <- e$time - e$dob
e$startage <- e$start - e$dob
  
e <- survSplit(e, 3:10*5, "age", "type", "startage", episode="agegr")
e$time <- e$dob + e$age
e$start <- e$dob + e$startage
e <- subset(e, agegr %in% 1:7)
e$type[e$type == 0] <- "cens"
e$dur <- e$time - e$start

epis <- e


## Construct Q matrix index
Qlist.mm <- expand.grid(fertelig=c(FALSE, TRUE), agegr=1:7)  # model matrix for Q objects
Qlist.mm$Qidx <- 1:nrow(Qlist.mm)

epis <- merge(epis, Qlist.mm)


## Split into datasets
epis.hivp <- subset(epis, hivst != "HIV-")
epis.hivp$type[epis.hivp$type %in% c("entry", "exit")] <- "cens"
epis.hivp <- split(epis.hivp, epis.hivp$id)
epis.hivp <- epis.hivp[order(epis.hivp$site, epis.hivp$id, epis.hivp$startage),]
epis.hivp$last <- FALSE
epis.hivp$last[cumsum(aggregate(type~site+id, epis.hivp, NROW)$type)] <- TRUE  # flag last episode for each individual (relies on order)

epis.hivn <- subset(epis, hivst == "HIV-")
dat.hivn <- Reduce(merge, list(aggregate(dur ~ Qidx, epis.hivn, sum),
                               setNames(aggregate(type=="birth" ~ Qidx, epis.hivn, sum), c("Qidx", "births")),
                               setNames(aggregate(type=="death" ~ Qidx, epis.hivn, sum), c("Qidx", "deaths"))))


events.entry <- rbind(with(subset(epis, type=="entry"),
                           data.frame(site     = site,
                                      age      = age,
                                      type     = "entry",
                                      agegr    = agegr,
                                      fertelig = FALSE)),
                      data.frame(site=  NA,
                                 age   = 3:9*5,
                                 type   = "cens",
                                 agegr  = 1:7,
                                 fertelig = FALSE))

events.entry <- merge(events.entry, Qlist.mm)
events.entry <- events.entry[order(events.entry$age),]

epis.entry <- events.entry[-1,]
epis.entry$dur <- events.entry$age[-1] - events.entry$age[-nrow(events.entry)]
epis.entry$entry <- epis.entry$type == "entry"


save(epis, epis.hivp, epis.hivn, dat.hivn, epis.entry, Qlist.mm, file="alpha-mmpp-data.RData")

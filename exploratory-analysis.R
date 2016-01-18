load("alpha-fertility-data.RData")

library(survival)

## Restrict residences

res <- subset(res, sex == "Female")

## pre-ART only: date < 2005
res <- survSplit(res, c(1989, 2005, 2015), "exit", "death", "entry", episode="period")
res <- subset(res, period %in% 1:2)
res$period <- factor(res$period, 1:2, c("<2005", ">2005"))

## fertile ages: 15-49
res$entryage <- res$entry - res$dob
res$exitage <- res$exit - res$dob

res <- survSplit(res, c(15, 50), "exitage", "death", "entryage")
res$entry <- res$dob + res$entryage
res$exit <- res$dob + res$exitage
res <- subset(res, entryage >= 15 & entryage < 50)



## Add HIV data

firstpos <- aggregate(testdate ~ site+idno, subset(hiv, result), min)
lastneg <- aggregate(testdate ~ site+idno, subset(hiv, !result), max)
names(firstpos)[3] <- "firstpos"
names(lastneg)[3] <- "lastneg"

res <- merge(res, firstpos, all.x=TRUE)
res <- merge(res, lastneg, all.x=TRUE)


table(res$lastneg == res$firstpos)

table(!is.na(res$firstpos))
table(!is.na(res$lastneg))

## records with no HIV data 
with(res, table(site, is.na(firstpos) & is.na(lastneg)))

## inconsistent HIV data
with(res, table(site, lastneg >= firstpos))

res <- subset(res, is.na(lastneg) | is.na(firstpos) | lastneg < firstpos)
               
## Split into age groups
res <- survSplit(res, 3:10*5, "exitage", "death", "entryage", episode="agegr")
res$entry <- res$dob + res$entryage
res$exit <- res$dob + res$exitage
res$agegr <- factor(5*(res$agegr+2))


## Split into HIV status
res.hivn <- subset(res, entry < lastneg)
res.hivp <- subset(res, exit > firstpos)
res.unkn <- subset(res, (is.na(lastneg) | lastneg < exit) & (is.na(firstpos) | firstpos > entry))

res.hivn$hivst <- "HIV-"
res.hivn$death <- with(res.hivn, ifelse(lastneg < exit, FALSE, death))
res.hivn$exit <- with(res.hivn, pmin(lastneg, exit))

res.hivp$hivst <- "HIV+"
res.hivp$entry <- with(res.hivp, pmax(firstpos, entry))

res.unkn$hivst <- "Unknown"
res.unkn$entry <- with(res.unkn, pmax(entry, lastneg, na.rm=TRUE))
res.unkn$death <- with(res.unkn, ifelse(!is.na(firstpos) & firstpos < exit, FALSE, death))
res.unkn$exit <- with(res.unkn, pmin(exit, firstpos, na.rm=TRUE))


## Check that PYs add up
with(res.hivn, sum(exit - entry)) +
with(res.hivp, sum(exit - entry)) + 
with(res.unkn, sum(exit - entry))

with(res, sum(exit-entry))

res <- rbind(res.hivn, res.hivp, res.unkn)

res$entryage <- res$entry - res$dob
res$exitage <- res$exit - res$dob



## Link fertility events to episodes
fert$fert_id <- seq_len(nrow(fert))

fert <- merge(fert, res, all.x=TRUE)

with(fert, table(deliverydate >= entry & deliverydate < exit))

fert <- subset(fert, deliverydate >= entry & deliverydate < exit)

sum(duplicated(fert$fert_id))  # confirm no duplicated fertility events

aggr.births <- with(fert, table(agegr, hivst, site))
aggr.pys <- with(res, tapply(exit - entry, list(agegr, hivst, site), sum))

aggr.births / aggr.pys

fert$births <- 1
res$pys <- res$exit - res$entry

fert.rate <- merge(aggregate(births ~ site + agegr + period + hivst, fert, sum),
                   aggregate(pys ~ site + agegr + period + hivst, res, sum))
fert.rate$rate <- with(fert.rate, births/pys)


####  Is being in an unknown episode that ends in HIV+ associated with fertility?

with(fert, table(hivst, exit == firstpos))
fert$hivst2 <- fert$hivst
fert$hivst2[fert$hivst=="Unknown" & fert$exit == fert$firstpos] <- "pre-positive"

res$hivst2 <- res$hivst
res$hivst2[res$hivst=="Unknown" & res$exit == res$firstpos] <- "pre-positive"

fert.rate2 <- merge(aggregate(births ~ site + agegr + period + hivst2, fert, sum),
                    aggregate(pys ~ site + agegr + period + hivst2, res, sum))


fert.rate2$rate <- with(fert.rate2, births/pys)


tmp <- read.dta("~/Documents/Research/hiv-subfertility/ALPHA/fertility_jeff_stata12.dta")

subset(tmp, study_name == "Rakai")
subset(tmp, study_name == "Kisesa")

save(res, fert, file="fert-episode-data.RData")


############################
####  Prepare episodes  ####
############################

names(fert)

fert$fertelig <- TRUE

entry <- with(subset(res, entryage < 15),
              data.frame(site  = site,
                         idno  = idno,
                         time  = entryage,
                         type  = ifelse(
                         
                         

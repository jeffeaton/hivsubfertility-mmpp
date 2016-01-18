################################################################################
####                                                                        ####
####  This script creates pooled residence, HIV, and fertility data from    ####
####  ALPHA specs 6.1, 6.2, and 7.                                          ####
####                                                                        ####
####  Created by Jeff Eaton on 16 January 2016.                             ####
####                                                                        ####
################################################################################


library(foreign)
library(readstata13)
library(lubridate)

decimal_date_na <- function(x){ y <- rep(NA, length(x));  y[!is.na(x)] <- decimal_date(x[!is.na(x)]); return(y) }

setwd("~/Documents/Research/hiv-subfertility/mmpp/")


kisesa.res <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA6_spec1_Kisesa_clean.dta")
masaka.res <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA6_spec1_Masaka_clean.dta")
rakai.res <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA6_spec1_Rakai_clean.dta")

kisesa.hiv <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA6_spec2b_Kisesa_clean.dta")
rakai.hiv.wide <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA6_spec2b_Rakai_clean_wide.dta")
masaka.hiv.wide <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA6_spec2b_Masaka_clean_wide.dta")

kisesa.fert <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA7_spec2_Kisesa_clean.dta")
masaka.fert <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA7_spec2_Masaka_clean.dta")
rakai.fert.wide <- read.dta13("~/Documents/Data/ALPHA/fertility/v2016-01-14/ALPHA7_spec2_Rakai_clean_wide.dta")


######################
####  Residences  ####
######################

kisesa.res$sex <- factor(kisesa.res$sex, 1:2, c("Male", "Female"))
masaka.res$sex <- factor(as.integer(masaka.res$sex), 1:2, c("Male", "Female"))
rakai.res$sex <- factor(rakai.res$sex, 1:2, c("Male", "Female"))

kisesa.res$dob <- as.Date(kisesa.res$dob, origin="1960-01-01")

exittype.levels <- c("present in study site", "death", "out-migration", "lost to follow-up")
kisesa.res$exit_type <- factor(as.integer(kisesa.res$exit_type), 1:4, exittype.levels)
masaka.res$exit_type <- factor(as.integer(masaka.res$exit_type), 1:4, exittype.levels)
rakai.res$exit_type <- factor(as.integer(rakai.res$exit_type), 1:4, exittype.levels)

res.cols <- c("idno", "sex", "dob", "entry_date", "exit_date", "exit_type", "study_name")
res <- rbind(masaka.res[,res.cols], rakai.res[,res.cols], kisesa.res[,res.cols])

res$site <- factor(res$study_name)

## round date to nearest day
res$entry_date <- as.Date(as.integer(res$entry_date), origin="1970-01-01")
res$exit_date <- as.Date(as.integer(res$exit_date), origin="1970-01-01")

res$dob <- decimal_date_na(res$dob)
res$entry <- decimal_date_na(res$entry_date)
res$exit <- decimal_date_na(res$exit_date)
res$death <- res$exit_type=="death"

res$res_id <- seq_len(nrow(res))
res <- res[c("res_id", "site", "idno", "sex", "dob", "entry", "exit", "death")]


## Remove invalid exit dates
res <- subset(res, !is.na(exit))

## if multiple deaths, keep the latest one
deaths <- subset(res, death)
dup.deaths <- merge(deaths, deaths[(duplicated(deaths[c("site", "idno")])), c("site", "idno")])
dup.deaths

deaths <- deaths[order(deaths$site, deaths$idno, deaths$exit),]
deaths <- subset(deaths, !duplicated(deaths[c("site", "idno")], fromLast=TRUE))

## eliminate any residence episodes that occur after death
res <- merge(res, deaths[c("res_id", "site", "idno", "exit", "death")],
             by=c("site", "idno"), all.x=TRUE, suffixes=c("", ".death"))
res <- subset(res, is.na(exit.death) | exit <= exit.death)

table(res$death[res$res_id!= res$res_id.death])
res$death[res$res_id!= res$res_id.death] <- FALSE


## Combine overlapping residences

res <- res[order(res$site, res$idno, res$entry, res$death),]

res.events <- with(res, rbind(data.frame(site, idno, sex, dob, time=entry, type="entry", death=FALSE),
                              data.frame(site, idno, sex, dob, time=exit, type="exit", death=death)))
res.events <- res.events[do.call(order, res.events[c("site", "idno", "time")]),]
res.events$code <- c(1,-1)[match(res.events$type, c("entry", "exit"))]
res.events$cumexit <- cumsum(res.events$code)
res.events$cumentry <- rev(cumsum(rev(res.events$code)))

res.entry <- subset(res.events, cumentry == 0)
res.exit <- subset(res.events, cumexit == 0)
res.entry$entry <- res.entry$time
res.exit$exit <- res.exit$time

res.dedup <- cbind(res.entry[c("site", "idno", "sex", "dob", "entry")], res.exit[c("exit", "death")])

 ## Check de-duplicated residences
with(res.dedup, table(entry <= exit))

res <- res.dedup

###############
####  HIV  ####
###############

rakai.hiv <- lapply(mapply(c, "study_name", "idno",
                           paste("hiv_test_date", 1:10, sep=""),
                           paste("hiv_test_result", 1:10, sep=""),
                           paste("source_of_test_information", 1:10, sep=""), SIMPLIFY=FALSE),
                    function(x) setNames(rakai.hiv.wide[,x], c("study_name", "idno", "hiv_test_date", "hiv_test_result", "source_of_test_information")))
rakai.hiv <- do.call(rbind, lapply(rakai.hiv, subset, !is.na(hiv_test_date)))

masaka.hiv <- lapply(mapply(c, "study_name", "idno",
                           paste("hiv_test_date", 1:79, sep=""),
                           paste("hiv_test_result", 1:79, sep=""),
                           paste("source_of_test_information", 1:79, sep=""), SIMPLIFY=FALSE),
                     function(x) setNames(masaka.hiv.wide[,x], c("study_name", "idno", "hiv_test_date", "hiv_test_result", "source_of_test_information")))
masaka.hiv <- do.call(rbind, lapply(masaka.hiv, subset, !is.na(hiv_test_date)))


rakai.hiv$hiv_test_date <- as.Date(rakai.hiv$hiv_test_date, origin="1960-01-01")
kisesa.hiv$hiv_test_date <- as.Date(kisesa.hiv$hiv_test_date, origin="1960-01-01")

source.levels <- c("part of a population based study",
                   "part of a special research study",
                   "clinical record",
                   "self-reported by respondent",
                   "va report proxy respondent")
kisesa.hiv$source_of_test_information <- factor(kisesa.hiv$source_of_test_information, 1:5, source.levels)
masaka.hiv$source_of_test_information <- factor(as.integer(masaka.hiv$source_of_test_information), 1:5, source.levels)
rakai.hiv$source_of_test_information <- factor(as.integer(rakai.hiv$source_of_test_information), 1:5, source.levels)

hivresult.levels <- c("negative", "positive", "indeterminate", "not reported")
kisesa.hiv$hiv_test_result <- factor(kisesa.hiv$hiv_test_result, c(0, 1, 9, 3), hivresult.levels)
rakai.hiv$hiv_test_result <- factor(as.integer(rakai.hiv$hiv_test_result), 1:4, hivresult.levels)

hiv.cols <- c("study_name", "idno", "hiv_test_date", "hiv_test_result", "source_of_test_information")
hiv <- rbind(masaka.hiv[hiv.cols], rakai.hiv[hiv.cols], kisesa.hiv[hiv.cols])

## Recode variables
hiv$site <- factor(hiv$study_name)
hiv$testdate <- decimal_date_na(hiv$hiv_test_date)

hiv$result <- hiv$hiv_test_result
hiv$result[!hiv$result %in% c("negative", "positive")] <- NA
hiv$result <- hiv$result == "positive"

hiv$type <- factor(as.integer(hiv$source_of_test_information), 1:5, c("survey", "special study", "clinic", "self report", "VA"))

hiv <- subset(hiv, type == "survey" & !is.na(result), c(site, idno, testdate, result))


#####################
####  Fertility  ####
#####################

rakai.fert <- lapply(lapply(paste("m_delivery_date", 1:16, sep=""), c, c("idno", "study_name")),
                     function(x) setNames(rakai.fert.wide[,x], c("m_delivery_date", "idno", "study_name")))
rakai.fert <- do.call(rbind, lapply(rakai.fert, subset, !is.na(m_delivery_date)))

kisesa.fert <- kisesa.fert[,c("study_name", "idno", "m_delivery_date")]
masaka.fert <- masaka.fert[,c("study_name", "idno", "m_delivery_date")]

fert <- rbind(masaka.fert, rakai.fert, kisesa.fert)

## Round date to nearest day
fert$m_delivery_date <- as.Date(as.integer(fert$m_delivery_date), origin="1970-01-01")



## Recode variables
fert$site <- factor(fert$study_name)
fert$deliverydate <- decimal_date_na(fert$m_delivery_date)
fert <- fert[c("site", "idno", "deliverydate")]

## Keep births during residence episode, and only women
fert$fert_id <- seq_len(nrow(fert))

fert <- merge(fert, res, by=c("site", "idno"), all.x=TRUE)
table(fert$sex, useNA="always")
with(fert, table(deliverydate >= entry,  deliverydate <= exit))

fert <- subset(fert, sex == "Female" & deliverydate >= entry & deliverydate <= exit)
sum(duplicated(fert$fert_id))  # confirm: no births included twice (should be because residences are non-overlapping)

fert <- fert[c("site", "idno", "deliverydate")]

save(res, fert, hiv, file="alpha-fertility-data.RData")

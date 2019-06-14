###########################
# social integration report
# load data
# author: sebastian daza
############################

library(data.table)
library(lubridate)

# auxiliary functions
any_values = function(x, values) as.numeric(any(x %in% values, na.rm = TRUE))
amount_higher_zero = function(x, values) as.numeric(any(x > values))
reverse = function(x) (max(x, na.rm = TRUE) + 1) - x
table = function (...) base::table(..., useNA = 'ifany')

#################
# first week
#################

d1s = fread('data/180829_1_primerasemana.csv')
setnames(d1s, names(d1s), tolower(names(d1s)))
setnames(d1s, 'folio2', 'reg_folio')

# replace all missing values for NA
d1s = d1s[, lapply(.SD, function(x) ifelse(x < 0, NA, x)), .SDcols = names(d1s)]

# temporary housing
d1s[, temp_housing := apply(.SD, 1, any_values, c(1, 9, 11:13)),
    .SDcols = names(d1s) %like% '^p47$']

table(d1s$temp_housing)

d1s[, spent_night := apply(.SD, 1, any_values, c(1, 9, 11:13)),
    .SDcols = names(d1s) %like% 'p20_c_dia[0-9]$']

table(d1s$spent_night)

d1s[, temporary_housing := apply(.SD, 1, any_values, 1),
    .SDcols = c('temp_housing', 'spent_night')]

table(d1s$temporary_housing)

# money from familiy (including partner)
d1s[, money_family := apply(.SD, 1, any_values, 1),
    .SDcols = c('p43_3_a', 'p43_4_a')]

table(d1s$money_family)

# living with family
d1s[, living_with_family := apply(.SD, 1, any_values, 4:6),
    .SDcols = names(d1s) %like% '^p47$']

table(d1s$living_with_family)

d1s[, family_support := apply(.SD, 1, any_values, 1),
    .SDcols = c('living_with_family', 'money_family')]

table(d1s$family_support)

# formal job
d1s[, work_formal := p26]
d1s[is.na(work_formal) & !is.na(p22), work_formal := 0]
d1s[, work_informal := p35]
d1s[is.na(work_informal) & !is.na(p22), work_informal := 0]

d1s[, work := apply(.SD, 1, any_values, 1),
     .SDcols = c('p22', 'work_formal', 'work_informal')]

table(d1s$work)

# receiving money from public services
d1s[, money_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p43_1_a', 'p43_2_a')]

table(d1s$money_pp)

# contact with institution or programs

d1s[, contact_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p16_2', 'p16_3', 'p16_4', 'p16_6',
                 'p16_10')]

table(d1s$contact_pp)

d1s[, public_assistance := apply(.SD, 1, any_values, 1),
    .SDcols = c('money_pp', 'contact_pp')]

# select columns
d1s = d1s[, .(reg_folio,
              money_family, living_with_family, family_support,
              temp_housing, spent_night, temporary_housing,
              work_formal, work_informal, work,
              contact_pp, money_pp, public_assistance)]


###############
# 2 months
###############

d2m = fread('data/180118_2_dosmeses.csv')
setnames(d2m, names(d2m), tolower(names(d2m)))
setnames(d2m, 'folio2', 'reg_folio')

d2m = d2m[, lapply(.SD, function(x) ifelse(x < 0, NA, x)), .SDcols = names(d2m)]

# any job
d2m[, work_formal := apply(.SD, 1, any_values, 1),
    .SDcols = names(d2m) %like% 'trabajo_[0-9]_ocurrencia_[0-9]+']

d2m[is.na(work_formal) & p21==0, work_formal := 0]
d2m[is.na(p21) & !is.na(work_formal), work_formal := NA]
table(d2m$work_formal, d2m$p21)

d2m[, work_informal := apply(.SD, 1, any_values, 1),
    .SDcols = names(d2m) %like% 'tbjo_cp_[0-9]_[0-9]+']

d2m[is.na(work_informal) & p28 == 0, work_informal := 0]
d2m[is.na(p28) & !is.na(work_informal), work_informal := NA]
table(d2m$work_informal, d2m$p28)

d2m[, work := apply(.SD, 1, any_values, 1),
    .SDcols = c('work_formal', 'work_informal')]

table(d2m$work)

# temporary housing
d2m[, temp_housing := apply(.SD, 1, any_values, c(1, 9, 11:12)),
    .SDcols= names(d2m) %like% 'lugar_[0-9]_tipo$']

table(d2m$temp_housing)

# spent the night in a risky place
d2m[, spent_night := apply(.SD, 1, any_values, c(1:4, 8)),
    .SDcols = names(d2m) %like% 'paso_noche_[0-9]$']

d2m[, temporary_housing := apply(.SD, 1, any_values, 1),
    .SDcols = c('temp_housing', 'spent_night')]

# money from familiy (including partner)
d2m[, money_family := apply(.SD, 1, any_values, 1),
    .SDcols = c('p34_a_3', 'p34_a_4')]

table(d2m$money_family)

# living with family
d2m[, living_with_family := apply(.SD, 1, any_values, 4:6),
    .SDcols = names(d2m) %like% 'lugar_[0-9]_tipo$']

table(d2m$living_with_family)

d2m[, family_support := apply(.SD, 1, any_values, 1),
    .SDcols = c('living_with_family', 'money_family')]

table(d2m$family_support)

# receiving money from public services
d2m[, money_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p34_a_1', 'p34_a_2')]

table(d2m$money_pp)

# contact public services
d2m[, contact_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p7_a_2', 'p7_a_3', 'p7_a_4', 'p7_a_6',
                'p7_a_10')]

d2m[, public_assistance := apply(.SD, 1, any_values, 1),
    .SDcols = c('money_pp', 'contact_pp')]

# select columns
d2m = d2m[, .(reg_folio,
              money_family, living_with_family, family_support,
              temp_housing, spent_night, temporary_housing,
              work_formal, work_informal, work,
              contact_pp, money_pp, public_assistance)]


###############
# 6 months
###############

d6m = fread('data/180829_3_seismeses.csv')
setnames(d6m, names(d6m), tolower(names(d6m)))
setnames(d6m, 'folio2', 'reg_folio')

d6m = d6m[, lapply(.SD, function(x) ifelse(x < 0, NA, x)), .SDcols = names(d6m)]

# any job
# I am using all the available months
d6m[, work_formal := apply(.SD, 1, any_values, 1),
    .SDcols = names(d6m) %like% 'trabajo_rem_[0-9]_[0-9]+']

d6m[is.na(work_formal) & p17 == 0, work_formal := 0]
d6m[is.na(p17) & !is.na(work_formal), work_formal := NA]
table(d6m$work_formal, d6m$p17)

table(d6m$work_formal)

d6m[, work_informal := apply(.SD, 1, any_values, 1),
    .SDcols = names(d6m) %like% 'trabajo_cp_[0-9]_[0-9]+']

d6m[is.na(work_informal) & p23 == 0, work_informal := 0]
d6m[is.na(p23) & !is.na(work_informal), work_informal := NA]
table(d6m$work_informal, d6m$p23)

d6m[, work := apply(.SD, 1, any_values, 1),
    .SDcols = c('work_formal', 'work_informal')]

table(d6m$work)

# temporary housing
d6m[, temp_housing := apply(.SD, 1, any_values, c(1, 9, 11:12)),
    .SDcols= names(d6m) %like% 'lugar_[0-9]_tipo$']

table(d6m$temp_housing)

# spent the night in a risky place
d6m[, spent_night := apply(.SD, 1, any_values, 1),
    .SDcols= names(d6m) %like% 'dormiste_[0-9]+_[1-4|8]$']

d6m[, temporary_housing := apply(.SD, 1, any_values, 1),
    .SDcols = c('temp_housing', 'spent_night')]

table(d6m$temp_housing)

# money from familiy (including partner)
d6m[, money_family := apply(.SD, 1, any_values, 1),
    .SDcols = c('p28_a_3', 'p28_a_4')]

table(d6m$money_family)

# living with family
d6m[, living_with_family := apply(.SD, 1, any_values, 4:6),
    .SDcols = names(d6m) %like% 'lugar_[0-9]_tipo$']

table(d6m$living_with_family)

d6m[, family_support := apply(.SD, 1, any_values, 1),
    .SDcols = c('living_with_family', 'money_family')]

table(d6m$family_support)

# receiving money from public services
d6m[, money_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p28_a_1', 'p28_a_2')]

table(d6m$money_pp)

# contact with public services
d6m[, contact_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p7_a_2', 'p7_a_3', 'p7_a_4', 'p7_a_6',
                'p7_a_10')]

d6m[, public_assistance := apply(.SD, 1, any_values, 1),
    .SDcols = c('money_pp', 'contact_pp')]

# select columns
d6m = d6m[, .(reg_folio,
              money_family, living_with_family, family_support,
              temp_housing, spent_night, temporary_housing,
              work_formal, work_informal, work,
              contact_pp, money_pp, public_assistance)]


###############
# 12 months
###############

d12m = fread('data/180829_4_docemeses.csv')
setnames(d12m, names(d12m), tolower(names(d12m)))
setnames(d12m, 'folio2', 'reg_folio')

d12m = d12m[, lapply(.SD, function(x) ifelse(x < 0, NA, x)), .SDcols = names(d12m)]

# any job
# I am using all the available months
d12m[, work_formal := apply(.SD, 1, any_values, 1),
    .SDcols = names(d12m) %like% 'trabajo_[0-9]_mes_[0-9]+']

d12m[is.na(work_formal) & p22 == 0, work_formal := 0]
d12m[is.na(p22) & !is.na(work_formal), work_formal := NA]
table(d12m$work_formal, d12m$p22)

table(d12m$work_formal)

d12m[, work_informal := apply(.SD, 1, any_values, 1),
    .SDcols = names(d12m) %like% 'trabajo_cp_[0-9]_mes_[0-9]+']

d12m[is.na(work_informal) & p29 == 0, work_informal := 0]
d12m[is.na(p29) & !is.na(work_informal), work_informal := NA]
table(d12m$work_informal, d12m$p29)

d12m[, work := apply(.SD, 1, any_values, 1),
    .SDcols = c('work_formal', 'work_informal')]

table(d12m$work)

# temporary housing
d12m[, temp_housing := apply(.SD, 1, any_values,  c(1, 9, 11:12)),
    .SDcols= names(d12m) %like% 'lugar_[0-9]_cod$']

table(d12m$temp_housing)

# spent the night in a risky place
d12m[, spent_night := apply(.SD, 1, any_values, c(1:4, 8)),
    .SDcols= names(d12m) %like% 'dormiste_meses_[0-9]+']

table(d12m$spent_night)

d12m[, temporary_housing := apply(.SD, 1, any_values, 1),
    .SDcols = c('temp_housing', 'spent_night')]

table(d12m$temporary_housing)

# money from familiy (including partner)
d12m[, money_family := apply(.SD, 1, any_values, 1),
    .SDcols = c('p34_a_3', 'p34_a_4')]

table(d12m$money_family)

# living with family
d12m[, living_with_family := apply(.SD, 1, any_values, 4:6),
    .SDcols = names(d12m) %like% 'lugar_[0-9]_cod$']

table(d12m$living_with_family)

d12m[, family_support := apply(.SD, 1, any_values, 1),
    .SDcols = c('living_with_family', 'money_family')]

table(d12m$family_support)

# receiving money from public services
d12m[, money_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p34_a_1', 'p34_a_2')]

table(d12m$money_pp)

# contact with public services
d12m[, contact_pp := apply(.SD, 1, any_values, 1),
    .SDcols = c('p8_a_3', 'p8_a_4', 'p8_a_5', 'p8_a_7',
                'p8_a_11')]

d12m[, public_assistance := apply(.SD, 1, any_values, 1),
    .SDcols = c('money_pp', 'contact_pp')]

# select columns
d12m = d12m[, .(reg_folio,
                money_family, living_with_family, family_support,
                temp_housing, spent_night, temporary_housing,
                work_formal, work_informal, work,
                contact_pp, money_pp, public_assistance)]


#####################
# outcome
#####################

outcome = rbindlist(list(d1s, d2m, d6m, d12m), idcol = 'time')
setorder(outcome, reg_folio, time)

# descriptives

outcome[, .(
            money_family = mean(money_family, na.rm = TRUE),
            living_with_family = mean(living_with_family, na.rm = TRUE),
            family_support = mean(family_support, na.rm = TRUE),
            temp_housing = mean(temp_housing, na.rm = TRUE),
            spent_night = mean(spent_night, na.rm = TRUE),
            temporary_housing = mean(temporary_housing, na.rm = TRUE),
            work = mean(work, na.rm = TRUE),
            work_formal = mean(work_formal, na.rm = TRUE),
            work_informal = mean(work_informal, na.rm = TRUE),
            money_pp = mean(money_pp, na.rm = TRUE),
            money_pp = mean(money_pp, na.rm = TRUE),
            contact_pp = mean(contact_pp, na.rm = TRUE),
            public_assistance = mean(public_assistance, na.rm = TRUE)
            )
        , time]
outcome[, .N, time]


#####################
# baseline variables
#####################

bs = fread('data/base_general_long.csv')
setnames(bs, names(bs), tolower(names(bs)))
bs = bs[reg_ola == 0 & reg_muestra == 1]
bs = bs[, lapply(.SD, function(x) ifelse(x < 0, NA, x)), .SDcols = names(bs)]

# add  latent classes
lclass  = fread('data/Folio_clase.csv')
nvars = c('reg_folio', 'prob1', 'prob2', 'prob3', 'probmax', 'class')
setnames(lclass, names(lclass), nvars)
lclass = lclass[, .(reg_folio, class)]

setkey(lclass, reg_folio)
setkey(bs, reg_folio)
bs = lclass[bs]
table(bs$class)

# age
setnames(bs, 'hdv_1', 'age')

# partner before prison
bs[, previous_partner := ifelse(par_6 == 1 | par_4==1 | par_4==2, 1, 0)]
bs[is.na(previous_partner) & par_4 == 0, previous_partner := 0]
table(bs$previous_partner)

# primary school dropout
bs[, only_primary := ifelse(hdv_7 <= 8, 1, 0)]
table(bs$only_primary)

bs[hdv_7 == 0, edu := 'none']
bs[hdv_7 %in% 1:7, edu := 'primary incomplete']
bs[hdv_7 == 8, edu := 'primary complete']
bs[hdv_7 %in% 9:11, edu := 'secondary incomplete']
bs[hdv_7 == 12, edu := 'secondary complete']
bs[hdv_7 > 12, edu := 'some terciary']

table(bs$edu)

# self efficacy
self_efficacy_vars = paste0('car_1_', 10:16)

bs[, self_efficacy_vars[c(2, 5)] := lapply(.SD, reverse),
   .SDcols = self_efficacy_vars[c(2, 5)]]

bs[, self_efficacy := scale(apply(.SD, 1, mean, na.rm=TRUE)),
    .SDcols = self_efficacy_vars]

# desire for help
help_vars = paste0('spg_8_', 1:4)

bs[, c(help_vars) := lapply(.SD, reverse),
   .SDcols = help_vars]

bs[, desire_change := scale(apply(.SD, 1, mean, na.rm=TRUE)),
    .SDcols = help_vars]

# work before prison
bs[, any_previous_work := apply(.SD, 1, any_values, 1),
    .SDcols = c('eaf_6', 'eaf_43')]

table(bs$any_previous_work)

# type of crime
table(bs$del_10)
bs[del_10 %in% c(1:6, 8:11, 17:18, 21:22), crime := 'property']
bs[del_10 == 7, crime := 'theft']
bs[del_10 %in% c(12:14, 19:20), crime := 'person']
bs[del_10 %in% c(15:16), crime := 'drug']
bs[del_10 == 23, crime := 'other']
bs[del_10 == 24 & reg_folio %in% c(30039, 30303), crime := 'other']
bs[reg_folio == 50129, crime := 'property']

table(bs$crime)

# previous sentences
setnames(bs, 'del_5_1', 'previous_sentences')
bs[is.na(previous_sentences) & del_4 == 0, previous_sentences := 0]

# time in prison (months)
table(bs$del_6_1) # months
table(bs$del_6_2) # years
table(bs$del_6_3) # days

bs[, del_6_2 := del_6_2 * 12]
bs[, del_6_3 := del_6_3 / 30.5]

bs[, total_previous_months_in_prison := apply(.SD, 1, sum, na.rm=TRUE),
    .SDcols = paste0('del_6_', 1:3)]

table(bs$total_previous_months_in_prison)

# add current time
bs[, current_time_in_prison := interval(ymd(reg_fprivacion), ymd(reg_fegreso)) %/% months(1)]
bs[, total_months_in_prison := apply(.SD, 1, sum, na.rm=TRUE),
    .SDcols = c('total_previous_months_in_prison', 'current_time_in_prison')]

table(bs$reg_fegreso)
table(bs$reg_fprivacion)

# report date of prison
# remove missing using mid day and month
bs[is.na(hdv_4_1), hdv_4_1 := 15]
bs[is.na(hdv_4_2), hdv_4_2 := 7]

bs[, report_date_prison := paste0(hdv_4_1, '-', hdv_4_2, '-', hdv_4_3)]

bs[total_months_in_prison == 0, .(reg_folio, total_months_in_prison, del_6_1, del_6_2, del_6_3,
                                  reg_fecha, reg_fprivacion, reg_fegreso, hdv_4_1, hdv_4_2, hdv_4_3,
                                  report_date_prison)]

bs[total_months_in_prison == 0, current_time_in_prison := interval(dmy(report_date_prison), ymd(reg_fecha)) %/% months(1)]
bs[total_months_in_prison == 0, total_months_in_prison := apply(.SD, 1, sum, na.rm=TRUE),
    .SDcols = c('total_previous_months_in_prison', 'current_time_in_prison')]

bs[total_months_in_prison == 0, .(reg_folio, total_months_in_prison, del_6_1, del_6_2, del_6_3,
                                  reg_fecha, reg_fprivacion, reg_fegreso, hdv_4_1, hdv_4_2, hdv_4_3,
                                  report_date_prison)]

hist(bs$total_months_in_prison)

# mental health
bs[, mental_health := scale(apply(.SD, 1, mean, na.rm = TRUE)),
    .SDcols = names(bs) %like% '^sal_31']

# hijos
setnames(bs, 'hij_1', 'nchildren')
bs[, any_children := ifelse(nchildren > 0, 1, 0)]

# early crime
bs[, early_crime := ifelse(del_15 < 15, 1, 0)]
table(bs$early_crime)

# last sentence extension
bs[, del_11_2 := del_11_2 / 12]
bs[, del_11_3 := del_11_3 / 365.25]
bs[, sentence_length := apply(.SD, 1, sum, na.rm = TRUE),
    .SDcols = paste0('del_11_', 1:3)]

bs[sentence_length == 0, sentence_length := NA]

# drug abuse and dependence
# most frequent drug
dep = paste0('dro_6_', c(1, 2, 4, 5, 6, 7))
abuse = paste0('dro_6_', c(8, 9, 10, 11))

bs[, c(dep, abuse) := lapply(.SD, function (x) ifelse(is.na(x), 0, x)),
    .SDcols = c(dep, abuse)]

bs[, dep_1 := apply(.SD, 1, sum, na.rm = TRUE),
    .SDcols = dep]
bs[, abuse_1 := apply(.SD, 1, sum, na.rm = TRUE),
    .SDcols = abuse]

# more problematic drug
dep = paste0('dro_7_', c(1, 2, 4, 5, 6, 7))
abuse = paste0('dro_7_', c(8, 9, 10, 11))

bs[, c(dep, abuse) := lapply(.SD, function (x) ifelse(is.na(x), 0, x)),
    .SDcols = c(dep, abuse)]

bs[, dep_2 := apply(.SD, 1, sum, na.rm = TRUE),
    .SDcols = dep]
bs[, abuse_2 := apply(.SD, 1, sum, na.rm = TRUE),
    .SDcols = abuse]

bs[, drug_dep := apply(.SD, 1, any_values, 3:6),
    .SDcols = names(bs) %like% '^dep_[1-2]$']

bs[, drug_abuse := apply(.SD, 1, any_values, 1:4),
    .SDcols = names(bs) %like% '^abuse_[1-2]$']

bs[, drug_depabuse := apply(.SD, 1, any_values, 1),
    .SDcols = c('drug_dep', 'drug_abuse')]

# family support and conflict
fsupport = paste0('sfa_8_', 1:4)
fconflict = paste0('sfa_8_', 5:7)

bs[, c(fsupport) := lapply(.SD, reverse),
   .SDcols = c(fsupport)]

bs[, c(fconflict) := lapply(.SD, reverse),
   .SDcols = c(fconflict)]

bs[, family_conflict := scale(apply(.SD, 1, mean, na.rm = TRUE)),
    .SDcols = c(fconflict)]

# select columns
bs = bs[, .(reg_folio, class, age, edu, only_primary, any_previous_work,
    nchildren, any_children, crime, early_crime,
    self_efficacy, desire_change, previous_partner,
    previous_sentences, mental_health, drug_depabuse,
    sentence_length, total_months_in_prison, family_conflict
    )]

# center variables
cvars = c('age', 'sentence_length', 'nchildren', 'previous_sentences')
bs[, paste0('c_', cvars) := lapply(.SD, scale, scale=FALSE), .SDcols = cvars]

# index to expand database
indx = bs[ , list(reg_folio = reg_folio, time = 1:4), by = 1:nrow(bs)
          ][, nrow := NULL]

# merge with outcome
setkey(bs, reg_folio)
setkey(indx, reg_folio)

# expand bs
bs = bs[indx]

# merge with outcome
setkey(outcome, reg_folio, time)
setkey(bs, reg_folio, time)

df = outcome[bs]
summary(df)

# explore some cases
ids = unique(df$reg_folio)
df[reg_folio == sample(ids, 1)]

saveRDS(df, 'output/clean_data.rd')

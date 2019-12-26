###########################∏∏
# social integration report
# analysis
# author: sebastian daza
############################

# devtools::install_github("leifeld/texreg")
library(data.table)
library(brms)
library(sdazar)
library(mice)
library(miceadds)
library(micemd)
library(texreg)
library(future)
library(stringr)
library(ggplot2)

source('src/utils.R')

# load data
df = readRDS('output/clean_data.rd')
names(df)

fvars = c('edu', 'crime', 'time', 'class')
df[, c(fvars) := lapply(.SD, factor), .SDcols = fvars]

# binary variables
fvars = c('previous_partner', 'any_previous_work', 'drug_depabuse', 'only_primary')
df[, c(fvars) := lapply(.SD, factor, labels = c('No', 'Sí')), .SDcols = fvars]

# class
df[, class := factor(class, labels = c('Clase 1', 'Clase 2', 'Clase 3'))]

# imputation
mdf = df[, .(reg_folio, time, money_family, living_with_family,
              temp_housing, spent_night, work_informal, work_formal,
              contact_pp, money_pp, age, only_primary, class,
              nchildren, previous_partner,
              crime, any_previous_work,
              self_efficacy, desire_change,
              previous_sentences, sentence_length, total_months_in_prison,
              mental_health, drug_depabuse, family_conflict)]

# explore correlations
cor(mdf[, .SD, .SDcols = sapply(mdf, is.numeric)], use = 'complete.obs')

cor(mdf[, .(previous_sentences, total_months_in_prison, sentence_length)], use = 'complete.obs')

cor(mdf[, .(self_efficacy, desire_change,
  family_conflict, mental_health)], use = 'complete.obs')

predM = mice::make.predictorMatrix(data = mdf)
predM[, "reg_folio"] = -2
impMethod = mice::make.method(data = mdf)
l2vars = c('previous_sentences', 'family_conflict', 'self_efficacy')
predM[l2vars, 'time'] = 0
impMethod[l2vars] <- '2lonly.function'
imputationFunction <- list('previous_sentences' = 'pmm' ,
                           'family_conflict' = 'pmm',
                           'self_efficacy' = 'pmm')
cluster_var <- list('previous_sentences' = 'reg_folio',
                    'family_conflict' = 'reg_folio',
                    'self_efficacy' = 'reg_folio')

#impMethod
number_imputations = 20

print('::::::: running imputations')
imp = mice::mice(mdf, predictorMatrix = predM, m = number_imputations, maxit = 30,
                 method = impMethod,
                 imputationFunction = imputationFunction,
                 cluster_var = cluster_var,
                 print = FALSE)

# explore quality of the imputation
plot(imp)
densityplot(imp)

# model to create wave plots (0)

file_model_0 = 'output/brms_model_0.rd'
if(!file.exists(file_model_0)) {

    print('::::::: running baseline model')
    plan(multiprocess)
    m0 = brm_multiple(mvbind(money_family, living_with_family, temp_housing, spent_night,
                      work_formal, work_informal, money_pp, contact_pp) ~
                      time,
                      data = imp,
                      family = bernoulli(),
                      control = list(adapt_delta=0.90),
                      chains = 2)

    saveRDS(m0, file = file_model_0)

} else { m0 = readRDS(file_model_0) }

check_convergence_mi(m0)

# create vector with dependent variables
depvars = c('moneyfamily', 'livingwithfamily', 'temphousing',
  'spentnight', 'workformal', 'workinformal', 'moneypp', 'contactpp')

# loop to create plots
fitted_values = list()
newdata = data.frame(time = factor(1:4))

for (i in depvars) {
  print(paste0('fitting values for ', i))
  pp = fitted(m0, newdata = newdata, resp = i,
              summary = TRUE, robust = TRUE)
  fitted_values[[i]] = data.table(pp)[, dep := i][, time := 1:4]
}

fitted_values = rbindlist(fitted_values)

fitted_values[, time := factor(time,
              levels = c(1:4),
              labels =c('Primera semana', 'Dos meses', 'Seis meses', 'Doce meses'))]

# explore fitted values summary
fitted_values

# create labels for dependent variables
cnames = c('Dinero familiares', 'Vive con familiares', 'Vivienda temporal', 'Noche en lugar de riesgo',
           'Trabajo formal', 'Trabajo informal', 'Subsidios', 'Contacto instituciones')

fitted_values[, dep := factor(dep, labels = cnames, levels = depvars)]

depvars_list = list(
  'family' = cnames[1:2],
  'housing'  = cnames[3:4],
  'work'  = cnames[5:6],
  'public'  = cnames[7:8]
  )

# function to create descriptive plot by dependent variable
create_plots_outcome = function(depvars, plotname) {

  savepdf(paste0('output/', plotname))
  print(
  ggplot(fitted_values[dep %in% depvars],
        aes(x=time, y=Estimate, group=dep, color=dep)) +
  geom_pointrange(aes(ymin=Q2.5, ymax=Q97.5), position = position_dodge(width=0.2)) +
  labs(x='\nOla\n', y='Probabilidad\n',
       caption = paste0("Nota: Intervalos de credibilidad (95%), ",
                        number_imputations, " imputaciones")) +
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0,.80)) +
  scale_color_manual(values=c("#f03b20", "#2c7fb8")) +
  theme_classic() +
  theme(legend.position="top", legend.title=element_blank(),
    plot.caption = element_text(hjust = 0, size = 6, face = "italic"))
  )
  dev.off()

}

for (i in seq_along(depvars_list)) {
  create_plots_outcome(depvars_list[[i]], names(depvars_list[i]))
}

# remove object to clear memory
remove(m0)

# initial descriptive model (1)

file_model_1 = 'output/brms_model_1.rd'

if(!file.exists(file_model_1)) {
    print('::::::: running model 1')
    plan(multiprocess)
    m1 = brm_multiple(mvbind(money_family, living_with_family, temp_housing,
                      spent_night, work_formal, work_informal, money_pp, contact_pp) ~
                      time + age + nchildren + only_primary + previous_partner +
                      self_efficacy + desire_change + any_previous_work +
                      #  family_conflict +
                      mental_health + drug_depabuse + total_months_in_prison +
                      previous_sentences + (1|p|reg_folio),
                      data = imp,
                      family = bernoulli(),
                      control = list(adapt_delta=0.95),
                      chains = 2)
    saveRDS(m1, file = file_model_1)


} else { m1 = readRDS(file_model_1) }

check_convergence_mi(m1, high = 1.1)

# screenreg(m1, omit.coef = paste0(depvars[-1], collapse='|'),
#   include.r2=TRUE)

# create latex table model 1
cmap = list('Intercept' = 'Constante',
            'time2' = 'Dos meses',
            'time3' = 'Seis meses',
            'time4' = 'Doce meses',
            'age' = 'Edad',
            'only_primarySí' = 'Educación básica o menos',
            'nchildren' = 'Número de hijos',
            'previous_partnerSí' = 'Pareja antes de la cárcel',
            'any_previous_workSí' = 'Trabajo antes de la cárcel',
            'mental_health' = 'Problemas de salud mental',
            'drug_depabuseSí' = 'Dependencia / abuso drogas',
            'self_efficacy' = 'Autoeficacia',
            'desire_change' = 'Disposición al cambio',
            'family_conflict' = 'Escala conflicto familiar',
            'classClase2' = 'Clase 2',
            'classClase3' = 'Clase 3',
            'crime' = 'Delito condena',
            'sentence_length' = 'Tiempo condena',
            'previous_sentences' = 'Número de condenas previas',
            'total_months_in_prison' = 'Tiempo total en la cárcel (meses)')

dep_regular_exp = c('^moneyfamily_', '^livingwithfamily_', '^temphousing_',
  '^spentnight_', '^workformal_', '^workinformal_',
  '^moneypp_', '^contactpp_')

list_texreg = create_texreg_multivariate(m1, dep_regular_exp,
              include.r2 = TRUE, rhat_max = 1.1)

ndeps = length(dep_regular_exp)

tab = texreg(list_texreg,
          custom.coef.map = cmap,
          custom.model.names = cnames,
          groups = list('Ola (ref = primera semana)' = 2:4),
          ci.test = NULL,
          float.pos = "htp",
          caption = paste0('Modelo Bayesiano multivariable (',
                           ndeps, ' variables dependientes)'),
          booktabs = TRUE,
          use.packages = FALSE,
          dcolumn = TRUE,
          caption.above = TRUE,
          scalebox = 0.70,
          # fontsize = 'scriptsize',
          label = "integracion_social_m1",
          sideways = TRUE,
          digits = 2,
          custom.note = paste0("Intervalos de credibilidad 95\\%. Coeficientes corresponden a un modelo con ", ndeps, " variables dependientes.
          Efectos aleatorios y correlaciones entre variables dependientes son omitidos.")
          )

# clean up table
tab = str_replace(tab, 'Num\\. obs\\.  reg\\\\_folio', 'Número mujeres')
tab = str_replace(tab, 'Num\\. obs\\.', 'Número observaciones')
top = "\\\\toprule
\\\\addlinespace
& \\\\multicolumn{2}{c}{Familia} &  \\\\multicolumn{2}{c}{Precariedad Residencial} &
\\\\multicolumn{2}{c}{Trabajo} &
\\\\multicolumn{2}{c}{Ayuda Institucional} \\\\\\\\
\\\\addlinespace
\\\\addlinespace
& \\\\multicolumn{1}{c}{Dinero familiares} & \\\\multicolumn{1}{c}{Vive con familiares} & \\\\multicolumn{1}{c}{Vivienda temporal} &
\\\\multicolumn{1}{c}{Noche en lugar de riesgo} & \\\\multicolumn{1}{c}{Trabajo formal} &
\\\\multicolumn{1}{c}{Trabajo informal} & \\\\multicolumn{1}{c}{Subsidios} & \\\\multicolumn{1}{c}{Contacto instituciones} \\\\\\\\
\\\\addlinespace
\\\\addlinespace
\\\\midrule"

tab = str_replace(tab, '\\\\toprule\\n.+\\n\\\\midrule', top)
cat(tab,  file = 'output/integracion_social_m1.tex')

# function to explore marginal effects
marginal_plots = function(model, covariate, depvar_location,
                          names_depvars, map_covariates) {
    output = marginal_effects(model, covariate)
    g = plot(output , plot = FALSE)[[depvar_location]] +
        theme_classic() +
        labs(x = paste0('\n', map_covariates[[covariate]]),
             y = paste0(names_depvars[depvar_location],'\n')) +
        scale_color_grey() +
        scale_fill_grey() +
        theme(axis.text.x = element_text(size = 22),
              axis.text.y = element_text(size = 22),
              axis.title.x = element_text(size = 22),
              axis.title.y = element_text(size = 22))
    print(g)
  }

# save marginal plots

# new cmap for axis labels
cmap = list('Intercept' = 'Constante',
            'time2' = 'Dos meses',
            'time3' = 'Seis meses',
            'time4' = 'Doce meses',
            'age' = 'Edad',
            'only_primary' = 'Educación básica o menos',
            'nchildren' = 'Número de hijos',
            'previous_partner' = 'Pareja antes de la cárcel',
            'any_previous_work' = 'Trabajo antes de la cárcel',
            'mental_health' = 'Problemas de salud mental',
            'drug_depabuse' = 'Dependencia / abuso drogas',
            'self_efficacy' = 'Autoeficacia',
            'desire_change' = 'Disposición al cambio',
            'family_conflict' = 'Escala conflicto familiar',
            'classClase2' = 'Clase 2',
            'classClase3' = 'Clase 3',
            'crime' = 'Delito condena',
            'sentence_length' = 'Tiempo condena',
            'previous_sentences' = 'Número de condenas previas',
            'total_months_in_prison' = 'Tiempo total en la cárcel (meses)')

output_index = c('moneyfamily',
                 'livingwithfamily', 'livingwithfamily', 'livingwithfamily',
                 'temphousing', 'temphousing',
                 'spentnight', 'spentnight',
                 'workformal', 'workformal', 'workformal',
                 'workinformal', 'workinformal',
                 'moneypp', 'moneypp',
                 'contactpp')

variables_prediction = c('previous_sentences',
                        'age', 'nchildren', 'mental_health',
                        'drug_depabuse', 'previous_sentences',
                        'mental_health', 'drug_depabuse',
                        'age', 'self_efficacy', 'previous_sentences',
                        'age', 'any_previous_work',
                        'nchildren', 'total_months_in_prison',
                        'drug_depabuse')

if (! length(output_index) == length(variables_prediction)) {
  stop ('Number of outcome variables do not match predictors in marginal plots')
}

for (i in seq_along(output_index)) {
    print(paste0(':::::::: plotting ', output_index[i], ' - ', variables_prediction[i]))
    marginal_plots(m1, variables_prediction[i],  which(depvars == output_index[i]), cnames, cmap)
    ggsave(paste0('output/', output_index[i], '_', gsub('_', '', variables_prediction[i]), '.pdf'),
           width = 9, height = 6)
}

# remove model to clear memory
remove(m1)

# model with classes and some demographic controls (2)

file_model_2 = 'output/brms_model_2.rd'

if(!file.exists(file_model_2)) {

    print('::::::: running model 2')
    plan(multiprocess)
    m2 = brm_multiple(mvbind(money_family, living_with_family, temp_housing,
                      spent_night, work_formal, work_informal, money_pp, contact_pp) ~
                      time + age + nchildren + only_primary + previous_partner +
                      class + (1|p|reg_folio),
                      data = imp,
                      family = bernoulli(),
                      control = list(adapt_delta = 0.97),
                      chains = 2)

    saveRDS(m2, file = file_model_2)

} else { m2 = readRDS(file_model_2) }

check_convergence_mi(m2, high = 1.1)

list_texreg = create_texreg_multivariate(m2, dep_regular_exp,
              include.r2 = TRUE, rhat_max = 1.1)

ndeps = length(dep_regular_exp)

# create latex table model 2
tab = texreg(list_texreg,
          custom.coef.map = cmap,
          custom.model.names = cnames,
          groups = list('Ola (ref = primera semana)' = 2:4,
                        'Perfil (ref = Clase 1)' = 5:6),
          ci.test = NULL,
          float.pos = "htp",
          caption = paste0('Modelo Bayesiano multivariable (',
                           ndeps, ' variables dependientes)'),
          booktabs = TRUE,
          use.packages = FALSE,
          dcolumn = TRUE,
          caption.above = TRUE,
          scalebox = 0.70,
          # fontsize = 'scriptsize',
          label = "integracion_social_m2",
          sideways = TRUE,
          digits = 2,
          custom.note = paste0("Intervalos de credibilidad 95\\%. Coeficientes corresponden a un modelo con ", ndeps, " variables dependientes.
          Efectos aleatorios y correlaciones entre variables dependientes son omitidos. Modelo además ajusta por edad, número de hijos, y educación
          básica o menos.")
          )

# clean up table
tab = str_replace(tab, 'Num\\. obs\\.  reg\\\\_folio', 'Número mujeres')
tab = str_replace(tab, 'Num\\. obs\\.', 'Número observaciones')
top = "\\\\toprule
\\\\addlinespace
& \\\\multicolumn{2}{c}{Familia} &  \\\\multicolumn{2}{c}{Precariedad Residencial} &
\\\\multicolumn{2}{c}{Trabajo} &
\\\\multicolumn{2}{c}{Ayuda Institucional} \\\\\\\\
\\\\addlinespace
\\\\addlinespace
& \\\\multicolumn{1}{c}{Dinero familiares} & \\\\multicolumn{1}{c}{Vive con familiares} & \\\\multicolumn{1}{c}{Vivienda temporal} &
\\\\multicolumn{1}{c}{Noche en lugar de riesgo} & \\\\multicolumn{1}{c}{Trabajo formal} &
\\\\multicolumn{1}{c}{Trabajo informal} & \\\\multicolumn{1}{c}{Subsidios} & \\\\multicolumn{1}{c}{Contacto instituciones} \\\\\\\\
\\\\addlinespace
\\\\addlinespace
\\\\midrule"

tab = str_replace(tab, '\\\\toprule\\n.+\\n\\\\midrule', top)
cat(tab,  file = 'output/integracion_social_m2.tex')

# create marginal plots for classes
cmap = list('Intercept' = 'Constante',
            'time2' = 'Dos meses',
            'time3' = 'Seis meses',
            'time4' = 'Doce meses',
            'classClase2' = 'Clase 2',
            'classClase3' = 'Clase 3')

output_index = c('temphousing',
                 'spentnight',
                 'workinformal',
                 'contactpp')

variables_prediction = c('class', 'class', 'class', 'class')

if (! length(output_index) == length(variables_prediction)) {
  stop ('Number of outcome variables do not match predictors in marginal plots')
}

for (i in seq_along(output_index)) {
    print(paste0(':::::::: plotting ', output_index[i], ' - ', variables_prediction[i]))
    marginal_plots(m2, variables_prediction[i],  which(depvars == output_index[i]), cnames, cmap)
    ggsave(paste0('output/', output_index[i], '_', gsub('_', '', variables_prediction[i]), '.pdf'),
           width = 9, height = 6)
}

remove(m2)

# end analysis

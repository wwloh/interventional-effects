rm(list=ls())
libraries_check <- c("data.table","foreign")
for (libs in libraries_check) {
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# load dataset
path <- file.path("Self-efficacy and fatigue among health care workers during COVID-19 outbreak A moderated mediation model of posttraumatic stress disorder symptoms and negative coping.sav")
Data <- data.table(read.spss(path, to.data.frame=TRUE))
setkey(Data)
dim(Data)
rm(path)

# specify variables ===========================================================
Data <- Data[,c(paste0("Q",c(1:4,69,70)), # covariates
                "Copingstyle_meanneg", # negative coping
                "GSES原始总分", # treatment 
                "PCLC_reexperiencing", # mediators (subscales of PTSD symptoms)
                "PCLC_avoidance",
                "PCLC_hyperarousal", 
                "Fatigue_new" # outcome
                ), with=FALSE]
setkey(Data)
# check for missing data
any(rowSums(Data[, lapply(.SD, is.na)])>0)
dim(Data)
setnames(Data,c("age","gender","years_working","marital","edu","tech",
                "neg_cop","A",paste0("M",1:3),"Y"))
summary(Data)
## convert factors to integers
Data[, age.f := as.character(levels(age))[age]]
Data[, age := as.integer(age.f)]
Data[sample(nrow(Data),10),list(age.f,age)] # random sample to check
Data[, age.f := NULL]
Data[, gender.f := as.character(levels(gender))[gender]]
Data[, gender := as.integer(gender.f=="女")]
table(Data[,list(gender.f,gender)]) # check
Data[, gender.f := NULL]
Data[, years_working.f := as.character(levels(years_working))[years_working]]
Data[, years_working := as.integer(sapply(Data[, years_working.f], function(x) 
  strsplit(x,split="-")[[1]][1])>1)]
table(Data[,list(years_working.f,years_working)]) # check
Data[, years_working.f := NULL]
Data[, marital.f := as.character(levels(marital))[marital]]
Data[, marital := as.integer(marital.f=="未婚")]
table(Data[,list(marital.f,marital)]) # check
Data[, marital.f := NULL]
Data[, edu.f := as.character(levels(edu))[edu]]
Data[, edu := as.integer(edu.f=="大学本科")]
table(Data[,list(edu.f,edu)]) # check
Data[, edu.f := NULL]
Data[, tech.f := as.character(levels(tech))[tech]]
Data[, tech := sapply(Data[, tech.f], function(x) {
  if (x=="中级") return (1L)
  if (x=="副高") return (2L)
  return (0L)
})]
table(Data[,list(tech.f,tech)]) # check
Data[, tech.f := NULL]
# standardized negative coping
Data[, neg_cop := scale(Data$neg_cop,center=TRUE,scale=TRUE)]
# dichotomize treatment
Data[, A := as.integer(A >= median(A))]
# dichotomize outcome: 1 (0) is presence (absence) of fatigue
Data[, Y := as.integer(Y >= 7)]
Data <- cbind("id"=1:nrow(Data),Data)
setkey(Data)
summary(Data)

m_labels <- c("reexperiencing","avoidance","hyperarousal")

# check (conditional) association between each mediator and treatment
for (s in 1:3) {
  print(summary(lm(paste0("M",s,"~A+",paste0(names(Data)[2:7],collapse="+")),
                   data=Data))  )
}

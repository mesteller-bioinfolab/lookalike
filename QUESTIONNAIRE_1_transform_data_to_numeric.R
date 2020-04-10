setwd("data/")


sample_info_csv <- read.table("INPUT_curated_questionaire.csv", header = T, sep = ",", stringsAsFactors = F)

rownames(sample_info_csv) <- sample_info_csv$ID_PEBC


# I remove physical traits and of place of origin 
physical_or_origin_traits <- c("Country",
                               "Gender",
                               "Language",
                               "Race",
                               "Current.Recidence.Category",
                               "Weight",
                               "Height",
                               "Hair.Color",
                               "Hair.Shape",
                               "Eyes.Color",
                               "Age_turning_in_2016",
                               "Country_born",
                               "City_village_born",
                               "Country_residence",
                               "City_village_residence",
                               "Lived_elsewhere")

sample_info_csv <- sample_info_csv[, !colnames(sample_info_csv) %in% physical_or_origin_traits]


### ORDERED FACTORS: (ordering and assigning numbers to ordered factors)

### Glasses
sample_info_csv$Use.Glasses. <- as.numeric(factor(sample_info_csv$Use.Glasses., levels = c("No", "Sometimes", "Yes"), ordered = T))

### Education
sample_info_csv$Education <- as.factor(sample_info_csv$Education )
levels(sample_info_csv$Education) #  "Basic School"           "Elementary School"      "High School"            
#  "Master"                 "Other"                  "Other (Upper level training)" "PhD"                    "University"     
levels(sample_info_csv$Education) <- c(1,1,2,4,2.5,2.5,5,3)
sample_info_csv$Education <- as.numeric(as.character(sample_info_csv$Education))


### Smoking
sample_info_csv$Smoking <- as.factor(sample_info_csv$Smoking )
levels(sample_info_csv$Smoking) # "Former" "Heavy"  "Light"  "Never" 
levels(sample_info_csv$Smoking) <- c(3,4,2,1)
sample_info_csv$Smoking <- as.numeric(as.character(sample_info_csv$Smoking))

### "Diet" # From more plant-based to less plant-based
sample_info_csv$Diet <- as.factor(sample_info_csv$Diet )
levels(sample_info_csv$Diet) # "Not Vegetarian"   "Not Vegetarian (former vegetarian)" "Vegan"  "Vegetarian"  "Vegetarian (mostly)"    
levels(sample_info_csv$Diet) <- c(5, 4, 1, 2,3)
sample_info_csv$Diet <- as.numeric(as.character(sample_info_csv$Diet))

### Alcohol
sample_info_csv$Drinking.Alcohol <- as.factor(sample_info_csv$Drinking.Alcohol )
levels(sample_info_csv$Drinking.Alcohol) #  "Never"           "Often"           "Sometimes"       "Sometimes-Often"
levels(sample_info_csv$Drinking.Alcohol) <- c(1,3,2,2.5)
sample_info_csv$Drinking.Alcohol <- as.numeric(as.character(sample_info_csv$Drinking.Alcohol))

### Soft-drinks
sample_info_csv$Drinking.Soft.drinks <- as.factor(sample_info_csv$Drinking.Soft.drinks )
levels(sample_info_csv$Drinking.Soft.drinks) #  "Never"           "Often"           "Sometimes"      
levels(sample_info_csv$Drinking.Soft.drinks) <- c(1,3,2)
sample_info_csv$Drinking.Soft.drinks <- as.numeric(as.character(sample_info_csv$Drinking.Soft.drinks))

# Work indoor-outdoor
sample_info_csv$Employment_in_outdoor <- as.factor(sample_info_csv$Employment_in_outdoor )
levels(sample_info_csv$Employment_in_outdoor) # "Indoor"         "Indoor/Outdoor" "Outdoor"
levels(sample_info_csv$Employment_in_outdoor) <- c(0,1,2)
sample_info_csv$Employment_in_outdoor <- as.numeric(as.character(sample_info_csv$Employment_in_outdoor))

## RH
sample_info_csv$RH <- as.factor(sample_info_csv$RH )
levels(sample_info_csv$RH) # "-", "+"
levels(sample_info_csv$RH) <- c(0,1)
sample_info_csv$RH <- as.numeric(as.character(sample_info_csv$RH))

## Handwriting
sample_info_csv$Handwriting <- as.factor(sample_info_csv$Handwriting )
levels(sample_info_csv$Handwriting) # Left, Right
levels(sample_info_csv$Handwriting) <- c(0,1)
sample_info_csv$Handwriting <- as.numeric(as.character(sample_info_csv$Handwriting))



###### Transform LOGICAL or Yes/No to numeric
for(i in 1:ncol(sample_info_csv)) {
  if(identical(sort(unique(sample_info_csv[, i])), c("No", "Yes"))) {
    sample_info_csv[, i] <- ifelse(sample_info_csv[, i] == "Yes", yes = 1, no = 0)
  }
  if(class(sample_info_csv[, i]) == "logical") sample_info_csv[, i] <- as.character(sample_info_csv[, i])
  
  if(identical(sort(unique(sample_info_csv[, i])), c("FALSE", "TRUE"))) {
    sample_info_csv[, i] <- ifelse(sample_info_csv[, i] == "TRUE", yes = 1, no = 0)
  }
  if(identical(sort(unique(sample_info_csv[, i])), c("FALSE"))) {
    sample_info_csv[, i][sample_info_csv[, i] == "FALSE"] <- 0
    sample_info_csv[, i] <- as.numeric(sample_info_csv[, i])
  }
}




### SPLITTING NON-ORDERED FACTORS into separate columns

### Why.Use.Glasses
unique(sample_info_csv$Why.Use.Glasses)
for(x in c("Myopia", "Presbyopia", "Hyperopia", "Fatigue", "Astigmatism")) {
  sample_info_csv[, x] <- as.numeric(grepl(x, sample_info_csv$Why.Use.Glasses))
}
sample_info_csv <- sample_info_csv[!colnames(sample_info_csv) %in% "Why.Use.Glasses"] # Remove original column

### Employment
for(x in  c("Medical leave", "Employee",  "Autonomous", "Unemployed", "Retired"  )) {
  sample_info_csv[, x] <- as.numeric(grepl(x, sample_info_csv$Employment))
}

sample_info_csv <- sample_info_csv[!colnames(sample_info_csv) %in% "Employment"] # Remove original column


### Employment_physical_office_travel
for(x in  c("Office", "Travel", "Physical"  )) {
  sample_info_csv[, x] <- as.numeric(grepl(x, sample_info_csv$Employment_physical_office_travel))
}
sample_info_csv <- sample_info_csv[!colnames(sample_info_csv) %in% "Employment_physical_office_travel"] # Remove original column

### Employment_category 
for(x in  c("Executive", "Salaried", "Own business"  )) {
  sample_info_csv[, x] <- as.numeric(grepl(x, sample_info_csv$Employment.Category))
}
sample_info_csv <- sample_info_csv[!colnames(sample_info_csv) %in% "Employment.Category"] # Remove original column

## Family_status
for(x in  c("Divorced/separated/widowed", "Married./.In.couple", "Single"  )) {
  sample_info_csv[, x] <- as.numeric(grepl(x, sample_info_csv$Family_status))
}
sample_info_csv <- sample_info_csv[!colnames(sample_info_csv) %in% "Family_status"] # Remove original column

## Exercice indoor-outdoor
for(x in  c("Indoor", "Outdoor")) {
  sample_info_csv[, x] <- as.numeric(grepl(x, sample_info_csv$Exercise..Indoor.Outdoor.))
}
sample_info_csv <- sample_info_csv[!colnames(sample_info_csv) %in% "Exercise..Indoor.Outdoor."] # Remove original column


## Blood group
sample_info_csv$Blood_group_A <- as.numeric(grepl("A", sample_info_csv$Blood.Group))
sample_info_csv$Blood_group_B <- as.numeric(grepl("B", sample_info_csv$Blood.Group))
sample_info_csv$Blood_group_A[is.na(sample_info_csv$Blood.Group)] <- NA
sample_info_csv$Blood_group_B[is.na(sample_info_csv$Blood.Group)] <- NA

sample_info_csv <- sample_info_csv[!colnames(sample_info_csv) %in% "Blood.Group"] # Remove original column



# I change NAs for median or for mode 

# But first I remove column "Last.cigarrete..years". Because to non-smokers I should actually asign an Inf. 
# I also remove Time married, because otherwise I have to put a number for married, single and widowed people. 
sample_info_csv <- sample_info_csv[, !colnames(sample_info_csv) %in% c("Time.Married...In.couple..years.", "Last.cigarrette..years.")]

# I change NAs for mode

mode_fun =function(x) {
  q=table(x)
  q=sort(q, decreasing =  TRUE)
  return(as.numeric(names(q[1])))
}


for(i in 3:ncol(sample_info_csv)) {
  sample_info_csv[, i][is.na(sample_info_csv[, i])] <- mode_fun(sample_info_csv[, i])
}




write.table(sample_info_csv, sep = ",", quote = F, row.names = F,
            file = "OUTPUT_lookalike_form_numeric.csv")





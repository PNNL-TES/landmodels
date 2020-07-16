#*****************************************************************************************************************
# This data contains 257 sites, 
# I want to know which factors (absLatitude, absLongitude, Biome, Ecosystem_type, Leaf_habit, MAT_Del, MAP_Del)
# can explain the MBE spatial variability, I used the random forest modeling, the results showed that MAP, absLatidue,
# MAT, abslongitude, and Ecosystem type are most important
# However, I am wondering a better way to do it and to visualization is using QRF, 
# Similar as the results in Figure 4 of your 2020 GBC paper
# In this case, I want MAP, absLatitude, MAT, and absLongitude in the a,b,c,d panels
# And in the bottom panel will be the Ecosystem types
# How long it will take for you to do this analysis, do you have time to do it?
#*****************************************************************************************************************

# Random forest model
data <- read.csv(here::here("outputs", "sub_res_site.csv"))

colnames(data)

# MBE: mean bias error
# absLatitude: absolute latitude, distance from equator
# absLongitude: absolute longitude, distance from 0 degree longitude 
# Biome: Biome type
# Ecosystem_type: ecosystem type
# Leaf_habit: everngreen, deciduous, mixed
# MAT_Del: mean annual temperature
# MAP_Del: mean annual precipitation
rf <- randomForest(MBE ~ absLatitude + absLongitude + Biome + Ecosystem_type + Leaf_habit + MAT_Del + MAP_Del,
                   data=data,
                   ntree = 100,
                   mtry = 2,
                   importance = TRUE,
                   proximity = TRUE)

summary(rf)
importance(rf, type = 1)
varImpPlot(rf)

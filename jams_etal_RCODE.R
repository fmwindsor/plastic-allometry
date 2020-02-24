## Supplementary Material - Data and code                                                            ##
## Jâms et al. (2020) Estimating the size distribution of plastics ingested by animals. Nat. Comms. ##
## Developed for R by Fred Windsor and Ifan Jâms                                                    ##
## Data exploration removed for conciseness, final analyses provided only                            ##
## For correspondence: Fred Windsor - fredric.windsor@ncl.ac.uk                                      ##

#### Software ####

# RStudio 1.1.463
# R 3.6.1 "Action of the Toes"


#### Packages and setup ####

# Check the specifications of R and RStudio (head off any potential versioning issues)
sessionInfo() # check the session information
search() # see which packages are attached

# Clear workspace
rm(list=ls())

# Load packages
library(tidyverse) # version 1.2.1
library(ggsci) # version 2.9
library(viridis) # version 0.5.1
library(Metrics) # version 0.1.4


#### Data input ####

# Read in the data file
dframe1 <- read.csv(file.choose()) # select "Jams_etal_DATA.csv"
names(dframe1) # check the column headings


#### Data analysis ####

### Modelling the allometric relationship

# Log10-log10 linear regression
model1 <- lm(log10(longest_ingested_plastic_mm) ~ log10(animal_body_length_mm), data = dframe1, na.action = na.omit)

# Model validation and results
plot(model1) # check the performance of the model
summary(model1) # get the results from the model

# Predictions for plotting
pdat <- expand.grid(animal_body_length_mm = seq(9,10340,0.5)) # Create a vector of numbers within the range of the animal size data
pred <- predict(model1, newdata = pdat, se.fit = T, level = 0.95, interval = "confidence") # Predict values for the vector using the model (95% CIs)
pred1 <- predict(model1, newdata = pdat, se.fit = T, level = 0.99, interval = "confidence") # Predict values for the vector using the model (99% CIs)
predframe <- data.frame(preds = pred$fit, se = pred$se.fit, animal_body_length_mm = pdat$animal_body_length_mm) # Store the size data in the same dataframe
predframe1 <- data.frame(preds = pred1$fit, se = pred1$se.fit, animal_body_length_mm = pdat$animal_body_length_mm) # Store the size data in the same dataframe

# Rename a column for plotting the n for each observation
dframe1<-rename(dframe1, n_pooledorganisms = longest.ingested.plastic.size.data_.is.for.how.many.pooled.organisms._i.e..the.number.of.organisms.that.contained.plastic.dicounting.those.that.did.not.include.plastic)
names(dframe1)

# Relevel the dataset
dframe1$group <- ordered(dframe1$group, levels = c("Invertebrates", "Reptiles","Mammals","Fishes"))

# Assign numeric values for the renamed column
dframe1$n_pooledorganisms <- as.numeric(as.character(dframe1$n_pooledorganisms))

# Figure 1a
plot_Fig1a <- ggplot(aes(x = log10(animal_body_length_mm), y = log10(longest_ingested_plastic_mm)), data = dframe1) +
  geom_ribbon(aes(x = log10(animal_body_length_mm), ymin = preds.lwr, ymax = preds.upr), fill = "black", inherit.aes = FALSE, alpha = 0.12, data = predframe) +
  geom_ribbon(aes(x = log10(animal_body_length_mm), ymin = preds.lwr, ymax = preds.upr), fill = "black", inherit.aes = FALSE, alpha = 0.12, data = predframe1) +
  scale_color_npg() + 
  geom_point(aes(color=group, size=n_pooledorganisms)) +
  scale_size(range = c(2.5, 11))+
  geom_line(aes(x = log10(animal_body_length_mm), y = preds.fit), data = predframe, cex = 0.8) +
  theme_bw() +  
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.background = element_rect(colour = "black", size = 0.6)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.box = "vertical", panel.grid = element_blank(), legend.background = element_rect(colour = "black", fill = 'NA')) +
  theme(text = element_text(size = 14), axis.text = element_text(colour = "black")) + 
  xlab("Animal body length (mm)") + 
  ylab("Longest ingested plastic (mm)") + 
  scale_x_continuous(labels = c(10,100,1000,10000), breaks = c(1,2,3,4)) +
  scale_y_continuous(labels = c(0.1,1,10,100,1000,10000), breaks = c(-1,0,1,2,3,4)) + 
  coord_cartesian(xlim = c(1,4), ylim = c(-1,3)) +
  ggtitle("") +
  annotation_logticks(scaled = TRUE, base = 10)
plot_Fig1a

# Export plot
ggsave("Fig1a.png", plot = last_plot(), scale = 1, width = 11, height = 14, units = c("cm"), dpi = 500, limitsize = TRUE)

# Get a map of the world for plotting
world <- map_data("world")

# Figure 1b 
plot_Fig1b <- ggplot() + 
  geom_polygon(aes(long, lat, group = group), data = world, size = 0.25, colour = "white", fill = "white") + 
  theme(text = element_text(size = 14), axis.text = element_text(colour = "black")) + 
  scale_color_npg() +
  geom_point(aes(longitude, latitude, color = group, size=n_pooledorganisms), data = dframe1) + 
  scale_size(range = c(2.5, 11))+
  theme(legend.position = "bottom", legend.title = element_blank(), legend.key=element_blank(), legend.box = "vertical", legend.background = element_rect(colour = "white", size = 0.4, fill = "white")) + 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.05), "cm"), axis.line = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + 
  theme(panel.background = element_rect(colour = "black", size = 1, fill = "lightgrey")) +
  ggtitle("") +  
  scale_x_continuous(breaks = c(-100,0,100)) + 
  scale_y_continuous(breaks = c(-50,0,50)) +
  ylab("Latitude") + 
  xlab("Longitude") + 
  coord_cartesian(xlim = c(-162, 175), ylim = c(-77, 84))
plot_Fig1b 

# Export plot
ggsave("Fig1b.png", plot = last_plot(), scale = 1, width = 6, height = 8.2, units = c("cm"), dpi = 500, limitsize =TRUE)


### Validating the allometric relationship

## Generating the predictions based on parameterisation data

# Set up the global variables
nperms = 1000 # number of random subsets of data
output = list(NULL) # setup a empty list to store the results

for (i in 1:nperms) { # sequentially select subsets of the data for analysis
  
  # Create subsets for parameterisation and validation 
  sub <- sample(nrow(dframe1),(round(nrow(dframe1)*0.1,0))) # randomly sample ~10% of the data
  val <- dframe1[sub,] # create a validation dataset (the validation data)
  obs <- dframe1[-sub,] # remove the data from the intial database (the parameterisation data)
  
  # Re-run the model for the subset of data
  lm <- lm(log10(longest_ingested_plastic_mm) ~ log10(animal_body_length_mm), data = obs) # model constructed based on the parameterisation data
  
  # Make predictions based on the model above
  pdat <- expand.grid(animal_body_length_mm = val$animal_body_length_mm) # set validation data for predictions 
  pred <- predict(lm, newdata = pdat, se.fit = T, level = 0.95, interval = "confidence") # use linear model to predict values for validation data
  predframe <- data.frame(preds = (10^pred$fit), se = pred$se.fit, animal.body.length.mm = pdat$animal_body_length_mm) # store results in a sensible way
  
  # Save the results of each iteration
  output[[i]] <- cbind(predframe[-4], val$longest_ingested_plastic_mm)
  
}

# Bind the results together into a dataset 
outputs <- bind_rows(output, .id = "ID")


## Statistical analysis of the predicted vs. observed data

# Root mean square of error for observations vs. predictions
rmse <- rmse(log10(outputs$`val$longest_ingested_plastic_mm`), log10(outputs$preds.fit))
rmse

# Compare the observed and predicted using a linear equation
comparison <- lm(log10(`val$longest_ingested_plastic_mm`) ~ log10(preds.fit), data = outputs)
plot(comparison)
summary.lm(comparison)
anova(comparison, test = "Chisq")

# Calculate the proportion of observed values falling within the confidence intervals of the equation
outputs$in_CIs <- rep(0, nrow(outputs)) # create a column to store the results
for (i in 1:nrow(outputs)) { 
  if (outputs[i,]$`val$longest_ingested_plastic_mm` > outputs[i,]$preds.lwr && outputs[i,]$`val$longest_ingested_plastic_mm` < outputs[i,]$preds.upr){ 
    outputs[i,]$in_CIs <- 1} # if the values fall within the confidence intervals then populate the record with a 1, if not the default is 0
}

# Proportion of the actual data that fall within the CIs of the predicted values
(sum(outputs$in_CIs)/length(outputs$in_CIs))*100


### Analysing the meta-data from the studies

## Depth range of sampled data

# Select data
dframe2 <- select(dframe1, species_binomial, group, max_depth_range_.m.)
dframe2 <- arrange(dframe2, species_binomial)
dframe2 <- distinct(dframe2, species_binomial, .keep_all= TRUE)
dframe2$max_depth_range_.m._numeric <- as.numeric(as.character(dframe2$max_depth_range_.m.))
dframe2 <- arrange(dframe2, max_depth_range_.m._numeric)

# Remove NAs from the datasets
dframe2<- dframe2 %>% drop_na(max_depth_range_.m._numeric)

# Add rows in for three global plastic models 
dframe2 <- add_row(dframe2, species_binomial = "Cózar et al. 2014", max_depth_range_.m._numeric = 25, group = "Global plastic distribution models")
dframe2<- add_row(dframe2, species_binomial = "Eriksen et al. 2014", max_depth_range_.m._numeric = 25, group = "Global plastic distribution models")
dframe2 <- add_row(dframe2, species_binomial = "van Sebille et al. 2015", max_depth_range_.m._numeric = 25, group = "Global plastic distribution models")

# Order the factor group for plotting 
dframe2$group <- ordered(dframe2$group, levels = c("Global plastic distribution models", "Invertebrates", "Mammals", "Reptiles", "Fishes"))
  
# Figure 2                                                
plot_Fig2 <- ggplot(dframe2, aes(x=reorder(species_binomial, -max_depth_range_.m._numeric), y=max_depth_range_.m._numeric, fill=group)) +
  geom_bar(stat="identity") +
  scale_y_reverse() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis(discrete = T, direction = -1) +
  labs(title = "", x = "", y = "Maximum depth range (m)", fill = "") + 
  theme(legend.position = c(0.705, 0.365)) +
  ggtitle("")
plot_Fig2

# Export plot
ggsave("Fig2.png", plot = last_plot(), scale = 1, width = 11, height = 11, units = c("cm"), dpi = 500, limitsize =TRUE)


## Magnification used in the studies

# Manipulate a new dataframe
dframe3<-select(dframe1, species_binomial, animal_body_length_mm, shortest_ingested_plastic_mm, magnification_used)
dframe3$magnification_used <- ordered(dframe3$magnification_used, levels = c( "Not reported or inferable", "Naked eye", "Stereomicroscope", "Stereomicroscope and compound microscope"))

# Predicting the relationship between the variables 
model2 <- lm(log10(shortest_ingested_plastic_mm) ~ log10(animal_body_length_mm), data= dframe3)
plot(model2)
summary(model2)

# Making predictions based on the model
pdat <- expand.grid(animal_body_length_mm = seq(9,10340,0.5))
pred <- predict(model2, newdata = pdat, se.fit = T, level = 0.95, interval = "confidence")
pred1 <- predict(model2, newdata = pdat, se.fit = T, level = 0.99, interval = "confidence")
predframe <- data.frame(preds = pred$fit, se = pred$se.fit, animal_body_length_mm = pdat$animal_body_length_mm)
predframe1 <- data.frame(preds = pred1$fit, se = pred1$se.fit, animal_body_length_mm = pdat$animal_body_length_mm)

# Figure 3
plot_Fig3 <- ggplot(aes(x = log10(animal_body_length_mm), y = log10(shortest_ingested_plastic_mm)), data = dframe3) +
  geom_ribbon(aes(x = log10(animal_body_length_mm), ymin = preds.lwr, ymax = preds.upr), data = predframe1, fill = "grey90", inherit.aes = FALSE) +
  geom_ribbon(aes(x = log10(animal_body_length_mm), ymin = preds.lwr, ymax = preds.upr), data = predframe, fill = "grey70", inherit.aes = FALSE) +
  geom_point(aes(x = log10(animal_body_length_mm), y = log10(shortest_ingested_plastic_mm), fill = magnification_used), cex = 3, pch = 21, colour = "black") +
  geom_line(aes(x = log10(animal_body_length_mm), y = preds.fit), data = predframe, cex = 0.8) +
  theme_bw() +  
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.background = element_rect(colour = "black", size = 0.6), panel.grid = element_blank()) +
  theme(legend.position = "", legend.title = element_blank(), legend.background = element_rect(colour = "black", fill = 'NA')) +
  theme(text = element_text(size = 14), axis.text = element_text(colour = "black")) + 
  xlab("Animal length (mm)") + 
  ylab("Smallest ingested plastic (mm)") + 
  scale_x_continuous(labels = c(10,100,1000,10000), breaks = c(1,2,3,4)) +
  scale_y_continuous(labels = c(0.1,1,10,100,1000), breaks = c(-1,0,1,2,3)) + 
  coord_cartesian(xlim = c(1,4), ylim = c(-1,3)) +
  scale_fill_viridis(discrete = T) +
  ggtitle("") +
  annotation_logticks(scaled = TRUE, base = 10)+
  theme_classic()+
  theme(legend.position = c(0.27, 0.885))+
  theme(legend.title = element_blank())
plot_Fig3

# Export plot
ggsave("Fig3.png", plot = last_plot(), scale = 1, width = 17, height = 15, units = c("cm"), dpi = 500, limitsize =TRUE)
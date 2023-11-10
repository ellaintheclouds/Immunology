#0 Import ----------------------------------------------------------------------
elisa_data <- read.csv("Data frame.csv", row.names = 1)
colnames(elisa_data) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                          "11", "12")

#1 Data analysis----------------------------------------------------------------
# Calculate the mean and standard deviation of the absorbency value for each of
# the standard concentrations and the patient samples in your data set and
# present the results in a table along with the original data. 

abs_mean <- apply(elisa_data, 2, mean)
abs_sd <- apply(elisa_data, 2, sd)

elisa_data[5 , ] <- abs_mean
elisa_data[7 , ] <- abs_sd

# Subtracting blank from mean
elisa_data[6, ] <- elisa_data[5, ] - elisa_data[5, 9]

row.names(elisa_data) <- c("A", "B", "C", "D", "mean", "mean_blank", "sd")

# Take some of this data frame to create a data frame for standard curve
standard_data <- elisa_data[1:7, 1:9]
patient_data <- elisa_data[1:7, 10:12]

# Add known concentrations to the standard data
standard_data[8,] <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0)
row.names(standard_data)[row.names(standard_data) == "8"] <- "conc"

#2 Plot-------------------------------------------------------------------------
# Plot the results from your standards on a graph of absorbance at 405nm against
# log antibody concentration (μg/ml)
# install.packages("ggplot2")
library(ggplot2)

# Adding log conc
#standard_data[6, 9] <- 1e-100 # Preventing issues with log(0)
standard_data[9, ] <- as.numeric(lapply(standard_data[8, ], log10))
row.names(standard_data)[row.names(standard_data) == "9"] <- "log_conc"

# Transforming data frame
curve_data <- standard_data[5:9, ]
curve_data_t <- data.frame(t(curve_data))
curve_data_t

curve_data_t$log_conc_restrict = ifelse(curve_data_t$log_conc < 1, 
                                        curve_data_t$log_conc, NA)

standard_curve <- 
  ggplot(curve_data_t, aes(log_conc, mean_blank)) +
  geom_point() + 
  geom_errorbar(aes(ymin = mean_blank - sd, ymax = mean_blank + sd, 
                    width = 0.1)) + 
  geom_smooth(aes(x = log_conc, y = mean_blank), formula = (y ~ exp(x)), 
              method = 'lm',se = FALSE, colour = "cornflowerblue", 
              fullrange = TRUE) + 
  scale_x_continuous(breaks = seq(-4.5, 3.5, by = 0.5)) + 
  scale_y_continuous(breaks = seq(0, 2, by = 0.25)) + 
  xlab("Log(IgG Concentration (μgml-1))") + ylab("Absorbance Units (405nm)") + 
  theme_bw()

ggsave("Exponential Standard Curve.png", plot = standard_curve, width = 10,
       height = 8)

#3 Curve modelling--------------------------------------------------------------
# Creating a dataframe with the required columns
df <- data.frame(concentrations = curve_data_t[, "conc"], # model dislikes log
                 measurements = curve_data_t[,"mean_blank"], 
                 sd = curve_data_t[, "sd"])

# Creating model
model4pl <- function(Concentration, Background, Mid, Slope, Bmax) {
  Bmax + ((Background - Bmax) / (1 + ((Concentration/Mid)^Slope)))
}

fit <- nls(measurements ~ model4pl(concentrations, Background, Mid, Slope, Bmax),
           data = df,
           start = c(Background=0, Mid=1.5, Slope=1, Bmax=1.75),
           control = nls.control(maxiter=1000, warnOnly=TRUE))
fit

cor(df$measurements, predict(fit))

#4 Re-plotting with modelled curve----------------------------------------------
modelled_curve <- ggplot(df, aes(concentrations, measurements)) +
  geom_point() + 
  geom_errorbar(aes(ymin = measurements - sd, ymax = measurements + sd, 
                    width = 0.1)) + 
  scale_x_log10() + scale_y_continuous(breaks = seq(0, 2, by = 0.25)) + 
  xlab("Log(IgG Concentration (μgml-1))") + ylab("Absorbance Units (405nm)") + 
  stat_function(data = df, fun  = model4pl,
                args = list(Mid = coef(fit)["Mid"],
                            Background = coef(fit)["Background"],
                            Slope = coef(fit)["Slope"],
                            Bmax = coef(fit)["Bmax"]), 
                colour = "cornflowerblue") + 
  theme_bw()

ggsave("Modelled Standard Curve.png", plot = modelled_curve, width = 10,
       height = 8)
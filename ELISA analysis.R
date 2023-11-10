#0 Setup------------------------------------------------------------------------
# Install packages and load libraries
# install.packages("ggplot2")
library(ggplot2)

# Import data obtained from ELISA practical
elisa_data <- read.csv("Data frame.csv", row.names = 1)

# Tidy column names
colnames(elisa_data) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                          "11", "12")

#1 Data analysis----------------------------------------------------------------
# Calculate mean and standard deviation of both control and patient samples

abs_mean <- apply(elisa_data, 2, mean)
abs_sd <- apply(elisa_data, 2, sd)

elisa_data[5 , ] <- abs_mean
elisa_data[7 , ] <- abs_sd

# Subtracting blank from mean
elisa_data[6, ] <- elisa_data[5, ] - elisa_data[5, 9]

# Providing row names
row.names(elisa_data) <- c("A", "B", "C", "D", "mean", "mean_blank", "sd")

# Creating a subset of data required for plotting standard curve
standard_data <- elisa_data[1:7, 1:9]
patient_data <- elisa_data[1:7, 10:12]

# Add known concentrations to the standard data
standard_data[8,] <- c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0)
row.names(standard_data)[row.names(standard_data) == "8"] <- "conc"

#2 Plot-------------------------------------------------------------------------
# Plot the results from your standards on a graph of absorbance at 405nm against
# log antibody concentration (μg/ml)

# Calculating log(concentration)
standard_data[9, ] <- as.numeric(lapply(standard_data[8, ], log10))
row.names(standard_data)[row.names(standard_data) == "9"] <- "log_conc"

# Switching rows and columns (to allow for plotting)
curve_data <- standard_data[5:9, ]
curve_data_t <- data.frame(t(curve_data))

# Ignore:(Archived code: to plot linear curve for restricted section)
#curve_data_t$log_conc_restrict = ifelse(curve_data_t$log_conc < 1, 
#                                        curve_data_t$log_conc, NA)

# Plotting a linear "standard curve"
standard_curve <- 
  ggplot(curve_data_t, aes(log_conc, mean_blank)) +
  geom_point() + 
  geom_errorbar(aes(ymin = mean_blank - sd, ymax = mean_blank + sd, 
                    width = 0.1)) + 
  geom_smooth(aes(x = log_conc, y = mean_blank), method = 'lm',se = FALSE, 
              colour = "cornflowerblue", fullrange = TRUE) + 
  scale_x_continuous(breaks = seq(-4.5, 3.5, by = 0.5)) + 
  scale_y_continuous(breaks = seq(0, 2, by = 0.25)) + 
  xlab("Log(IgG Concentration (μgml-1))") + ylab("Absorbance Units (405nm)") + 
  theme_bw()

#ggsave("Standard Curves/Standard Curve Linear.png", plot = standard_curve, 
# width = 10, height = 8)

#3 Curve modelling--------------------------------------------------------------
# Example curve modelling accessed on the following website: 
# https://jonathanrd.com/20-04-07-fitting-elisa-data-with-r/

# Creating another subset data frame with the required columns
df <- data.frame(concentrations = curve_data_t[, "conc"], # model dislikes log
                 measurements = curve_data_t[,"mean_blank"], 
                 sd = curve_data_t[, "sd"])

# Creating a model
model4pl <- function(Concentration, Background, Mid, Slope, Bmax) {
  Bmax + ((Background - Bmax) / (1 + ((Concentration/Mid)^Slope)))
}

# Iteration (took many attempts to find suitable start values)
fit <- nls(measurements ~ model4pl(concentrations, Background, Mid, Slope, Bmax),
           data = df,
           start = c(Background=0, Mid=1.5, Slope=1, Bmax=1.75),
           control = nls.control(maxiter=1000, warnOnly=TRUE))
print(fit)

cor(df$measurements, predict(fit))

#4 Re-plotting with modeled curve----------------------------------------------
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

#ggsave("Standard Curves/Modelled Standard Curve.png", plot = modelled_curve, 
#width = 10, height = 8)

#5 Calculating concentrations---------------------------------------------------
#(a mathematical method of "reading off the graph"
CalcConc <- function(Background, Mid, Slope, Bmax, y) {
  as.numeric(Mid * ((Background - Bmax)/(y - Bmax) - 1)^(1/Slope))
}

CalcConc( coef(fit)["Background"],
          coef(fit)["Mid"],
          coef(fit)["Slope"],
          coef(fit)["Bmax"],
          y = 1.40 )

# Putting these output into a dataframe
patient_calcultaions <- data.frame("Patient 1" = c(0.32, 0.1684709), 
                                   "Patient 2" = c(0.21, 0.027622), 
                                   "Patient 3" = c(1.40, 188.4006))

row.names(patient_calcultaions) <- c("Abs", "Log_Conc")

# Calculating anti-log to find concentration in μgmL-1
antilog <- function(x){10^x}
patient_calcultaions["Conc", ] <- apply(patient_calcultaions["Log_Conc", ], 2, 
                                        antilog)

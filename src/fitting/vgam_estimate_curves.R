# Usage:
#   ref_curve_generator.R input_file output_file degrees_of_freedom
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    print("No arguments supplied, use R CMD BATCH \"--args input_file='filename' output_file='filename' degrees_of_freedom=2 save_plot=1 plot_file='filename'\" ref_curve_generator.R")
    quit()
} else {
    for (i in 1:length(args)) {
        eval( parse(text = args[[i]]) )
    }
}

print(input_file)
print(output_file)
print(degrees_of_freedom)
print(plot_file)
print(densities_file)

library(VGAM)

#
# Read the data from the csv file given as first argument.
#
table_all <- read.csv(file = input_file, head = TRUE, sep = ",")
table <- table_all[c(2,3)]

#
# Perform the fit. The degrees of freedom are given by the third argument.
#
#fit <- vgam(value ~ s(progress, df = as.integer(degrees_of_freedom)), lms.yjn2(percentiles = c(5,25,50,75,95), lsigma = 'identity'), table, trace = TRUE)
fit <- vgam(value ~ s(progress, df = as.integer(degrees_of_freedom)), lms.yjn2(percentiles = c(1,5,25,50,75,95,99)), table, trace = TRUE)
#fit <- vgam(value ~ s(progress, df = c(as.integer(degrees_of_freedom), 1)), lms.yjn2(percentiles = c(5,25,50,75,95), lsigma = 'identity', dfmu.init = 2) , table, trace = TRUE)
#fit <- vgam(value ~ s(progress, df = 2), lms.bcn(percentiles = c(5,25,50,75,95), dfmu.init = 2, lmu = 'identity', zero = c(1,3)), table, trace = TRUE)
#fit <- vgam(value ~ s(progress, df = c(as.integer(degrees_of_freedom), 2)), lms.bcn(percentiles=c(5,25,50,75,95), dfmu.init = 2), table, trace = TRUE)
#fit <- vgam(value ~ s(progress, df = c(3,1)), lms.bcn(percentiles = c(5,25,50,75,95), zero = 1), table, trace = TRUE)

#
# Concatenate the age, lms coeffs and fitted percentiles.
#
fit_frame <- data.frame(fit@x, predict(fit), fit@fitted.values, fit@misc$yoffset)
#fit_frame <- data.frame(fit@x, predict(fit), fit@fitted.values, cdf(fit), table)

#
# Write the output data to csv file with filename given by the second argument
#
df <- write.csv(fit_frame, file = output_file)

#
# Save probability densities to csv file
#
min_progression <- -42
max_progression <- 51
min_value <- min(table[,2])
max_value <- max(table[,2])
value_offset <- (max_value - min_value) * 0.2
value_step <- (max_value - min_value) / 250
values <- seq( min_value - value_offset, max_value+value_offset, by = value_step )
values_frame <- data.frame(values)
for (prog in seq(min_progression, max_progression)) {
	deplot(fit, x0 = prog, y = values, show.plot = FALSE) -> aa
	values_frame <- data.frame( values_frame, aa@post$deplot$density )
}
df <- write.csv(t(values_frame), file = densities_file)


#
# Save plot
#
w <- 12
h <- 8
if (save_plot == 0) {
    x11(width = w, height = h)
} else {
    pdf(plot_file, width = w, height = h)
}
qtplot(fit, percentiles = c(5, 25, 50, 75, 95), main = "Quantiles",xlim = c(-40, 40), las = 1, ylab = "Metric value", lwd = 2, lcol = 4)

if (save_plot == 0) {
    z <- locator()
} else {
    dev.off()
}

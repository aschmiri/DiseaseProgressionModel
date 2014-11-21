# Usage:
#   ref_curve_generator.R input_file output_file degrees_of_freedom
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    print("No arguments supplied, use R CMD BATCH \"--args input_file='filename' output_file='filename' densities_file='filename' degrees_of_freedom=2 plot_file='filename'\" ref_curve_generator.R")
    quit()
} else {
    for (i in 1:length(args)) {
        eval(parse(text = args[[i]]))
    }
}

print(input_file)
print(output_file)
print(degrees_of_freedom)
print(zero)

library(VGAM)
sessionInfo()
#packageVersion("VGAM")

#
# Read the data from the csv file given as first argument.
#
table_all <- read.csv(file = input_file, head = TRUE, sep = ",")
table <- table_all[c(2,3)]

#
# Perform the fit. The degrees of freedom are given by the third argument.
#
if (zero=="None") {
    fit <- vgam(value ~ s(progress, df = as.integer(degrees_of_freedom)), lms.yjn2(), data = table, trace = TRUE)
} else {
    fit <- vgam(value ~ s(progress, df = as.integer(degrees_of_freedom)), lms.yjn2(zero = as.integer(zero)), data = table, trace = TRUE)
}

#
# Concatenate the age, lms coeffs and fitted percentiles.
#
fit_frame <- data.frame(fit@x, predict(fit), fit@fitted.values, fit@misc$yoffset)

#
# Write the output data to csv file with filename given by the second argument
#
df <- write.csv(fit_frame, file = output_file)

#
# Save plot
#
#plot_file <- "/vol/medic01/users/aschmidt/projects/DiseaseProgressionModel/models/plot.pdf"
#save_plot <- 1
#min_progression <- min(table[,1])
#max_progression <- max(table[,1])
##if (exists(plot_file)) {
#	print(plot_file)
#	
#	w <- 12
#	h <- 8
#	if (save_plot == 0) {
#	    x11(width = w, height = h)
#	} else {
#	    pdf(plot_file, width = w, height = h)
#	}
#	qtplot(fit, percentiles = c(5, 25, 50, 75, 95), main = "Quantiles", xlim = c(min_progression, max_progression), las = 1, ylab = "Metric value", lwd = 2, lcol = 4)
#	
#	if (save_plot == 0) {
#	    z <- locator()
#	} else {
#	    dev.off()
#	}
##}

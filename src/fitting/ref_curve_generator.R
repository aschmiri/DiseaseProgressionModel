# Usage:
#   ref_curve_generator.R input_file output_file degrees_of_freedom
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0){
    print("No arguments supplied, use R CMD BATCH \"--args input_file='filename' output_file='filename' degrees_of_freedom=2 save_plot=1 plot_filename='filename'\" ref_curve_generator.R")
    quit()
} else {
    for (i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}

print(input_file)
print(output_file)
print(degrees_of_freedom)
print(plot_filename)

library(VGAM)

# Read the data from the csv file given as first argument.
bv <- read.csv(file=input_file,head=TRUE,sep=",")
bv_mod <- bv[c(2,3)];

# Perform the fit. The degrees of freedom are given by the third argument.
fit = vgam(volume ~ s(scan_age, df=as.integer(degrees_of_freedom)), lms.yjn2(percentiles=c(5,25,50,75,95),lsigma='identity'), bv_mod, trace=TRUE)
#fit = vgam(volume ~ s(scan_age, df=c(as.integer(degrees_of_freedom),1)), lms.yjn2(percentiles=c(5,25,50,75,95),lsigma='identity',dfmu.init=2), bv_mod, trace=TRUE)
#fit = vgam(volume ~ s(scan_age, df=2), lms.bcn(percentiles=c(5,25,50,75,95),dfmu.init=2,lmu='identity',zero=c(1,3)), bv_mod, trace=TRUE)
#fit = vgam(volume ~ s(scan_age, df=c(as.integer(degrees_of_freedom),2)), lms.bcn(percentiles=c(5,25,50,75,95), dfmu.init = 2), bv_mod, trace=TRUE)
#fit = vgam(volume ~ s(scan_age, df=c(3,1)), lms.bcn(percentiles=c(5,25,50,75,95),zero = 1), bv_mod, trace=TRUE)

# Concatenate the age, lms coeffs and fitted percentiles.
fit_c <- data.frame(fit@x, predict(fit), fit@fitted.values, fit@misc$yoffset)
#fit_c <- data.frame(fit@x, predict(fit), fit@fitted.values, cdf(fit), bv_mod)

# Write the output data to csv file with filename given by the second argument
df <- write.csv(fit_c, file=output_file)

w = 12
h = 8
if(save_plot==0) {
    x11(width=w,height=h)
} else {
    pdf(plot_filename,width=w,height=h)
}
qtplot(fit, percentiles = c(5, 25, 50, 75, 95), main = "Quantiles",xlim = c(43, 100), las = 1, ylab = "volume", lwd = 2, lcol = 4)
if(save_plot==0) {
    z <- locator()
} else {
    dev.off()
}


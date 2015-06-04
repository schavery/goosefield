# R for plotting the results.
d <- read.table("results.csv", sep=",", header = TRUE)
fft <- d[d$type == "fft",]
naive <- d[d$type == "naïve",]
hart <- d[d$type == "hartley",]
foil <- d[d$type == "foil",]

# do the degree graph
x <- c(fft$degree)
y_fft <- c(fft$time)
y_naive <- c(naive$time)
y_hart <- c(hart$time)
y_foil <- c(foil$time)

plot(x,	y_fft, log="xy", xlab="Degree",
	ylab="Time (s)", main="Polynomial Multiplication Timing",
	ylim=c(min(y_foil),max(y_naive)), col=3, pch=1,
	sub="Using cached M*")

points(x, y_naive, col=4, pch=3)

# reset x
x <- c(foil$degree)
points(x, y_hart, col=5, pch=4) #?
points(x, y_foil, col=2, pch=2)

# need to adjust position
legend(10, 250, legend=c("FFT", "Naïve", "FOIL", "Hartley"),
	col=c(3,4,2,5), pch=c(1,3,2,4))

# this might work
dev.copy(png, "degrees.png")
dev.off()

# now do the sparsity graph
sfft <- foil[fft$degree == 1000,]
snaive <- naive[naive$degree == 1000,]
sfoil <- foil[foil$degree == 1000,]
shart <- hart[hart$degree == 1000,] # heh

x <- rev(c(sfft$sparsity))
y_fft <- rev(c(sfft$time))
y_naive <- rev(c(snaive$time))
y_hart <- rev(c(shart$time))
y_foil <- rev(c(sfoil$time))

plot(x, y_fft, log="x", xlab="Density (% non-zero elements)",
	ylab="Time (s)", main="Polynomial Multiplication Timing",
	ylim=c(min(y_foil), max(y_naive)), col=3,
	sub="Degree = 1000; Using cached M*")

points(x, y_foil, col=2, pch=2)
points(x, y_naive, col=4, pch=3)
points(x, y_hart, col=5, pch=4) #?

legend(2, 55, legend=c("FFT", "Naïve", "FOIL", "Hartley"),
	col=c(3,4,2,5), pch=c(1,3,2,4))

dev.copy(png, "sparse.png")
dev.off()

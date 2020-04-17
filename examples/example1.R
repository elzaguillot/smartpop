data=read.table('example1all.txt',h=T)

png('example1.png',width=500,height=500)
boxplot(data$mt.thetapi~as.factor(data$popsize),xlab='Population size',ylab=expression(paste(theta[pi])))
dev.off()

png('example1bis.png',width=500,height=500)
boxplot(data$mt.thetawatterson~as.factor(data$popsize),xlab='Population size',ylab=expression(paste(theta[Watterson])))
dev.off()

png('example1_hap.png',width=500,height=500)
boxplot(data$mt.nbhap~as.factor(data$popsize),xlab='Population size',ylab='Number of haplotypes')
dev.off()

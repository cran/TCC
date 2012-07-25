
MAplot <- function(datalist, FDR_threshold=0.01){
  data <- datalist$counts
  data.cl <- datalist$group
  norm_f_RPM <- 1e6/colSums(data)
  norm_f_TbT <- datalist$norm_f_TbT
  RPM <- sweep(data, 2, norm_f_RPM, "*")
  data <- RPM
  meanA <- log2(apply(data[,data.cl==1], 1, mean))
  meanB <- log2(apply(data[,data.cl==2], 1, mean))
  x_axis <- (meanA + meanB)/2
  y_axis <- meanB - meanA 
  plot(x_axis, y_axis, xlab="A = (log2(B)+log2(A))/2", ylab="M = log2(B)-log2(A)", pch=20, cex=.3)
  grid(col="gray", lty="dotted")
  points(x_axis[datalist$data$FDR < FDR_threshold], y_axis[datalist$data$FDR< FDR_threshold], col=2, pch=20, cex=0.3)
  baseline_TbT <- log2(mean(norm_f_TbT[data.cl==2])/mean(norm_f_TbT[data.cl==1]))
  abline(h=baseline_TbT, col="red", lwd=1)
}

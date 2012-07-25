
MAplot <- function(datalist, FDR_threshold=0.01){
  data <- datalist$counts
  data.cl <- datalist$group
  norm_f_TbT <- datalist$norm_f_TbT
  x_axis <- datalist$Mval
  y_axis <- datalist$Aval
  plot(x_axis, y_axis, xlab="A = (log2(B)+log2(A))/2", ylab="M = log2(B)-log2(A)", pch=20, cex=.3)
  grid(col="gray", lty="dotted")
  points(x_axis[datalist$data$FDR < FDR_threshold], y_axis[datalist$data$FDR< FDR_threshold], col=2, pch=20, cex=0.3)
  baseline_TbT <- log2(mean(norm_f_TbT[data.cl==2])/mean(norm_f_TbT[data.cl==1]))
  abline(h=baseline_TbT, col="red", lwd=1)
}

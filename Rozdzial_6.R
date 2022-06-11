##################################
#METODA NAJWIEKSZEJ WIAROGODNOSCI#
##################################

Newtons_method <- function(data, derivative, s_derivative, start_point, epsilon = 0.0001) { #metoda newtona
  i = 0
  best_est <- start_point(data)
  while(abs(derivative(best_est, data) / s_derivative(best_est)) >= epsilon) {
    best_est <- best_est - derivative(best_est, data) / s_derivative(best_est)
    if(i == 15) {
      print("Brak zbieznosci!")
      return(NA)
    }
    i = i + 1
  }
  return(best_est)
}

alfa_value <- function(data) {
  return(- log(mean(data)) + 1/length(data)*sum(log(data)))
}

derivative_of_l <- function(x, data) {
  return(psigamma(x, deriv = 0) - log(x) - alfa_value(data))
}

second_derivative_of_l <- function(x) {
  return(psigamma(x, deriv = 1) - 1/x)
}


first_est_of_k <- function(data) { 
  return( (-3 - sqrt(9-12*alfa_value(data)))/(12*alfa_value(data)) )
}

est_of_theta_mle <- function(data, k) {
  return(mean(data)/k)
}

#################
#METODA MOMENTOW#
#################

first_moment_est <- function(data) {
  return(mean(data))
} 

second_moment_est <- function(data) {
  return(1/length(data)*sum(data^2))
}

k_mom_est <- function(data) {
  return(first_moment_est(data)^2/(second_moment_est(data) - first_moment_est(data)^2))
}

theta_mom_est <- function(data) {
  return((second_moment_est(data) - first_moment_est(data)^2)/first_moment_est(data))
}

############################
#ESTYMATORY BIAS, VAR I MSE#
############################

generate <- function(shape, scale = 1, trial_size = 100, n = 10000, epsilon = 10^(-6)) {
  trials <- lapply(1:n, function(x) rgamma(trial_size, shape = shape, rate = 1/scale))
  mle_estimators_shape <- sapply(1:n, function(x) Newtons_method(trials[[x]], derivative_of_l, second_derivative_of_l, first_est_of_k, epsilon))
  mle_estimators_scale <- sapply(1:n, function(x) est_of_theta_mle(trials[[x]], mle_estimators_shape[x]))
  mom_estimators_shape <- sapply(1:n, function(x) k_mom_est(trials[[x]]))
  mom_estimators_scale <- sapply(1:n, function(x) theta_mom_est(trials[[x]]))
  
  Estimators_vec <- c(mean(mle_estimators_shape), mean(mom_estimators_shape), mean(mle_estimators_scale), mean(mom_estimators_scale))
  Estimators <- matrix(Estimators_vec, 2, 2)
  colnames(Estimators) <- c("shape", "scale")
  rownames(Estimators) <- c("mle", "mom")
  
  biases_vec <- c(1/n*sum(mle_estimators_shape - shape), 1/n*sum(mom_estimators_shape - shape), 1/n*sum(mle_estimators_scale - scale), 1/n*sum(mom_estimators_scale - scale))
  biases <- matrix(biases_vec, 2, 2)
  colnames(biases) <- c("shape", "scale")
  rownames(biases) <- c("mle", "mom")
  
  MSE_vec <- c(1/n*sum((mle_estimators_shape - shape)^2), 1/n*sum((mom_estimators_shape - shape)^2), 1/n*sum((mle_estimators_scale - scale)^2), 1/n*sum((mom_estimators_scale - scale)^2))
  MSEs <- matrix(MSE_vec, 2, 2)
  colnames(MSEs) <- c("shape", "scale")
  rownames(MSEs) <- c("mle", "mom")
  
  VAR_vec <- c(var(mle_estimators_shape), var(mom_estimators_shape), var(mle_estimators_scale), var(mom_estimators_scale))
  Variances <- matrix(VAR_vec, 2, 2)
  colnames(Variances) <- c("shape", "scale")
  rownames(Variances) <- c("mle", "mom")
  
  MSEs_est <- biases^2 + Variances
  
  distance_shape <- sqrt(sum((mle_estimators_shape - mom_estimators_shape)^2))
  distance_scale <- sqrt(sum((mle_estimators_scale - mom_estimators_scale)^2)) 
  
  return(list("Estimators" = Estimators, "Biases" = biases, "MSEs" = MSEs, "Variances" = Variances, "MSEs_est" = MSEs_est, "shape_distance" = distance_shape, "scale_distance" = distance_scale))
}


summary <- function(trial_size = 100, shape = 5, scale = 1, size = 10000, epsilon = 10^(-6)) {
  data <- generate(shape, scale, trial_size, size, epsilon)
  cat("--- SHAPE:", shape, "\t SCALE:", scale," ---\n")
  cat("Mean of Estimators:\n")
  print(data$Estimators)
  cat("\n Biases:\n")
  print(data$Biases)
  cat("\nVariances:\n")
  print(data$Variances)
  cat("\nMSEs:\n")
  print(data$MSEs)
  cat("\nMSEs estimated from variances and biases:\n")
  print(data$MSEs_est)
  cat("\nDistance between shape estimators in L^2-norm: ")
  cat(data$shape_distance)
  cat("\nDistance between scale estimators in L^2-norm: ")
  cat(data$scale_distance)
}

#########
#WYKRESY#
#########
library("ggplot2")
library("ggpubr")

N <- 10000
n <- 50

shape1 <- c(0.5, 1, 5, 10, 15, 20)
scale1 <- 1
estimators1 <- lapply(shape1, function(x) generate(shape = x, scale = scale1, trial_size = n, n = N))
MLE_bias1 <- sapply(estimators1, function(x) x$Biases['mle','shape'])
MLE_MSE1 <- sapply(estimators1, function(x) x$MSEs['mle','shape'])
MLE_var1 <- sapply(estimators1, function(x) x$Variances['mle','shape'])
MOM_bias1 <- sapply(estimators1, function(x) x$Biases['mom','shape'])
MOM_MSE1 <- sapply(estimators1, function(x) x$MSEs['mom','shape'])
MOM_var1 <- sapply(estimators1, function(x) x$Variances['mom','shape'])
data1 <- data.frame("kształt" = c(shape1, shape1), Estymator = factor(c(rep("ENW", length(estimators1)), rep("Metoda momentów", length(estimators1)))), "bias" = c(MLE_bias1, MOM_bias1), "MSE" = c(MLE_MSE1, MOM_MSE1), "Var" = c(MLE_var1, MOM_var1))

plot1_1 <- ggplot(data1, aes(x = kształt, y = bias, colour = Estymator)) +
      geom_line(size = 0.5) + 
      geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
      scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
      theme(legend.position = "none")

plot1_2 <- ggplot(data1, aes(x = kształt, y = MSE, colour = Estymator)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  theme(legend.position = "none")

plot1_3 <- ggplot(data1, aes(x = kształt, y = Var, colour = Estymator)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9"))

ggarrange(plot1_1, plot1_2, plot1_3, ncol=3, common.legend = TRUE, legend = "top") #RYSUNEK 5#

###############

shape2 <- c(0.5, 3, 100)
scale2 <- c(0.5, 1, 5, 10, 15, 25, 30, 35, 40)
parameters2 <- rbind(rep(scale2, 3), rep(shape2, each = 9))
estimators2 <- lapply(1:27, function(x) generate(shape = parameters2[2,x], scale = parameters2[1,x], trial_size = n, n = N))
MLE_bias2 <- sapply(estimators2, function(x) x$Biases['mle','scale'])
MLE_MSE2 <- sapply(estimators2, function(x) x$MSEs['mle','scale'])
MLE_var2 <- sapply(estimators2, function(x) x$Variances['mle','scale'])
MOM_bias2 <- sapply(estimators2, function(x) x$Biases['mom','scale'])
MOM_MSE2 <- sapply(estimators2, function(x) x$MSEs['mom','scale'])
MOM_var2 <- sapply(estimators2, function(x) x$Variances['mom','scale'])
data2 <- data.frame("skala" = rep(c(scale2, scale2, scale2),3),
                    "Kształt" = factor(rep(parameters2[2,],3)),
                    Estymator = factor(rep(c(rep("ENW", length(estimators2)), rep("Metoda momentów", length(estimators2))),3)),
                    "bias" = c(MLE_bias2, MOM_bias2),
                    "MSE" = c(MLE_MSE2, MOM_MSE2),
                    "Var" = c(MLE_var2, MOM_var2))

plot2_1 <- ggplot(data2, aes(x = skala, y = bias, colour = Estymator, linetype = Kształt)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  theme(legend.position = "none") +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"))

plot2_2 <- ggplot(data2, aes(x = skala, y = MSE, colour = Estymator, linetype = Kształt)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  theme(legend.position = "none") +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"))

plot2_3 <- ggplot(data2, aes(x = skala, y = Var, colour = Estymator, linetype = Kształt)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"))

ggarrange(plot2_1, plot2_2, plot2_3, ncol=3, common.legend = TRUE, legend = "top") #RYSUNEK 6#

##############

n <- 100

shape3 <- c(0.5, 1, 5, 10, 15, 20)
scale3 <- 1
estimators3 <- lapply(shape3, function(x) generate(shape = x, scale = scale3, trial_size = n, n = N))
MLE_bias3 <- sapply(estimators3, function(x) x$Biases['mle','shape'])
MLE_MSE3 <- sapply(estimators3, function(x) x$MSEs['mle','shape'])
MLE_var3 <- sapply(estimators3, function(x) x$Variances['mle','shape'])
MOM_bias3 <- sapply(estimators3, function(x) x$Biases['mom','shape'])
MOM_MSE3 <- sapply(estimators3, function(x) x$MSEs['mom','shape'])
MOM_var3 <- sapply(estimators3, function(x) x$Variances['mom','shape'])
data3 <- data.frame("kształt" = c(shape3, shape3), Estymator = factor(c(rep("ENW", length(estimators3)), rep("Metoda momentów", length(estimators3)))), "bias" = c(MLE_bias3, MOM_bias3), "MSE" = c(MLE_MSE3, MOM_MSE3), "Var" = c(MLE_var3, MOM_var3))

plot3_1 <- ggplot(data3, aes(x = kształt, y = bias, colour = Estymator)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  theme(legend.position = "none")

plot3_2 <- ggplot(data3, aes(x = kształt, y = MSE, colour = Estymator)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  theme(legend.position = "none")

plot3_3 <- ggplot(data3, aes(x = kształt, y = Var, colour = Estymator)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9"))

ggarrange(plot3_1, plot3_2, plot3_3, ncol=3, common.legend = TRUE, legend = "top") #RYSUNEK 7#

#################

shape4 <- c(0.5, 3, 100)
scale4 <- c(0.5, 1, 5, 10, 15, 25, 30, 35, 40)
parameters4 <- rbind(rep(scale4, 3), rep(shape4, each = 9))
estimators4 <- lapply(1:27, function(x) generate(shape = parameters4[2,x], scale = parameters4[1,x], trial_size = n, n = N))
MLE_bias4 <- sapply(estimators4, function(x) x$Biases['mle','scale'])
MLE_MSE4 <- sapply(estimators4, function(x) x$MSEs['mle','scale'])
MLE_var4 <- sapply(estimators4, function(x) x$Variances['mle','scale'])
MOM_bias4 <- sapply(estimators4, function(x) x$Biases['mom','scale'])
MOM_MSE4 <- sapply(estimators4, function(x) x$MSEs['mom','scale'])
MOM_var4 <- sapply(estimators4, function(x) x$Variances['mom','scale'])
data4 <- data.frame("skala" = rep(c(scale4, scale4, scale4),3),
                    "Kształt" = factor(rep(parameters4[2,],3)),
                    Estymator = factor(rep(c(rep("ENW", length(estimators4)), rep("Metoda momentów", length(estimators4))),3)),
                    "bias" = c(MLE_bias4, MOM_bias4),
                    "MSE" = c(MLE_MSE4, MOM_MSE4),
                    "Var" = c(MLE_var4, MOM_var4))

plot4_1 <- ggplot(data4, aes(x = skala, y = bias, colour = Estymator, linetype = Kształt)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  theme(legend.position = "none") +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"))

plot4_2 <- ggplot(data4, aes(x = skala, y = MSE, colour = Estymator, linetype = Kształt)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  theme(legend.position = "none") +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"))

plot4_3 <- ggplot(data4, aes(x = skala, y = Var, colour = Estymator, linetype = Kształt)) +
  geom_line(size = 0.5) + 
  geom_point(shape=21, color="black", fill="#f39730", size = 1.5) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#c10b4b", "#0c9fb9")) +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"))

ggarrange(plot4_1, plot4_2, plot4_3, ncol=3, common.legend = TRUE, legend = "top") #RYSUNEK 8#

###################
#ANALIZA WARIANCJI#
###################

var_k_MLE <- function(k) {
  return(k/(k*psigamma(k, deriv = 1) - 1))
}

var_t_MLE <- function(t, k) {
  return(psigamma(k, deriv = 1)*t^2/(k*psigamma(k, deriv = 1) - 1))
}

var_k_MOM <- function(k) {
  return(2*k*(k+1))
}

var_t_MOM <- function(k, t) {
  return(t^2*(3/k + 2))
}

#########
#WYKRESY#
#########

trial_size1 <- 200
domain1 <- seq(0.0001,150, 0.001)
variance_k <- var_k_MLE(domain1)
estimated_var_k <- trial_size1*sapply(seq(0.1,100.1,10), function(x) generate(x, trial_size = trial_size1)$Variances[1,1])
data1 <- data.frame("Kształt" = domain1,
                    "Wariancja" = variance_k)

plot5_1 <- ggplot(data1, aes(x = Kształt, y = Wariancja)) +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_k, col = 'black', size = 2, shape = 21, fill = "red") +
  geom_line(size = 0.6, col = "red") +
  coord_cartesian(xlim =c(0, 100), ylim = c(0, 20000))

plot5_1 #RYSUNEK 9#

###############

trial_size2 <- 200
shapes2 <- c(0.1, 0.5, 1, 10)
domain2 <- seq(0.0001, 150, 0.001)
variance_t <- as.vector(sapply(shapes2, function(x) var_t_MLE(domain2, x)))
estimated_var_t <- trial_size2*sapply(shapes2, function(y) sapply(seq(0.1,100.1,10), function(x) generate(shape = y, scale = x, trial_size = trial_size2)$Variances[1,2]))
data2 <- data.frame("Skala" = rep(domain2, length(shapes2)),
                    "Wariancja" = variance_t,
                    "Kształt" = factor(rep(shapes2, each = length(domain2))))

plot6_1 <- ggplot(data2, aes(x = Skala, y = Wariancja, colour = Kształt)) +
  geom_line(size = 0.6) +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,1], col = 'black', size = 2, shape = 21, fill = "#F8766D") +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,2], col = 'black', size = 2, shape = 21, fill = "#7CAE00") +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,3], col = 'black', size = 2, shape = 21, fill = "#00BFC4") +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,4], col = 'black', size = 2, shape = 21, fill = "#C77CFF") +
  coord_cartesian(xlim =c(0, 100), ylim = c(0, 119000))

plot6_1 #RYSUNEK 10#

##########

trial_size1 <- 200
domain1 <- seq(0.0001,150, 0.001)
variance_k <- var_k_MOM(domain1)
estimated_var_k <- trial_size1*sapply(seq(0.1,100.1,10), function(x) generate(x, trial_size = trial_size1)$Variances[2,1])
data1 <- data.frame("Kształt" = domain1,
                    "Wariancja" = variance_k)

plot7_1 <- ggplot(data1, aes(x = Kształt, y = Wariancja)) +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_k, col = 'black', size = 2, shape = 21, fill = "red") +
  geom_line(size = 0.6, col = "red") +
  coord_cartesian(xlim =c(0, 100), ylim = c(0, 20000))

plot7_1 #RYSUNEK 11#

###############

trial_size2 <- 200
shapes2 <- c(0.1, 0.5, 1, 10)
domain2 <- seq(0.0001, 150, 0.001)
variance_t <- as.vector(sapply(shapes2, function(x) var_t_MOM(x,domain2)))
estimated_var_t <- trial_size2*sapply(shapes2, function(y) sapply(seq(0.1,100.1,10), function(x) generate(shape = y, scale = x, trial_size = trial_size2)$Variances[2,2]))
data2 <- data.frame("Skala" = rep(domain2, length(shapes2)),
                    "Wariancja" = variance_t,
                    "Kształt" = factor(rep(shapes2, each = length(domain2))))

plot8_1 <- ggplot(data2, aes(x = Skala, y = Wariancja, colour = Kształt)) +
  geom_line(size = 0.6) +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,1], col = 'black', size = 2, shape = 21, fill = "#F8766D") +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,2], col = 'black', size = 2, shape = 21, fill = "#7CAE00") +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,3], col = 'black', size = 2, shape = 21, fill = "#00BFC4") +
  annotate("point", x = seq(0.1,100.1,10), y = estimated_var_t[,4], col = 'black', size = 2, shape = 21, fill = "#C77CFF") +
  coord_cartesian(xlim =c(0, 100), ylim = c(0, 310000))

plot8_1 #RYSUNEK 12#

###############

domain1 <- seq(0.0001,150, 0.001)
variance_k_mom <- var_k_MOM(domain1)
variance_k_mle <- var_k_MLE(domain1)
data1 <- data.frame("Kształt" = rep(domain1,2),
                    "Wariancja" = c(variance_k_mle, variance_k_mom),
                    "Estymator" = factor(c(rep("ENW", length(domain1)), rep("Metoda momentów", length(domain1)))))

plot9_1 <- ggplot(data1, aes(x = Kształt, y = Wariancja, colour = Estymator)) +
  geom_line(size = 0.6) +
  coord_cartesian(xlim =c(0, 60), ylim = c(0, 8000)) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#d41aaf", "#2b70d4"))

plot9_1 #RYSUNEK 13#

###############

domain2 <- seq(0.0001,50, 0.1)
shapes2 <- c(0.1, 0.5, 20)
variance_t_mom <- as.vector(sapply(shapes2, function(x) var_t_MOM(x, domain2)))
variance_t_mle <- as.vector(sapply(shapes2, function(x) var_t_MLE(domain2, x)))
data2 <- data.frame("Skala" = rep(domain2,6),
                    "Wariancja" = c(variance_t_mle, variance_t_mom),
                    "Estymator" = factor(c(rep("ENW", length(domain2)*3), rep("Metoda momentów", length(domain2)*3))),
                    "Kształt" = factor(c(rep(rep(shapes2, each = length(domain2)),2))))

plot10_1 <- ggplot(data2, aes(x = Skala, y = Wariancja, colour = Estymator, linetype = Kształt)) +
  geom_line(size = 0.7) +
  coord_cartesian(xlim =c(0, 25), ylim = c(0, 20000)) +
  scale_color_manual(breaks = c("ENW","Metoda momentów"), values = c("#d41aaf", "#2b70d4")) +
  scale_linetype_manual(values=c("twodash", "dotted", "solid"))

plot10_1 #RYSUNEK 14#
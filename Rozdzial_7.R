library("ggplot2")
library("ggpubr")

##############################################
#FUNKCJE DO ESTYMACJI PARAMETROW POCZATKOWYCH#
##############################################
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

first_est_of_k <- function(data) { 
  return( (-3 - sqrt(9-12*alfa_value(data)))/(12*alfa_value(data)) )
}

est_of_theta_mle <- function(data, k) {
  return(mean(data)/k)
}
##########################################
#FUNKCJE DO ESTYMACJI  PARAMETRU KSZTALTU#
##########################################

S_fun <- function(X, a) {
  log(sum(a)) - log(sum(X*a)) + sum(a*log(X))/sum(a)
}

estimate_of_shape_fun <- function(X, a) {
  S <- S_fun(X, a)
  return((-3 - sqrt(9 - 12*S))/(12*S))
}

derivative_of_l_EM <- function(k, X, a) {
  return(psigamma(k, deriv = 0) - log(k) - S_fun(X, a))
}

second_derivative_of_l <- function(k) {
  return(psigamma(k, deriv = 1) - 1/k)
}


Newtons_method_EM <- function(X, a, derivative, s_derivative, epsilon = 10^(-6)) {
  i = 0
  best_est <- estimate_of_shape_fun(X, a)
  while(abs(derivative(best_est, X, a) / s_derivative(best_est)) >= epsilon) {
    best_est <- best_est - derivative(best_est, X, a) / s_derivative(best_est)
    if(i == 50) {
      print("Brak zbieznosci!")
      return(NA)
    }
    i = i + 1
  }
  return(best_est)
}
####################
#GENEROWANIE DANYCH#
####################


gamma_mixture_pdf <- function(x, parameters) { #Gestosc mieszaniny rozkladow gamma w punkcie x
  epsilons <- parameters[,1]
  shapes <- parameters[,2]
  scales <- parameters[,3]
  
  value <- rowSums(sapply(1:nrow(parameters), function(y) dgamma(x, shape = shapes[y], scale = scales[y])*epsilons[y]))
  return(value)
}

mixture_gamma <- function(parameters, size) { #funkcja generujaca probe
  
  if(sum(parameters$epsilons) != 1) return(NA)
  
  probabilities <- runif(size)
  sums <- sapply(1:nrow(parameters), function(x) sum(parameters$epsilons[1:x]))
  X_parameters <- t(sapply(probabilities, function(x) parameters[length(sums[sums <= x]) + 1,]))
  X <- sapply(1:size, function(x) rgamma(1, shape = X_parameters[x, 2][[1]], scale = X_parameters[x, 3][[1]]))
  return(X)
}
#############
#ALGORYTM EM#
#############

EM_algorithm <- function(data, max_num_of_dist = 10, reps = 200, error = 10^(-6), newton_error = 10^(-4)) {
  
  results <- list()
  
  for(K in 1:max_num_of_dist) {
    
    epsilon <- list()
    shape <- list()
    scale <- list()
    
    ints <- quantile(data, seq(0,1, length.out=K+1))
    parts <- sapply(data, function(x) max(which(ints <= x)))
    parts[parts == K+1] <- K
    parts <- factor(parts)
    data_parts <- lapply(1:K, function(x) data[parts == x])
    
    shape[[1]] <- sapply(data_parts, function(x) Newtons_method(x, derivative_of_l, second_derivative_of_l, first_est_of_k, newton_error))
    scale[[1]] <- sapply(1:K, function(x) est_of_theta_mle(data_parts[[x]], shape[[1]][x]))
    epsilon[[1]] <- rep(1/K, K)

    D=1
    i=1
    while (D > error & i < reps) {
      denominator <- rowSums(data.frame(sapply(1:K, function(x) dgamma(data, shape = shape[[i]][x], scale = scale[[i]][x])*epsilon[[i]][x])))
      pi <- sapply(1:K, function(y) sapply(1:K, function(x) dgamma(data, shape = shape[[i]][x], scale = scale[[i]][x])*epsilon[[i]][x])[,y]/denominator)
      if(any(is.na(pi))) break
      shape[[i+1]] <- sapply(1:K, function(x) Newtons_method_EM(data, pi[,x], derivative_of_l_EM, second_derivative_of_l))
      scale[[i+1]] <- sapply(1:K, function(x) sum(data*pi[,x])/(shape[[i+1]][x]*sum(pi[,x])))
      epsilon[[i+1]] <- sapply(1:K, function(x) sum(pi[,x])/length(data))
      
      D <- (sum((epsilon[[i+1]] - epsilon[[i]])^2 + (shape[[i+1]] - shape[[i]])^2 + (scale[[i+1]] - scale[[i]])^2))/K^(7/2)
      i <- i+1
    }
    final_parameters <- matrix(c(epsilon[[i]], shape[[i]], scale[[i]]), K, 3, dimnames = list(1:K, c("epsilons", "shapes", "scales")))
    results <- append(results, list(list(epsilon = epsilon, shape = shape, scale = scale, iterations = i, final_parameters = final_parameters)))
  }
  
  return(append(results,  list(trial = data)))
}

#########
#WYKRESY#
#########

### 1 ###
shapes <- c(0.5, 4)
scales <- c(1, 2)
epsilons <- c(0.4, 0.6) #ma sie sumowac do 1
parameters <- data.frame(epsilons, shapes, scales) #ramka danych z parametrami
data <- mixture_gamma(parameters, size = 500) #proba z mieszaniny rozkladow gamma

results <- EM_algorithm(data, max_num_of_dist = 8, reps = 200)

colors <- c("Gęstość danych" = "#a9d717", "Gęstość z algorytmu EM" = "#d46611")

plot1_1 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 30), size = 1, aes(col = "Gęstość danych"), n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[1]]$final_parameters), xlim = c(0, 30), inherit.aes = FALSE, size = 1, aes(col = "Gęstość z algorytmu EM"), n = 500) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 0.35)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 1)")
plot1_2 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 30), size = 1, col =  colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[2]]$final_parameters), xlim = c(0, 30), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 0.35)) +
  ylab("Gęstość") +
  xlab("Dane (K = 2)")
plot1_3 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 30), size = 1, col = colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[3]]$final_parameters), xlim = c(0, 30), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 0.35)) +
  ylab("Gęstość") +
  xlab("Dane (K = 3)")

ggarrange(plot1_1, plot1_2, plot1_3, ncol=3, common.legend = TRUE, legend = "top") #RYSUNEK 15#

plot1_overfitting <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 30), size = 1, aes(col = "Gęstość danych"), n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[8]]$final_parameters), xlim = c(0, 30), inherit.aes = FALSE, size = 1, aes(col = "Gęstość z algorytmu EM"), n = 500) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 0.35)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 8)")

plot1_overfitting #RYSUNEK 16#

### 2 ###
shapes <- c(0.5, 4, 200, 10)
scales <- c(1, 2, 0.1, 1)
epsilons <- c(0.2, 0.2, 0.3, 0.3) #ma sie sumowac do 1
parameters <- data.frame(epsilons, shapes, scales) #ramka danych z parametrami
data <- mixture_gamma(parameters, size = 500) #proba z mieszaniny rozkladow gamma

results <- EM_algorithm(data, max_num_of_dist = 4, reps = 200)

colors <- c("Gęstość danych" = "#a9d717", "Gęstość z algorytmu EM" = "#d46611")

plot2_1 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, aes(col = "Gęstość danych"), n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[1]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, aes(col = "Gęstość z algorytmu EM"), n = 500) +
  coord_cartesian(xlim = c(0,25), ylim = c(0, 0.35)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 1)") +
  theme(text = element_text(size = 8))   
plot2_2 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, col =  colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[2]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,25), ylim = c(0, 0.35)) +
  ylab("Gęstość") +
  xlab("Dane (K = 2)") +
  theme(text = element_text(size = 8))  
plot2_3 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, col = colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[3]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,25), ylim = c(0, 0.35)) +
  ylab("Gęstość") +
  xlab("Dane (K = 3)") +
  theme(text = element_text(size = 8))  
plot2_4 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, col = colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[4]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,25), ylim = c(0, 0.35)) +
  ylab("Gęstość") +
  xlab("Dane (K = 4)") +
  theme(text = element_text(size = 8))  

ggarrange(plot2_1, plot2_2, plot2_3, plot2_4, ncol=4, common.legend = TRUE, legend = "top", font.label = list(size = 20, color = "red")) #RYSUNEK 17#

### 3 ###
shapes <- c(10, 50, 400)
scales <- c(1, 0.5, 0.2)
epsilons <- c(0.2, 0.2, 0.6) #ma sie sumowac do 1
parameters <- data.frame(epsilons, shapes, scales) #ramka danych z parametrami
data <- mixture_gamma(parameters, size = 500) #proba z mieszaniny rozkladow gamma

results <- EM_algorithm(data, max_num_of_dist = 4, reps = 200)

colors <- c("Gęstość danych" = "#a9d717", "Gęstość z algorytmu EM" = "#d46611")

plot3_1 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, aes(col = "Gęstość danych"), n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[1]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, aes(col = "Gęstość z algorytmu EM"), n = 500) +
  coord_cartesian(xlim = c(0,100), ylim = c(0, 0.07)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 1)") +
  theme(text = element_text(size = 8))   
plot3_2 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, col =  colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[2]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,100), ylim = c(0, 0.07)) +
  ylab("Gęstość") +
  xlab("Dane (K = 2)") +
  theme(text = element_text(size = 8))  
plot3_3 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, col = colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[3]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,100), ylim = c(0, 0.07)) +
  ylab("Gęstość") +
  xlab("Dane (K = 3)") +
  theme(text = element_text(size = 8))  
plot3_4 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 150), size = 1, col = colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[4]]$final_parameters), xlim = c(0, 150), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,100), ylim = c(0, 0.07)) +
  ylab("Gęstość") +
  xlab("Dane (K = 4)") +
  theme(text = element_text(size = 8))  

ggarrange(plot3_1, plot3_2, plot3_3, plot3_4, ncol=4, common.legend = TRUE, legend = "top", font.label = list(size = 20, color = "red")) #RYSUNEK 18#

### 4 ###
shapes <- c(0.5, 4)
scales <- c(1, 2)
epsilons <- c(0.4, 0.6) #ma sie sumowac do 1
parameters <- data.frame(epsilons, shapes, scales) #ramka danych z parametrami
data <- mixture_gamma(parameters, size = 100) #proba z mieszaniny rozkladow gamma

results <- EM_algorithm(data, max_num_of_dist = 8, reps = 200)

colors <- c("Gęstość danych" = "#a9d717", "Gęstość z algorytmu EM" = "#d46611")

plot4_1 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black",) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 30), size = 1, aes(col = "Gęstość danych"), n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[1]]$final_parameters), xlim = c(0, 30), inherit.aes = FALSE, size = 1, aes(col = "Gęstość z algorytmu EM"), n = 500) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 0.35)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 1)")
plot4_2 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 30), size = 1, col =  colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[2]]$final_parameters), xlim = c(0, 30), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 0.35)) +
  ylab("Gęstość") +
  xlab("Dane (K = 2)")
plot4_3 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black") +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0, 30), size = 1, col = colors[1], n = 500) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[3]]$final_parameters), xlim = c(0, 30), inherit.aes = FALSE, size = 1, col = colors[2], n = 500) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 0.35)) +
  ylab("Gęstość") +
  xlab("Dane (K = 3)")

ggarrange(plot4_1, plot4_2, plot4_3, ncol=3, common.legend = TRUE, legend = "top") #RYSUNEK 19#

### 5 ###
shapes <- c(0.2, 50, 400, 1000, 700)
scales <- c(1, 2, 2, 0.2, 0.5)
epsilons <- c(0.05, 0.15, 0.5, 0.2, 0.1) #ma sie sumowac do 1
parameters <- data.frame(epsilons, shapes, scales) #ramka danych z parametrami
data <- mixture_gamma(parameters, size = 800) #proba z mieszaniny rozkladow gamma

results <- EM_algorithm(data, max_num_of_dist = 8, reps = 200)

colors <- c("Gęstość danych" = "#a9d717", "Gęstość z algorytmu EM" = "#d46611")

plot5_1 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black", bins = 50) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0.1, 1000), size = 0.5, aes(col = "Gęstość danych"), n = 8000) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[4]]$final_parameters), xlim = c(0.1, 1000), inherit.aes = FALSE, size = 0.5, aes(col = "Gęstość z algorytmu EM"), n = 8000) +
  coord_cartesian(xlim = c(0,950), ylim = c(0, 0.01)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 4)") +
  theme(text = element_text(size = 8))   
plot5_2 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black", bins = 50) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0.1, 1000), size = 0.5, aes(col = "Gęstość danych"), n = 8000) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[5]]$final_parameters), xlim = c(0.1, 1000), inherit.aes = FALSE, size = 0.5, aes(col = "Gęstość z algorytmu EM"), n = 8000) +
  coord_cartesian(xlim = c(0,950), ylim = c(0, 0.01)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 5)") +
  theme(text = element_text(size = 8))   
plot5_3 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black", bins = 50) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0.1, 1000), size = 0.5, aes(col = "Gęstość danych"), n = 8000) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[6]]$final_parameters), xlim = c(0.1, 1000), inherit.aes = FALSE, size = 0.5, aes(col = "Gęstość z algorytmu EM"), n = 8000) +
  coord_cartesian(xlim = c(0,950), ylim = c(0, 0.01)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 6)") +
  theme(text = element_text(size = 8))   
plot5_4 <- ggplot(data.frame(data), aes(x = data)) + 
  geom_histogram(aes(y=..density..), colour = "black", bins = 50) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, parameters), xlim = c(0.1, 1000), size = 0.5, aes(col = "Gęstość danych"), n = 8000) +
  stat_function(fun = function(x) gamma_mixture_pdf(x, results[[7]]$final_parameters), xlim = c(0.1, 1000), inherit.aes = FALSE, size = 0.5, aes(col = "Gęstość z algorytmu EM"), n = 8000) +
  coord_cartesian(xlim = c(0,950), ylim = c(0, 0.01)) +
  labs(color = NULL) +
  scale_color_manual(values = colors) +
  ylab("Gęstość") +
  xlab("Dane (K = 7)") +
  theme(text = element_text(size = 8))   

ggarrange(plot5_1, plot5_2, plot5_3, plot5_4, ncol=2, nrow = 2, common.legend = TRUE, legend = "top", font.label = list(size = 20, color = "red")) #RYSUNEK 20#


library("ggplot2")

k1 <- c(0.5, 1, 2, 5, 15) #ksztalt
t1 <- c(1, 1, 1, 1, 1) #skala
parameters1 <- lapply(1:5, function(x) c(k1[x], t1[x]))
domain1 <- seq(0.001, 25, 0.001)

densities1 <- as.vector(sapply(parameters1, function(x) dgamma(domain1, shape = x[[1]], scale = x[[2]])))
data1 <- data.frame('gęstość' = densities1, x = rep(domain1, 5), 'kształt' = as.factor(rep(c(0.5, 1, 2, 5, 15), each = length(domain1))))

scale1_plot <- ggplot(data1, aes(x = x, y = gęstość, colour = kształt)) +
        geom_line(size = 1) +
        coord_cartesian(xlim = c(0, 23), ylim=c(0, 0.75)) +
        scale_colour_brewer(palette = "Set1")

scale1_plot #RYSUNEK 1#
#########

k2 <- c(2, 2, 2, 2, 2) #ksztalt
t2 <- c(0.5, 1, 2, 5, 15) #skala
parameters2 <- lapply(1:5, function(x) c(k2[x], t2[x]))
domain2 <- seq(0.001, 35, 0.001)

densities2 <- as.vector(sapply(parameters2, function(x) dgamma(domain2, shape = x[[1]], scale = x[[2]])))
data2 <- data.frame('gęstość' = densities2, x = rep(domain2, 5), 'skala' = as.factor(rep(c(0.5, 1, 2, 5, 15), each = length(domain2))))

shape1_plot <- ggplot(data2, aes(x = x, y = gęstość, colour = skala)) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(0, 33), ylim=c(0, 0.75)) +
  scale_colour_brewer(palette = "Set1")

shape1_plot #RYSUNEK 2#

#########

k3 <- c(2, 10, 20, 50, 100) #ksztalt
t3 <- c(0.5, 0.5, 0.5, 0.5, 0.5) #skala
parameters3 <- lapply(1:5, function(x) c(k3[x], t3[x]))
domain3 <- seq(0.001, 70, 0.001)

densities3 <- as.vector(sapply(parameters3, function(x) dgamma(domain3, shape = x[[1]], scale = x[[2]])))
data3 <- data.frame('gęstość' = densities3, x = rep(domain3, 5), 'kształt' = as.factor(rep(c(2, 10, 20, 50, 100), each = length(domain3))))

small_scale_plot <- ggplot(data3, aes(x = x, y = gęstość, colour = kształt)) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(0, 65), ylim=c(0, 0.75)) +
  scale_colour_brewer(palette = "Set1")

small_scale_plot #RYSUNEK 3#

#########

k4 <- rep(0.8, 5) #ksztalt
t4 <- c(1, 5, 10, 15, 30) #skala
parameters4 <- lapply(1:5, function(x) c(k4[x], t4[x]))
domain4 <- seq(0.0000001, 40, 0.001)

densities4 <- as.vector(sapply(parameters4, function(x) dgamma(domain4, shape = x[[1]], scale = x[[2]])))
data4 <- data.frame('gęstość' = densities4, x = rep(domain4, 5), 'skala' = as.factor(rep(c(1, 5, 10, 15, 30), each = length(domain4))))

small_shape_plot <- ggplot(data4, aes(x = x, y = gęstość, colour = skala)) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(0, 35), ylim=c(0, 0.5)) +
  scale_colour_brewer(palette = "Set1")

small_shape_plot #RYSUNEK 4#

# library(gamstackr)
# library(mgcv)
#
# data("iris")
#
# ex1 <- lm(Sepal.Length ~ Sepal.Width, data = iris)
# ex2 <- lm(Sepal.Length ~ Petal.Length, data = iris)
# ex3 <- lm(Sepal.Length ~ Petal.Width, data = iris)
#
# den1 <- matrix(dnorm(iris$Sepal.Length, predict(ex1), sd(ex1$residuals), log = TRUE))
# den2 <- matrix(dnorm(iris$Sepal.Length, predict(ex2), sd(ex2$residuals), log = TRUE))
# den3 <- matrix(dnorm(iris$Sepal.Length, predict(ex3), sd(ex3$residuals), log = TRUE))
#
# stack_fam = NestedStack(list(den1, den2, den3), list(id(), id(), id()))
#
# stack <- gam(list(Sepal.Length ~ Species + Sepal.Width,  ~Species), data = iris, family = stack_fam)
#
# W <- predict(stack, type = "response")
#
# barplot(t(unique(W[,1:3])), col = unique(iris$Species),
#         names.arg = unique(iris$Species), xlab = "Species", ylab = "Weight")
# legend("topright", legend = c("Expert 1", "Expert 2", "Expert 3"), fill = unique(iris$Species))

library(ggplot2)

data(trees)
#1
nrow(trees)
head(trees)

#2
?trees

#3
mod1 = lm(Girth~Height, data = trees)
summary(mod1)
ggplot(trees)+
  geom_point(aes(Height,Girth))+
  geom_abline(aes(intercept = mod1$coefficients[1],slope = mod1$coefficients[2]), color = 'red')+
  labs(title = 'Simple Linear Regression: Diameter~Height', x = 'Height (ft)', y = 'Diameter (in)')

#4
mod2 = lm(Girth~Volume, data = trees)
summary(mod2)
ggplot(trees)+
  geom_point(aes(Volume,Girth))+
  geom_abline(aes(intercept = mod2$coefficients[1],slope = mod2$coefficients[2]), color = 'purple')+
  labs(title = 'Simple Linear Regression: Diameter~Volume', x = 'Volume (cubic ft)', y = 'Diameter (in)')

#5
cov(trees$Height, trees$Girth)
cov(trees$Volume, trees$Girth)

#6
h_cor = cor(trees$Height, trees$Girth)
h_cor
h_cor^2
v_cor = cor(trees$Volume, trees$Girth)
v_cor
v_cor^2

#7
mod3 = lm(Girth~ Height+Volume, data = trees)
summary(mod3)

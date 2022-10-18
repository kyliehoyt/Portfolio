library(ggplot2)
library(gridExtra)
# 1


df = cars
nrow(df)
length(df)
head(df)
?cars

# 2
ggplot(df)+geom_point(aes(speed,dist))+labs(title = 'Original Data')

# box plot
dist_box = ggplot(df, aes(y=dist)) + geom_boxplot(outlier.color = 'red')

speed_box = ggplot(df, aes(y=speed)) + geom_boxplot(outlier.color = 'red')

grid.arrange(dist_box,speed_box,nrow = 1)

# 5
boxplot.stats(df$dist)


# 6
dist_pdf = ggplot(df) + geom_density(aes(x = dist))
speed_pdf = ggplot(df) + geom_density(aes(x = speed))
grid.arrange(dist_pdf,speed_pdf,nrow = 1)


# 7
cor(df$speed, df$dist)

# 8
linearMod = lm(dist ~ speed, data=cars)
linearMod

ggplot(df) + geom_point(aes(speed,dist)) + 
  geom_abline(slope = 3.932,intercept = -17.579) +
  labs(title = 'Linear Regression of Speed and Stopping Distance of Cars',x ='Speed (mph)',y = 'Stopping Distance (ft)')

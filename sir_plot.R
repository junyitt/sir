library(reshape2)
library(ggplot2)

setwd("C:/Users/User/git/sir/sir_cplusplus")
df <- read.csv("data_sir_deterministic.csv")
View(df)
df2 <- melt(df, id.vars = "time")

setwd("C:/Users/User/git/sir_cplusplus")
df1 <- read.csv("data_sir_rand.csv")[,1:6]
View(df1)
df2 <- melt(df, id.vars = "time")
ggplot(df2, aes(x=time, y = value, group = variable, color = variable)) + geom_point()


setwd("C:/Users/User/git/sir_cplusplus")
df <- read.csv("data_sir_stochastic.csv")
# View(df)
df2 <- melt(df, id.vars = "time")


ggplot(df2, aes(x=time, y = value, group = variable, color = variable)) + geom_point()

df <- read.csv("data_sir_immi_emi.csv")
df2 <- melt(df, id.vars = "time")
ggplot(df2, aes(x=time, y = value, group = variable, color = variable)) + geom_point()

library(tidyverse)
library(lubridate)
df <- read_csv("/Users/claytonlamb/Downloads/all.aar.clayton.csv")%>%
  mutate(time_decimal = hour(date.time) + minute(date.time) / 60 + second(date.time) / 3600)

options(scipen = 999)
##find a species with lots of data to work with
df%>%
  count(species)


ggplot(df, aes(x = time_decimal)) +
  geom_histogram() +
  labs(title = "Raw data: Mule deer") +
  facet_wrap(vars(species), scales="free_y")

##
df.select <- df%>%
  filter(species=="muledeer")

df.select$aar%>%mean(na.rm=TRUE)%>%log()

df.select%>%
  filter(lag.before<=24,
         lag.after<=24)%>%
  pull(aar)%>%
  mean(na.rm=TRUE)%>%
  log()


##plot


ggplot(df.select%>%
         filter(lag.before<=24,
                lag.after<=24), aes(x = time_decimal, y = log(aar))) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_smooth(se = FALSE, method = "loess") +
    facet_wrap(vars(class)) +
    labs(title = "Raw data: Mule deer")


ggplot(df.select, aes(x = time_decimal)) +
  geom_histogram() +
  labs(title = "Raw data: Mule deer")


ggplot(df%>%
         filter(lag.before<=48,
                lag.after<=48), aes(x = time_decimal, y = log(aar))) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = FALSE, method = "loess") +
  facet_wrap(vars(species)) +
  labs(title = "Raw data")


ggplot(df%>%filter(class=="mtn.bikers",
                   lag.before<=36,
                   lag.after<=36), aes(x = time_decimal, y = log(aar))) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = FALSE, method = "loess") +
  facet_wrap(vars(species)) +
  labs(title = "Raw data")
  
  


mod <- df.select%>%
  filter(lag.before<=24,
                   lag.after<=24)%>%
  mutate(log.aar=log(aar))%>%
  glm(log.aar~scale(time_decimal) + class, data=.)
  
summary(mod)

pred.dat <- tibble(time_decimal=seq(1,24, by=0.5), class="mtn.bikers")
pred.dat$pred <- predict(mod, newdata = pred.dat)
pred.dat$se <- predict(mod, newdata = pred.dat,se=TRUE)$se.fit



ggplot(pred.dat, aes(x = time_decimal, y = pred, ymin=pred-se, ymax=pred+se)) +
  geom_point() +
  geom_linerange()+
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "modelled response")




df.select%>%
  filter(lag.before<=48,
         lag.after<=48,
         time_decimal>=4,
         time_decimal<=10)%>%
  pull(aar)%>%
  mean(na.rm=TRUE)%>%
  log()

df.select%>%
  filter(lag.before<=48,
         lag.after<=48,
         time_decimal>=15,
         time_decimal<=20)%>%
  pull(aar)%>%
  mean(na.rm=TRUE)%>%
  log()

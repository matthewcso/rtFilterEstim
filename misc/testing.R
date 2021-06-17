source('../R/rt.R')
source('../R/filter.R')

library(zoo)
library(ggplot2)

incidence <- c()

for(n in c(1:250)){
    incidence <- append(incidence, 100*1.03**n)
}

estimated_rt <- filter_and_estimate(incidence, 
                                    n_resample=1)

df <- data.frame(incidence=incidence, rt_estimate = estimated_rt$mean, t = c(1:250), ideal = rep(log(1.03), 250))

ggplot(df) +
    geom_line(aes(x=t, y=rt_estimate, color='Estimated')) +
    geom_line(aes(x=t, y=ideal, color='Ideal'))+
    labs(x='time', y='Growth rate')

print(estimated_rt)
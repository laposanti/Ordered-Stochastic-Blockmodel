
library(ggplot2)

#lognormal
mu=1
sigma=1
lmu = log(mu**2/(sqrt(mu**2+sigma)))
lsigma = sqrt(log(1+ sigma/mu**2))
gen=rlnorm(1000000,lmu, lsigma)
plot(density(gen))
abline(v=mean(gen),col="red")
df= data.frame(mean = mean(gen),
five.perc = quantile(gen, 0.05),
ninentifive.perc= quantile(gen, .95))

dlnorm_norm_param = function(y,z.mu=1,z.sigma=1){
  lmu = log(z.mu**2/(sqrt(z.mu**2+z.sigma)))
  lsigma = sqrt(log(1+ z.sigma/z.mu**2))
  dlnorm(y,lmu, lsigma)
}



ggplot(df, aes(x = mean, y = mean, ymin = five.perc, ymax = ninentifive.perc)) +
  geom_errorbar(width = 0.05, color = "red") +
  geom_point(color = "blue", size=2) +
  labs(x = "Mean", y = "95th Percentile") +
  scale_x_continuous(limits = c(0, max(df$mean) + 1)) +
  theme_classic()

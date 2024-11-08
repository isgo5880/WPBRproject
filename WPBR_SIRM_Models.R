library(tidyverse)
library(odin)

#SIRM Model for 4 Age Groups
#0-5yrs: Seedling (s)
#5-10yrs: saPling (p)
#10-20yrs: Established (e)
#20-30yrs: Young tree (y)

#Age Groups Model
sirm_age <- odin::odin({
  
  deriv(Ss) <- -Ss*beta_s*(Is+Ip+Ie+Iy) - alpha*Ss
  deriv(Is) <- -Is*(gamma_s + alpha) + beta_s*(Is+Ip+Ie+Iy)*(Ss + Rs)
  deriv(Rs) <- -Rs*beta_s*(Is+Ip+Ie+Iy) + p_s*gamma_s - alpha*Rs
  deriv(Ms) <- (1 - p_s)*gamma_s*Is
  
  deriv(Sp) <- alpha*Ss - Sp*beta_p*(Is+Ip+Ie+Iy) - alpha*Sp
  deriv(Ip) <- -Ip*(gamma_p + alpha) + beta_p*(Is+Ip+Ie+Iy)*(Sp + Rp) + alpha*Is
  deriv(Rp) <- alpha*Rs - Rp*beta_p*(Is+Ip+Ie+Iy) + p_p*gamma_p - alpha*Rp
  deriv(Mp) <- (1 - p_p)*gamma_p*Ip
  
  deriv(Se) <- alpha*Sp - Se*beta_e*(Is+Ip+Ie+Iy) - omega*Se
  deriv(Ie) <- -Ie*(gamma_e + omega) + beta_e*(Is+Ip+Ie+Iy)*(Se + Re) + alpha*Ip
  deriv(Re) <- alpha*Rp - Re*beta_e*(Is+Ip+Ie+Iy) + p_e*gamma_e - omega*Re
  deriv(Me) <- (1 - p_e)*gamma_e*Ie
  
  deriv(Sy) <- omega*Se - Sy*beta_y*(Is+Ip+Ie+Iy)
  deriv(Iy) <- -Iy*gamma_y + Ie*omega + beta_y*(Is+Ip+Ie+Iy)*(Sy + Ry)
  deriv(Ry) <- omega*Re - Ry*beta_y*(Is+Ip+Ie+Iy) + p_y*gamma_y
  deriv(My) <- (1 - p_y)*gamma_y*Iy
    
  #Initial Conditions and Equations
  #Infected has been found to vary between 0 and 7.3%
  initial(Ss) <- (0.17)*(0.99)
  initial(Is) <- (0.17)*(0.01)
  initial(Rs) <- 0
  initial(Ms) <- 0
  
  initial(Sp) <- (0.11)*(0.99)
  initial(Ip) <- (0.11)*(0.01)
  initial(Rp) <- 0
  initial(Mp) <- 0
  
  initial(Se) <- (0.29)*(0.99)
  initial(Ie) <- (0.29)*(0.01)
  initial(Re) <- 0
  initial(Me) <- 0
  
  initial(Sy) <- (0.43)*(0.99)
  initial(Iy) <- (0.43)*(0.01)
  initial(Ry) <- 0
  initial(My) <- 0
  
  gamma_s <- user(0.5) #1/2 years^-1
  gamma_p <- user(0.2) #1/5 years^-1
  gamma_e <- user((1/7)) #1/7 years^-1
  gamma_y <- user(0.1) #1/10 years^-1
  
  #stress - pH, geography, elevation, climate change/pollution, shade/sunlight/growing conditions
  #stress <- 0.5 #stress is a level between 0(not stressed at all) and 1(practically dead with stress)
  #susceptibility <- 0.67 #susceptibility is a factor of age/height - the younger and shorter, the more susceptible
  #effectiveness depends on the presence of secondary hosts, climate(the presence of supercooled leaves), and geography
  #effectiveness <- 0.9 #a level between 0(rust is not present) and 1(incredibly effective)
  # starting level of sigma is 0.27
  #sigma <- stress*susceptibility*effectiveness #sigma is a value between 0 and 1
  sigma <- 0.3
  
  alpha <- 0.2 #1/5 years^-1
  omega <- 0.1 #1/10 years^-1
  beta_s <- 0.5 + 2*sigma
  beta_p <- 0.2 + 0.8*sigma
  beta_e <- (1/7) + (4/7)*sigma
  beta_y <- 0.1 + 0.4*sigma
  p_s <- 0.30 - 0.25*sigma
  p_p <- 0.40 - 0.30*sigma
  p_e <- 0.50 - 0.35*sigma
  p_y <- 0.60 - 0.40*sigma
})

model_age <- sirm_age$new()

t <- seq(from=0, to=50, by=0.1) 
sol_age <- as_tibble(data.frame(model_age$run(t)))

sol_age %>% 
  pivot_longer(-t) %>% 
  mutate(compartment=substr(name,1,1)) %>%
  mutate(agegrp=substr(name,2,2)) %>%
  ggplot(aes(x=t, y=value, col=factor(compartment,levels=c("S","I","R","M")), lty=factor(agegrp,levels=c("s","p","e","y")))) + 
  geom_line(linewidth=1) +
  scale_color_manual(values=c("S"="blue","I"="red","R"="green","M"="black")) + 
  theme_classic() + 
  theme(text=element_text(size=12), legend.title=element_blank()) 
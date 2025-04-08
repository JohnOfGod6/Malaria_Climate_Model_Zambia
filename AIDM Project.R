gc()
rm(list=ls())
rm()

#** load needed libraries*

library(deSolve)
library(zoo)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyr)

#** read in your climate data and select your country from the dataset*
#*
all_met_df=readRDS("C:/Users/user/Downloads/Group project AIDM/student_met_data.rds")
 
met_df=all_met_df%>%filter(country =="zambia")
View(met_df)

#** dataframe to store your scenario information *
scen_df=data.frame(country="zambia",scenario_name=c("climatology","scen_1.9","scen_4.5","scen_8.5"),
                   Temp_scen=c("climatology_Temp","dy_Temp_1.9","dy_Temp_4.5","dy_Temp_8.5"),
                   ppt_scen=c("climatology_ppt","dy_ppt_1.9","dy_ppt_4.5","dy_ppt_8.5"),
                   start_year=c(1988,1988,1988,1988),
                   end_year=c(2023,2035,2035,2035))

View(scen_df)

res_df=list() # list to store scenario results
for(v in 1:nrow(scen_df)) try({
  country=scen_df[v,]$country
  scenario_name=scen_df[v,]$scenario_name
  ppt_scen=scen_df[v,]$ppt_scen
  Temp_scen=scen_df[v,]$Temp_scen
  start_year=scen_df[v,]$start_year
  
  end_year=scen_df[v,]$end_year

ppt=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,][,ppt_scen]))  #** <- load you rainfall data, with the correct date range*
Temp=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,][,Temp_scen]))  #** <- load you Temperature data, with the correct date range*

# get the year and day of the year too
days=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,]$day))
years=unlist(unname(met_df[met_df$year>=start_year & met_df$year<=end_year,]$year))

#** create the mosquito population model equations from equation (1) in White et al  2011 here*
mpop <- function(t,y,parms) {
  
  beta <- beta[pmax(1,ceiling(t))] # daily egg laying rate
  d_E <- d_E[pmax(1,ceiling(t))] # duration of egg development
  d_L <- d_L[pmax(1,ceiling(t))] # duration of larvae development
  d_P <- d_P[pmax(1,ceiling(t))] # duration of pupae development
  mu0E <- mu0E[pmax(1,ceiling(t))] # baseline daily mortality rate of eggs
  K <- K[pmax(1,ceiling(t))] # carrying capacity
  mu0L <- mu0L[pmax(1,ceiling(t))] # baseline daily mortality rate of larvae
  muP <- muP[pmax(1,ceiling(t))] # daily mortality rate of pupae
  muM <- muM[pmax(1,ceiling(t))] # daily mortality rate of adult
  gamma <- gamma[pmax(1,ceiling(t))] # relative density dependent morality
  
  # transmission parameters
  bE<-bE[pmax(1, ceiling(t))] # mosquito contact rate
  bI<-bI[pmax(1, ceiling(t))] # mosquito EIP
  fIH<-fIH[pmax(1, ceiling(t))] # fraction of infectious humans
  
  # the states
  E=y["E"];L=y["L"];P=y["P"]
  SM=unname(y["SM"]) ;EM=unname(y["EM"]) ;IM=unname(y["IM"]);
  
#**build out your coupled model here*
  
  dE= beta * (SM+EM+IM) - (E / d_E) - mu0E*(1+(E+L)/K) * E #eggs
  dL=  (E / d_E) - (L / d_L) - mu0L*(1+gamma*(E+L)/K) * L  #larvae  
  dP=  (L / d_L) - (P / d_P) - muP * P #pupae
  dSM= P/(2*d_P) - bE * fIH * SM - muM * SM #susceptible mosquitoes
  dEM= bE * fIH * SM - bI * EM - muM * EM #exposed mosquitoes
  dIM= bI * EM - muM * IM #infected mosquitoes

    list(c(dE,dL,dP,dSM,dEM,dIM))
}
 
# set the duration of simulation in days. Use the length of rain data
numDay=length(ppt)

############################# carrying capacity  #############################
#** use the carrying capacity from equation (8) in White et al 2011* 
#** you will need tau, mean_rain and lambda to derive K the carrying capacity*

tau=7 #** <-insert the length of days over which rainfall contributes to the carrying capacity*

mean_rain=rollmean(ppt, k=tau, fill=NA, align="right") #** <-insert the rolling mean rainfall during the past tau days*

lambda =1e6 #** the site scaling factor*

K=lambda*mean_rain #** compute K the site-adjusted carrying-capacity*
K[is.na(K)]=mean(K, na.rm = T) # replace missing values in K with the average

############################## development times ##############################
#** instead of the d_E,d_L and d_P in White et al, use the temperature-regulated development here*
#** See the Table 2 exercise handout for temperature formulation for each immature stage *

# temperature-regulated egg duration
a_d_E=1.3346 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out *
b_d_E=13.1687 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_d_E=13.8354 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
d_d_E=5.6711 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *

d_E=a_d_E+(b_d_E/(1+((Temp/c_d_E)^d_d_E))) #** <- compute temperature-regulated development duration here *

# temperature-regulated egg to larave duration
a_d_L=15.405231 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_d_L=-8.461096 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_d_L=22.166679 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
d_d_L=-16.260036 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
d_L=a_d_L+(b_d_L/(1+((Temp/c_d_L)^d_d_L))) #** <- compute temperature-regulated development duration here *
  
# temperature-regulated larvae to pupae duration
a_d_P=10.740867 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_d_P=-9.549660 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_d_P=12.876982 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
d_d_P=-4.231138 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
d_P=a_d_P+(b_d_P/(1+((Temp/c_d_P)^d_d_P))) #** <- compute temperature-regulated development duration here *
  
#################################################################################

############################ daily mortality rates ############################
#** instead of the muE0,muL and muP in White et al, use the temperature and moisture regulated mortality here*
#** See the  Table 2 exercise assignment for temperature and moisture formulations for each immature stage *

############################ moisture regulated mortality rates ############################
#** But first calculate the dessication days D and the probability immature survival due to moisture following Parham et al 2012 *

#** for dessication, i.e., drying days, how many days has it been without rain (without rain defined as rainfall less than 1 mm)*
#** run the lines of code below to calculate D*
D=numeric()
D[1]=0

for(i in 2:length(ppt)){
  # get the last rain event
  rain_event=tail(which(ppt[1:(i-1)]>=1),1)
  
  if(length(rain_event)==1){
  # how many days since it last rained
  D[i]=length(which(ppt[(rain_event):i]<1))
  }else{
    D[i]=0  
  }
  }

#** Egg *
wE=0.405 #** <- insert egg sensitivity to dessication here*
p_E_R=2*exp(-wE*D)/(1+exp(-wE*D)) #** <- insert formulation relating dessication to Egg survival probability*

#** Larvae *
wL=.855 #** <- insert larvae sensitivity to dessication here*
p_L_R=2*exp(-wL*D)/(1+exp(-wL*D)) #** <- insert formulation relating dessication to larvae survival probability*

#** Pupae *
wP=0.602 #** <- insert larvae sensitivity to dessication here*
p_P_R=2*exp(-wP*D)/(1+exp(-wP*D)) #** <- insert formulation relating dessication to larvae survival probability*
  
################################################################################
################# temperature regulated survival probability ###################

#** Egg *
a_p_E_T=-0.25223 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out *
b_p_E_T=0.0936 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_E_T=-0.002026 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *

p_E_T=(a_p_E_T+(b_p_E_T*Temp)+(c_p_E_T*Temp^2)) #** <- compute temperature-regulated Egg survival probability due to temperature here *

#** Larvae *
a_p_L_T=0.3594 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_p_L_T=0.04983 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_L_T=-0.0010 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
p_L_T=(a_p_L_T+(b_p_L_T*Temp)+(c_p_L_T*Temp^2)) #** <- compute temperature-regulated pupae survival probability due to temperature here *
  
#** Pupae *

a_p_P_T=0.19067 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_p_P_T=0.02643 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_P_T=-5.3e-05 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
p_P_T=(a_p_P_T+(b_p_P_T*Temp)+(c_p_P_T*Temp^2)) #** <- compute temperature-regulated pupae survival probability due to temperature here *
  
  
#** Adult *
a_p_M_T=-0.000828 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
b_p_M_T= 0.0367#** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
c_p_M_T=0.522 #** <- insert the corresponding coefficient values from Table 2, in the exercise hand out  *
  
p_M_T=a_p_M_T* Temp^2 + b_p_M_T *Temp + c_p_M_T #** <- compute temperature-regulated adult survival probability due to temperature here *
  
################################### daily mortality rates #####################################

mu0E=-log(p_E_R*p_E_T) #** <- compute the daily mortality rate from probability of Egg survival due to temperature and dessication*
  
mu0L=-log(p_L_R*p_L_T) #** <- compute the daily mortality rate from probability of Larvae survival due to temperature and dessication*

#** set differential mortality due to density dependent according to Prior median in Table 1, White et al 2011*
gamma=rep(13.06, numDay) #** <- insert relative density dependent morality value here *

muP=-log(p_P_R*p_P_T)  #** <- compute the daily mortality rate from probability of Pupae survival due to temperature and dessication*

muM=-log(p_M_T) #** <- compute the daily mortality rate. See  Table 2 exercise handout for probability of adult survival due to temperature*


################################### Oviposition #####################################

#** To calculate the daily number of eggs laid per adult, you will need *
#** emax, GP (the length of the gonotrophy cycle), and muM (daily mortality rate)*
#** use the equation 3 in White el 2011 to derive a temperature-regulated daily eggs laid*
#** for this, use the temperature-regulated GP (see  Table 2 exercise handout) and temperature-dependent muM calculated above for adults*

emax=93.6 #** <- insert maximum number of eggs mosquitoes can lay during an oviposition. See table 1 prior median *

GP=1/pmin(1,pmax(0,0.017*Temp-0.165))  #** <-insert gonotrophy duration due to temperature. see  Table 2 exercise handout for formulation*
beta=emax*(muM/(exp(GP*muM)-1)) #** calculate beta, the daily number of eggs laid per adult using equation 3 in white et al 2011*
################################################################################


######################## Transmission parameters ###############################
bE=pmin(1,pmax(0,0.017*Temp-0.165))   # shapiro et al 2017, see supplement
bI=pmin(1,pmax(0,(0.000112*Temp*(Temp-15.384)*sqrt(35-Temp)))) # Mordecai et al 2013  and Understanding the link between malaria risk and climate


# assume a force of infection from human populaton that is constant= 0.2
fIH=rep(0.20,numDay)

# set the initial state conditions
#** use the following initial state conditions *
Mpop=1e5 # initial adult population
EM0=1e4 # exposed mosquitoes at time t = 0
IM0=0 # infected & infectious mosquitoes at time t = 0
SM0=Mpop-IM0-EM0 # susceptible mosquitoes at time t = 0
E0=0 # initial egg population
L0=P0=0 # initial larvae and pupae population

# the time sequence over which to simulate
times=seq(0,numDay, by=1)

# initialize the model
init_states=c(E=E0, L=L0,P=P0,SM=SM0,EM=EM0,IM=IM0)

# run the model. simulate mosquito population
sim_out=ode(y=init_states, times=times, func = mpop, parms = NULL,method=NULL)



#**Get the EIR number of infectious bites per human per day EIR=maZ*
m=rowSums(sim_out[,c("SM","EM","IM")])/(Mpop/5) # mosquito density per human. Here assume mosquitoes are 5 times the number of humans from t=0 
a=bE[1:numDay] # mosquito biting rate, i.e. the mosquito contact rate
Z=sim_out[,"IM"]/rowSums(sim_out[,c("SM","EM","IM")]) #Proportion of mosquitoes infectious
EIR=m[-(numDay+1)]*a*Z[-(numDay+1)] 


res_df[[v]]=data.frame(country=country,scenario_name=scenario_name,day=days, year=years,
                       E=sim_out[,"E"][-(numDay+1)],L=sim_out[,"L"][-(numDay+1)],P=sim_out[,"P"][-(numDay+1)],
                       SM=sim_out[,"SM"][-(numDay+1)],EM=sim_out[,"EM"][-(numDay+1)],IM=sim_out[,"IM"][-(numDay+1)],
                       EIR=EIR,
                       oviposition=beta,sporogony=1/bI,gonotrophy_length=1/bE,
                       egg_duration=d_E,larvae_duration=d_L,pupae_duration=d_P,
                       egg_survival=exp(-mu0E),larvae_survival=exp(-mu0L),
                       pupae_survival=exp(-muP),
                       adult_survival=exp(-muM))

},silent = T)


View(do.call(rbind, res_df))

plot(res_df[[v]])

################################################################################
#                                                                              #
#                              presentation questions                          #
#                                                                              #
################################################################################
 

#**Q1: Describe the current climate conditions of your country*
#**Use the plot of the annual trend and seasonality of temperature and rainfall to describe this.*
#**Use the 1988-2023 observed_Temp and observed_ppt as your temperature and rainfall variables from the climate dataset to do this.*

# Ensure 'day' is in Date format
met_df$day <- as.Date(met_df$day)

# Filter for relevant period and country
climate_df <- met_df %>%
  filter(country == "zambia", year >= 1988 & year <= 2023) %>%
  select(year, day, observed_Temp, observed_ppt)

# Annual Trends  
annual_climate <- climate_df %>%
  group_by(year) %>%
  summarise(mean_temp = mean(observed_Temp, na.rm = TRUE),
            total_ppt = mean(observed_ppt, na.rm = TRUE))

# Plot annual temperature trend
annual_trend_temp <- ggplot(annual_climate, aes(x = year, y = mean_temp)) +
  geom_line(color = "red") +
  geom_point(color = "red") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "darkred") +
  labs(title = "Annual Temperature Trend (Zambia: 1988–2023)",
       x = "Year", y = "Mean Temperature (°C)") +
  theme_minimal()

# Plot annual rainfall trend
annual_trend_ppt <- ggplot(annual_climate, aes(x = year, y = total_ppt)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "darkblue") +
  labs(title = "Annual Rainfall Trend (Zambia: 1988–2023)",
       x = "Year", y = "Total Rainfall (mm)") +
  theme_minimal()

#   Daily Seasonality 
# Extract day-of-year
climate_df <- climate_df %>%
  mutate(doy = as.integer(format(day, "%j")))

# Daily average over all years
seasonal_climate_daily <- climate_df %>%
  group_by(doy) %>%
  summarise(mean_temp = mean(observed_Temp, na.rm = TRUE),
            mean_ppt = mean(observed_ppt, na.rm = TRUE))

# Plot daily temperature
seasonal_temp <- ggplot(seasonal_climate_daily, aes(x = doy, y = mean_temp)) +
  geom_line(color = "red") +
  labs(title = "Daily Average Temperature (Seasonality)",
       x = "Day of Year", y = "Temperature (°C)") +
  theme_minimal()

# Plot daily rainfall
seasonal_ppt <- ggplot(seasonal_climate_daily, aes(x = doy, y = mean_ppt)) +
  geom_line(color = "blue") +
  labs(title = "Daily Average Rainfall (Seasonality)",
       x = "Day of Year", y = "Rainfall (mm)") +
  theme_minimal()

# Show all plots  
grid.arrange(annual_trend_temp, annual_trend_ppt,
             seasonal_temp, seasonal_ppt,ncol=2)

#**Q2: According to the climate data, what is the motivation for investigating how climate change may alter the risk of malaria transmission?*
#**To help you answer this, compare the recent annual trend of temperature and rainfall from 2018-2023 (taken as recent) to the climate average of temperature and rainfall (based on data from 1988-2023)*
#**Also, compare the seasonality of rainfall and temperature from 2018-2023 to the daily climatology of rainfall & temperature (based on data from 1988-2023)*
#**Discuss the differences between the annual trends & climate averages and also the differences between recent seasonality vs monthly climatology*
#**why would a change in the averages and seasonality be a concern to the risk of malaria transmission*

 
met_df$day <- as.Date(met_df$day)
met_df$month <- as.integer(format(met_df$day, "%m"))

df_climate <- met_df %>% filter(year >= 1988 & year <= 2023)
df_recent  <- met_df %>% filter(year >= 2018 & year <= 2023)

#Annual Averages
annual_recent <- df_recent %>%
  group_by(year) %>%
  summarise(
    observed_Temp = mean(observed_Temp, na.rm = TRUE),
    observed_ppt = mean(observed_ppt, na.rm = TRUE)
  )

climate_avg_temp <- mean(df_climate$observed_Temp, na.rm = TRUE)
climate_avg_ppt  <- mean(df_climate$observed_ppt, na.rm = TRUE)

# Annual Plots 
p1 <- ggplot(annual_recent, aes(x = year, y = observed_Temp)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = climate_avg_temp, linetype = "dashed", color = "red") +
  labs(title = "Annual Temperature (2018–2023)", y = "°C", x = "Year") +
  annotate("text", x = 2018.5, y = climate_avg_temp + 0.2, label = "Clim. Avg", color = "red") +
  theme_minimal()

p2 <- ggplot(annual_recent, aes(x = year, y = observed_ppt)) +
  geom_line(color = "darkgreen", size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = climate_avg_ppt, linetype = "dashed", color = "red") +
  labs(title = "Annual Rainfall (2018–2023)", y = "Rainfall (mm)", x = "Year") +
  annotate("text", x = 2018.5, y = climate_avg_ppt + 5, label = "Clim. Avg", color = "red") +
  theme_minimal()

grid.arrange(p1, p2, nrow = 1)

#Monthly Seasonality
monthly_climate <- df_climate %>%
  group_by(month) %>%
  summarise(
    temp_clim = mean(observed_Temp, na.rm = TRUE),
    rain_clim = mean(observed_ppt, na.rm = TRUE)
  )

monthly_recent <- df_recent %>%
  group_by(month) %>%
  summarise(
    temp_recent = mean(observed_Temp, na.rm = TRUE),
    rain_recent = mean(observed_ppt, na.rm = TRUE)
  )

seasonality <- left_join(monthly_climate, monthly_recent, by = "month")

#Temperature Seasonality Plot with Legend
p3 <- ggplot(seasonality, aes(x = month)) +
  geom_line(aes(y = temp_clim, color = "Climatology"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = temp_recent, color = "Recent"), size = 1.2) +
  scale_color_manual(values = c("Climatology" = "red", "Recent" = "blue")) +
  labs(title = "Monthly Temperature: Climatology vs Recent", y = "°C", x = "Month", color = "Legend") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal()

# Rainfall Seasonality Plot with Legend
p4 <- ggplot(seasonality, aes(x = month)) +
  geom_line(aes(y = rain_clim, color = "Climatology"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = rain_recent, color = "Recent"), size = 1.2) +
  scale_color_manual(values = c("Climatology" = "red", "Recent" = "darkgreen")) +
  labs(title = "Monthly Rainfall: Climatology vs Recent", y = "Rainfall (mm)", x = "Month", color = "Legend") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_minimal()

grid.arrange(p3, p4, nrow = 1)

#**Q3: Build out your coupled model of mosquito population and malaria transmission. And describe the system*

# Coupled mosquito population and malaria transmission model
mpop <- function(t, y, parms) {
  beta <- beta[pmax(1, ceiling(t))]
  d_E <- d_E[pmax(1, ceiling(t))]
  d_L <- d_L[pmax(1, ceiling(t))]
  d_P <- d_P[pmax(1, ceiling(t))]
  mu0E <- mu0E[pmax(1, ceiling(t))]
  K <- K[pmax(1, ceiling(t))]
  mu0L <- mu0L[pmax(1, ceiling(t))]
  muP <- muP[pmax(1, ceiling(t))]
  muM <- muM[pmax(1, ceiling(t))]
  gamma <- gamma[pmax(1, ceiling(t))]
  
  bE <- bE[pmax(1, ceiling(t))]
  bI <- bI[pmax(1, ceiling(t))]
  fIH <- fIH[pmax(1, ceiling(t))]
  
  E <- y["E"]; L <- y["L"]; P <- y["P"]
  SM <- y["SM"]; EM <- y["EM"]; IM <- y["IM"]
  
  dE <- beta * (SM + EM + IM) - (E / d_E) - mu0E * (1 + (E + L) / K) * E
  dL <- (E / d_E) - (L / d_L) - mu0L * (1 + gamma * (E + L) / K) * L
  dP <- (L / d_L) - (P / d_P) - muP * P
  dSM <- -bE * fIH * SM - muM * SM
  dEM <- bE * fIH * SM - bI * EM - muM * EM
  dIM <- bI * EM - muM * IM
  
  return(list(c(dE, dL, dP, dSM, dEM, dIM)))
}

# Initial state values
init_states <- c(E = 0, L = 0, P = 0, 
                 SM = 90000, EM = 10000, IM = 0)

# Time steps for simulation
times <- seq(0, numDay, by = 1)

# Run simulation
sim_out <- ode(y = init_states, times = times, func = mpop, parms = NULL)

# Calculate EIR: Entomological Inoculation Rate
m <- rowSums(sim_out[, c("SM", "EM", "IM")]) / (Mpop / 5)
a <- bE[1:numDay]
Z <- sim_out[, "IM"] / rowSums(sim_out[, c("SM", "EM", "IM")])
 


#**Q4: Describe the mosquito lifecycle and malaria transmission parameters based on temperature and moisture regulation entomology in your country*
#**Be sure to plot & discuss how your local climate conditions of temperature and rainfall dictate the oviposition rate, the sporogonic cycle, the gonotrophic cycle, the development time taken for egg, larvae and pupae to develop, and the survival of egg, larvae and pupae and adults*
  
# Prepare your climatology data
climatology_data <- res_df[[1]] %>%
  group_by(day) %>%
  summarise(
    oviposition = mean(oviposition, na.rm = TRUE),
    sporogony = mean(sporogony, na.rm = TRUE),
    gonotrophy = mean(gonotrophy_length, na.rm = TRUE),
    egg_d = mean(egg_duration, na.rm = TRUE),
    larva_d = mean(larvae_duration, na.rm = TRUE),
    pupa_d = mean(pupae_duration, na.rm = TRUE)
  )

# Oviposition
gg_oviposition <- ggplot(climatology_data, aes(x = day, y = oviposition)) +
  geom_line(color = "darkgreen", size = 1.1) +
  labs(title = "Oviposition Rate", x = "Day of Year", y = "Eggs per Day (β)") +
  theme_minimal(base_size = 14)

# Sporogonic Cycle
gg_sporogony <- ggplot(climatology_data, aes(x = day, y = sporogony)) +
  geom_line(color = "purple", size = 1.1) +
  labs(title = "Sporogonic Cycle", x = "Day of Year", y = "Days (1/bI)") +
  theme_minimal(base_size = 14)

# Gonotrophic Cycle
gg_gonotrophy <- ggplot(climatology_data, aes(x = day, y = gonotrophy)) +
  geom_line(color = "orange", size = 1.1) +
  labs(title = "Gonotrophic Cycle", x = "Day of Year", y = "Days (1/bE)") +
  theme_minimal(base_size = 14)

# Immature Development Times
q4_dev_long <- climatology_data %>%
  pivot_longer(cols = c(egg_d, larva_d, pupa_d), names_to = "stage", values_to = "duration")

gg_development <- ggplot(q4_dev_long, aes(x = day, y = duration, color = stage)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = c("egg_d" = "red", "larva_d" = "blue", "pupa_d" = "black"),
                     labels = c("Egg", "Larva", "Pupa")) +
  labs(title = "Immature Development Time", x = "Day of Year", y = "Duration (days)", color = "Stage") +
  theme_minimal(base_size = 14)

# Arrange in 2x2 grid
grid.arrange(gg_oviposition, gg_sporogony, gg_gonotrophy, gg_development, nrow = 2)



# Prepare rainfall-based parameters (climatology)
q4_surv <- res_df[[1]] %>%
  group_by(day) %>%
  summarise(
    egg_survival = mean(egg_survival, na.rm = TRUE),
    larvae_survival = mean(larvae_survival, na.rm = TRUE),
    pupae_survival = mean(pupae_survival, na.rm = TRUE),
    adult_survival = mean(adult_survival, na.rm = TRUE)
  )

#   Egg Survival vs Rainfall
gg_egg_surv <- ggplot(q4_surv, aes(x = day, y = egg_survival)) +
  geom_line(color = "red", size = 1.1) +
  labs(title = "Egg Survival Probability", x = "Day of Year", y = "Survival Probability") +
  theme_minimal(base_size = 14)

# Plot   Larvae Survival vs Rainfall
gg_larvae_surv <- ggplot(q4_surv, aes(x = day, y = larvae_survival)) +
  geom_line(color = "blue", size = 1.1) +
  labs(title = "Larvae Survival Probability", x = "Day of Year", y = "Survival Probability") +
  theme_minimal(base_size = 14)

# Plot   Pupae Survival vs Rainfall
gg_pupae_surv <- ggplot(q4_surv, aes(x = day, y = pupae_survival)) +
  geom_line(color = "black", size = 1.1) +
  labs(title = "Pupae Survival Probability", x = "Day of Year", y = "Survival Probability") +
  theme_minimal(base_size = 14)

# Plot   Adult Survival vs Temperature
gg_adult_surv <- ggplot(q4_surv, aes(x = day, y = adult_survival)) +
  geom_line(color = "gray30", size = 1.1) +
  labs(title = "Adult Survival Probability", x = "Day of Year", y = "Survival Probability") +
  theme_minimal(base_size = 14)

# Arrange in 2x2 grid
grid.arrange(gg_egg_surv, gg_larvae_surv, gg_pupae_surv, gg_adult_surv, nrow = 2)



# use model simulations to answer the questions below
#**Q5:According to the daily climatological conditions (calculated based on 1988-2023 data), what is the expected annual EIR?*
#**To do answer this, run the model using the climatology_Temp and climatology_ppt as inputs for your Temp and ppt variables.*
#**Now Let the model run from 1988 to 2023, then calculate the annual EIR of the last 10 years (2014 to 2023) and then take the average of the annual EIRs over these 10 years.*
#**Use this as the expected annual EIR according to the climatology*
#*

# Extract simulation results for climatology scenario
climatology_data2 <- res_df[[1]]  # Already simulated from 1988–2023 using climatology_Temp and climatology_ppt

#   Filter data for last 10 years (2014 to 2023)
eir_last10yrs <- climatology_data2 %>%
  filter(year >= 2014 & year <= 2023)

#   Compute annual EIR  
annual_eir <- eir_last10yrs %>%
  group_by(year) %>%
  summarise(total_EIR = sum(EIR, na.rm = TRUE))

#   Compute the mean expected annual EIR across the 10 years
expected_annual_EIR <- mean(annual_eir$total_EIR, na.rm = TRUE)

# Output results
print(annual_eir)
cat("\nExpected Annual EIR (2014–2023 average):", round(expected_annual_EIR, 2), "\n")




#**Q6:Plot & examine the average daily EIR (i.e., daily seasonality) calculated from these last 10 years.*
#**What is the average and the range of number of infectious bites per day (the EIR) one would expect to receive in a season?*
#**Also according to the climatological EIR, when is transmission most/least likely to happen*
 

#  Create a "date" variable from the day of year (for nice month labels)
eir_last10yrs <- res_df[[1]] %>%
  filter(year >= 2014 & year <= 2023) %>%
  mutate(date = as.Date(day - 1, origin = "2023-01-01"))  # dummy year

# Group by day and calculate mean EIR across years
avg_daily_eir <- eir_last10yrs %>%
  group_by(date) %>%
  summarise(mean_EIR = mean(EIR, na.rm = TRUE))

#  Plot daily EIR with month labels
ggplot(avg_daily_eir, aes(x = date, y = mean_EIR)) +
  geom_line(color = "darkgreen", size = 1.2) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(
    title = "Average Daily EIR (2014–2023) under Climatology",
    subtitle = "Monthly trend of infectious bites per person",
    x = "Month", y = "Mean EIR (bites/person/day)"
  ) +
  theme_minimal(base_size = 14)


# Calculate summary statistics
mean_eir <- mean(avg_daily_eir$mean_EIR, na.rm = TRUE)
range_eir <- range(avg_daily_eir$mean_EIR, na.rm = TRUE)

peak_day <- avg_daily_eir$date[which.max(avg_daily_eir$mean_EIR)]
low_day  <- avg_daily_eir$date[which.min(avg_daily_eir$mean_EIR)]
 
cat(" Average daily EIR:", round(mean_eir, 3), "bites/person/day\n")
cat(" Range: ", round(range_eir[1], 3), "to", round(range_eir[2], 3), "\n")
cat(" Peak transmission around:", format(peak_day, "%B %d"), "\n")
cat(" Lowest transmission around:", format(low_day, "%B %d"), "\n")




# use future climate simulations to answer the questions below
#**Q7:Now run the model using the future climate scenario dy_Temp_1.9(for temperature), and dy_ppt_1.9(for rainfall). Let the model run from 1988 to 2035*
#**Compute and plot the average annual EIR over the years 2030-2035, compare it to the expected annual EIR from Q5.Discuss any differences in risk*
#**Plot & examine the average daily EIR(i.e., daily seasonality) under this future climate (2030-2035). Plot & discuss the seasonality in the EIR compared to the climatology EIR from Q6.*
#**What is the average and the range of infectious bites one would expect to receive in a season?*
#**Also under this future scenario, when is transmission most/least likely to happen*
 

future_1.9_data <- res_df[[which(scen_df$scenario_name == "scen_1.9")]]

# Filter for years 2030–2035
eir_future_1.9 <- future_1.9_data %>%
  filter(year >= 2030 & year <= 2035)

# Annual EIR totals
annual_eir_1.9 <- eir_future_1.9 %>%
  group_by(year) %>%
  summarise(total_EIR = sum(EIR, na.rm = TRUE))

# Average annual EIR for 2030–2035
expected_eir_1.9 <- mean(annual_eir_1.9$total_EIR, na.rm = TRUE)

# Create baseline for climatology
baseline_df <- data.frame(
  year = 2030:2035,
  total_EIR = expected_annual_EIR,
  Scenario = "Climatology Avg"
)

# Add scenario label
annual_eir_1.9$Scenario <- "RCP 1.9"

# Combine both
combined_eir <- bind_rows(annual_eir_1.9, baseline_df)

# Plot annual comparison
ggplot(combined_eir, aes(x = year, y = total_EIR, color = Scenario, linetype = Scenario)) +
  geom_line(size = 1.2) +
  geom_point(data = filter(combined_eir, Scenario == "RCP 1.9"), size = 2) +
  scale_color_manual(values = c("RCP 1.9" = "darkgreen", "Climatology Avg" = "blue")) +
  scale_linetype_manual(values = c("RCP 1.9" = "solid", "Climatology Avg" = "dashed")) +
  labs(
    title = "Annual EIR (2030–2035) under RCP 1.9",
    subtitle = paste("Avg EIR (RCP 1.9):", round(expected_eir_1.9, 1),
                     "| Climatology Avg:", round(expected_annual_EIR, 1)),
    x = "Year", y = "Total Annual EIR",
    color = "Scenario", linetype = "Scenario"
  ) +
  theme_minimal(base_size = 14)


# RCP 1.9 Daily EIR (2030–2035)
daily_eir_1.9 <- eir_future_1.9 %>%
  group_by(day) %>%
  summarise(mean_EIR = mean(EIR, na.rm = TRUE)) %>%
  mutate(date = as.Date(day - 1, origin = "2023-01-01"),  # dummy non-leap year
         Scenario = "RCP 1.9")

# Climatology Daily EIR (2014–2023)
avg_daily_eir <- res_df[[1]] %>%
  filter(year >= 2014 & year <= 2023) %>%
  group_by(day) %>%
  summarise(mean_EIR = mean(EIR, na.rm = TRUE)) %>%
  mutate(date = as.Date(day - 1, origin = "2023-01-01"),
         Scenario = "Climatology")

# Combine for comparison
combined_seasonality <- bind_rows(daily_eir_1.9, avg_daily_eir)

ggplot(combined_seasonality, aes(x = date, y = mean_EIR, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Climatology" = "blue",
                                "RCP 1.9" = "darkgreen")) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(
    title = "Comparison of Daily EIR Seasonality",
    subtitle = "Climatology(2014–2023) vs RCP 1.9(2030–2035)",
    x = "Month", y = "Mean EIR (bites/person/day)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)


# Stats for RCP 1.9
mean_daily_eir_1.9 <- mean(daily_eir_1.9$mean_EIR, na.rm = TRUE)
range_daily_eir_1.9 <- range(daily_eir_1.9$mean_EIR, na.rm = TRUE)
peak_day_1.9 <- daily_eir_1.9$date[which.max(daily_eir_1.9$mean_EIR)]
low_day_1.9  <- daily_eir_1.9$date[which.min(daily_eir_1.9$mean_EIR)]


cat("Avg daily EIR:", round(mean_daily_eir_1.9, 3), "bites/person/day\n")
cat("Range:", round(range_daily_eir_1.9[1], 3), "to", round(range_daily_eir_1.9[2], 3), "\n")
cat("Peak transmission:", format(peak_day_1.9, "%B %d"), "\n")
cat("Lowest transmission:", format(low_day_1.9, "%B %d"), "\n")
cat("Avg annual EIR (RCP 1.9):", round(expected_eir_1.9, 1), "\n")
cat("Climatology Avg Annual EIR (Q5):", round(expected_annual_EIR, 1), "\n")


#**Q8:Run the model now, but using the future climate scenario dy_Temp_4.5, and dy_ppt_4.5, let the model run from 1988 to 2035*
#**Compute and plot the average annual EIR over the years 2030-2035, compare it to the expected annual EIR from Q5.Discuss any differences in risk*
#**Plot & examine the average daily EIR(i.e., daily seasonality) under this future climate (2030-2035). Plot & discuss the seasonality in the EIR compared to the climatology EIR from Q6.*
#**What is the average and the range of infectious bites one would expect to receive in a season?*
#**Also under this future scenario, when is transmission most/least likely to happen*
 
# Get RCP 4.5 model result from res_df
future_4.5_data <- res_df[[which(scen_df$scenario_name == "scen_4.5")]]

# Filter years 2030–2035
eir_future_4.5 <- future_4.5_data %>%
  filter(year >= 2030 & year <= 2035)

# Annual EIR
annual_eir_4.5 <- eir_future_4.5 %>%
  group_by(year) %>%
  summarise(total_EIR = sum(EIR, na.rm = TRUE))

expected_eir_4.5 <- mean(annual_eir_4.5$total_EIR)

# Add labels
annual_eir_4.5$Scenario <- "RCP 4.5"
baseline_df <- data.frame(
  year = 2030:2035,
  total_EIR = expected_annual_EIR,
  Scenario = "Climatology Avg"
)

combined_eir <- bind_rows(annual_eir_4.5, baseline_df)

# Plot
ggplot(combined_eir, aes(x = year, y = total_EIR, color = Scenario, linetype = Scenario)) +
  geom_line(size = 1.2) +
  geom_point(data = filter(combined_eir, Scenario == "RCP 4.5"), size = 2) +
  scale_color_manual(values = c("RCP 4.5" = "darkorange", "Climatology Avg" = "blue")) +
  scale_linetype_manual(values = c("RCP 4.5" = "solid", "Climatology Avg" = "dashed")) +
  labs(
    title = "Annual EIR (2030–2035) under RCP 4.5",
    subtitle = paste("Avg EIR (RCP 4.5):", round(expected_eir_4.5, 1),
                     "| Climatology Avg:", round(expected_annual_EIR, 1)),
    x = "Year", y = "Total Annual EIR",
    color = "Scenario", linetype = "Scenario"
  ) +
  theme_minimal(base_size = 14)

# RCP 4.5 daily seasonality (2030–2035)
daily_eir_4.5 <- eir_future_4.5 %>%
  group_by(day) %>%
  summarise(mean_EIR = mean(EIR, na.rm = TRUE)) %>%
  mutate(date = as.Date(day - 1, origin = "2023-01-01"),
         Scenario = "RCP 4.5")

# Climatology from Q6 (2014–2023)
avg_daily_eir <- res_df[[1]] %>%
  filter(year >= 2014 & year <= 2023) %>%
  group_by(day) %>%
  summarise(mean_EIR = mean(EIR, na.rm = TRUE)) %>%
  mutate(date = as.Date(day - 1, origin = "2023-01-01"),
         Scenario = "Climatology")

# Combine both
combined_seasonality <- bind_rows(daily_eir_4.5, avg_daily_eir)



ggplot(combined_seasonality, aes(x = date, y = mean_EIR, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Climatology" = "blue",
                                "RCP 4.5" = "darkorange")) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(
    title = "Comparison of Daily EIR Seasonality",
    subtitle = "Climatology(2014–2023) vs RCP 4.5 (2030–2035)",
    x = "Month", y = "Mean EIR (bites/person/day)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)



mean_daily_eir_4.5 <- mean(daily_eir_4.5$mean_EIR, na.rm = TRUE)
range_daily_eir_4.5 <- range(daily_eir_4.5$mean_EIR, na.rm = TRUE)
peak_day_4.5 <- daily_eir_4.5$date[which.max(daily_eir_4.5$mean_EIR)]
low_day_4.5  <- daily_eir_4.5$date[which.min(daily_eir_4.5$mean_EIR)]


cat("Avg daily EIR:", round(mean_daily_eir_4.5, 3), "bites/person/day\n")
cat("Range:", round(range_daily_eir_4.5[1], 3), "to", round(range_daily_eir_4.5[2], 3), "\n")
cat("Peak transmission:", format(peak_day_4.5, "%B %d"), "\n")
cat("Lowest transmission:", format(low_day_4.5, "%B %d"), "\n")
cat("Avg annual EIR (RCP 4.5):", round(expected_eir_4.5, 1), "\n")
cat("Climatology Avg Annual EIR (Q5):", round(expected_annual_EIR, 1), "\n")


#**Q9:Again, run the model but using the future climate scenario dy_Temp_8.5, and dy_ppt_8.5. Let the model run from 1988 to 2035*
#**Compute & plot the average annual EIR over the years 2030-2035, compare it to the expected annual EIR from Q5.Discuss any differences in risk*
#**Plot & examine the average daily EIR(i.e., daily seasonality) under this future climate (2030-2035).Plot & discuss the seasonality in the EIR compared to the climatology EIR from Q6.*
#**What is the average and the range of infectious bites one would expect to receive in a season?*
#**Also under this future scenario, when is transmission most/least likely to happen*


# Extract model output for RCP 8.5 scenario
future_8.5_data <- res_df[[which(scen_df$scenario_name == "scen_8.5")]]

# Filter for target years
eir_future_8.5 <- future_8.5_data %>%
  filter(year >= 2030 & year <= 2035)


# Summarise annual EIR
annual_eir_8.5 <- eir_future_8.5 %>%
  group_by(year) %>%
  summarise(total_EIR = sum(EIR, na.rm = TRUE))

# Mean across 2030–2035
expected_eir_8.5 <- mean(annual_eir_8.5$total_EIR)

# Label scenarios
annual_eir_8.5$Scenario <- "RCP 8.5"
baseline_df <- data.frame(
  year = 2030:2035,
  total_EIR = expected_annual_EIR,
  Scenario = "Climatology Avg"
)

# Combine for plotting
combined_eir_8.5 <- bind_rows(annual_eir_8.5, baseline_df)

# Plot annual EIR
ggplot(combined_eir_8.5, aes(x = year, y = total_EIR, color = Scenario, linetype = Scenario)) +
  geom_line(size = 1.2) +
  geom_point(data = filter(combined_eir_8.5, Scenario == "RCP 8.5"), size = 2) +
  scale_color_manual(values = c("RCP 8.5" = "firebrick", "Climatology Avg" = "blue")) +
  scale_linetype_manual(values = c("RCP 8.5" = "solid", "Climatology Avg" = "dashed")) +
  labs(
    title = "Annual EIR (2030–2035) under RCP 8.5",
    subtitle = paste("Avg EIR (RCP 8.5):", round(expected_eir_8.5, 1),
                     "| Climatology Avg:", round(expected_annual_EIR, 1)),
    x = "Year", y = "Total Annual EIR",
    color = "Scenario", linetype = "Scenario"
  ) +
  theme_minimal(base_size = 14)


# RCP 8.5: Daily seasonality (2030–2035)
daily_eir_8.5 <- eir_future_8.5 %>%
  group_by(day) %>%
  summarise(mean_EIR = mean(EIR, na.rm = TRUE)) %>%
  mutate(date = as.Date(day - 1, origin = "2023-01-01"),
         Scenario = "RCP 8.5")

# Climatology: Daily EIR (already computed in Q6)
avg_daily_eir <- res_df[[1]] %>%
  filter(year >= 2014 & year <= 2023) %>%
  group_by(day) %>%
  summarise(mean_EIR = mean(EIR, na.rm = TRUE)) %>%
  mutate(date = as.Date(day - 1, origin = "2023-01-01"),
         Scenario = "Climatology")

# Combine both for comparison
combined_seasonality_8.5 <- bind_rows(daily_eir_8.5, avg_daily_eir)


ggplot(combined_seasonality_8.5, aes(x = date, y = mean_EIR, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Climatology" = "blue",
                                "RCP 8.5" = "firebrick")) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(
    title = "Comparison of Daily EIR Seasonality",
    subtitle = "Climatology (2014–2023) vs RCP 8.5(2030–2035)",
    x = "Month", y = "Mean EIR (bites/person/day)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 14)


mean_daily_eir_8.5 <- mean(daily_eir_8.5$mean_EIR, na.rm = TRUE)
range_daily_eir_8.5 <- range(daily_eir_8.5$mean_EIR, na.rm = TRUE)
peak_day_8.5 <- daily_eir_8.5$date[which.max(daily_eir_8.5$mean_EIR)]
low_day_8.5  <- daily_eir_8.5$date[which.min(daily_eir_8.5$mean_EIR)]
 
cat("Avg daily EIR:", round(mean_daily_eir_8.5, 3), "bites/person/day\n")
cat("Range:", round(range_daily_eir_8.5[1], 3), "to", round(range_daily_eir_8.5[2], 3), "\n")
cat("Peak transmission:", format(peak_day_8.5, "%B %d"), "\n")
cat("Lowest transmission:", format(low_day_8.5, "%B %d"), "\n")
cat("Avg annual EIR (RCP 8.5):", round(expected_eir_8.5, 1), "\n")
cat("Climatology Avg Annual EIR (Q5):", round(expected_annual_EIR, 1), "\n")

# Conclusions
#**Q10: Reflect over the future risk predicted under the different the three scenarios. Are the risk of transmission, as you would expect?.*
#**How similar are they across future conditions? Which scenario is more concerning for malaria transmission?*


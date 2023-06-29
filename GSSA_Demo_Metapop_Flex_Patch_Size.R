#########################################
# Demographics SIRS metapopulation model#
#########################################
rm(list = ls())

library(GillespieSSA)
library(tidyverse)

# Define Paramenters
patchPopSize <-     c(500, 200, 100, 100)    # Patch size
initial_infected <- c(  1,   0,   0,   0)    # Initial infected
U <- length(patchPopSize)                    # Number of patches
simName <- "SIRS metapopulation model"       # Simulation name
tf <- 500                                    # Final time

#Collect parameters
parms <- list(
  sigma = 1/20,                          # E to I rate
  gamma = 0.1,                           # I to R rate
  omega = 0.005,                         # R to S rate
  mu = 1/250,                            # Birth/death rate per person per day
  alpha = 1/100) 

#Transmission terms
beta = 0.8
within_pop_contact = 1
between_pop_contact = 0.005/U     # normalised by number of patches 

nextgen_matrix <- matrix(nrow = U, ncol = U, data = 0)

for(i in 1:U){
  for(j in 1:U){
    parms[[paste0("beta_",i,i)]] = within_pop_contact*beta
    nextgen_matrix[i,i] = within_pop_contact*beta*(1/parms$gamma)
    parms[[paste0("beta_",j,i)]] = between_pop_contact*beta
    nextgen_matrix[j,i] = between_pop_contact*beta*(1/parms$gamma)
    nextgen_matrix[i,j] = between_pop_contact*beta*(1/parms$gamma)
    parms[[paste0("beta_",j,j)]] = within_pop_contact*beta
    nextgen_matrix[j,j] = within_pop_contact*beta*(1/parms$gamma)
  }
  parms[[paste0("N", i)]] = patchPopSize[i]
}

eigenvalues <- eigen(nextgen_matrix, only.values = T)

R0 <- max(abs(eigenvalues$values)) 

EEI <- function(R0, I)



#Create the named initial state vector for the U-patch system.

x0 <- unlist(lapply(
  seq_len(U), 
  function(i){ 
    c(patchPopSize[i] - initial_infected[i], initial_infected[i], 0, 0, patchPopSize[i])
  }
))

names(x0) <- unlist(lapply(seq_len(U), function(i) paste0(c("S","E","I", "R", "N"), i)))


# Define the state change matrix for a single patch
nu <- matrix(c( -1,  0,  0, +1, +1, -1,  0,  0,  0,  0, # S
                +1, -1,  0,  0,  0,  0, -1,  0,  0,  0, # E
                 0, +1, -1,  0,  0,  0,  0, -1,  0, -1, # I
                 0,  0, +1, -1,  0,  0,  0,  0, -1,  0, # R 
                 0,  0,  0,  0, +1, -1 ,-1, -1, -1, -1), # N
             nrow=5,byrow=TRUE)

# Define propensity functions
# Mass-action
a <-
  unlist(lapply(
    seq_len(U),
    function(patch) {
      i <- patch
      patches <- 1:U
      #j <- if (patch == 1) U else patch - 1
      other_patches <- patches[-i]
      patch_beta <- c()
      for(k in (1:(U-1))){
        patch_beta[k] = paste0("+(beta_", other_patches[k],i, "*I", other_patches[k], "/N", other_patches[k], ")*S", i)
      }
      c(
        paste0("(beta_", i, i, "*I", i,"/N", i, ")*S",i, paste0(patch_beta, collapse="")), # Infection
        paste0("sigma*E", i),                                       # Becomes infecious
        paste0("gamma*I", i),                                       # Recovery from infection
        paste0("omega*R", i),       # Loss of immunity
        paste0("mu*N", i),                             # Births
        paste0("mu*S", i),                                             # Deaths (S)
        paste0("mu*E", i),                                             # Deaths (E)
        paste0("mu*I", i),                                             # Deaths (I)
        paste0("mu*R", i),                                             # Deaths (R)
        paste0("alpha*I", i)                                           # Deaths from infection
        
      )
    }
  ))

# Run simulations with the Direct method
 set.seed(2)
out <- ssa(
  x0 = x0,
  a = a,
  nu = nu,
  parms = parms,
  tf = tf,
  method = ssa.d(),
  simName = simName,
  verbose = FALSE,
  consoleInterval = 1
) 


# Built-in plotting function
# ssa.plot(out, by = 2, show.title = TRUE, show.legend = T) #all states

## Extra Plots
plot_data <- out$data %>%
  as_tibble() %>%
  pivot_longer(!t, names_to = "ID", values_to = "count") %>%
  separate(ID, 
           into = c("state", "patch"), 
           sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(state = factor(state, levels = c("S", "E", "I", "R", "N"))) %>%
  filter(state != "N")

ggplot(data = plot_data, aes(x=t, y=count, colour=state))+
  geom_line()+
  facet_wrap(~factor(patch, levels = unique(patch)) ,ncol = 1, scales = "free_y")+
  labs(x="Time",
       y="Frequency")+
  theme_bw()

## Table showing extinction/transmission info for each patch

extinct_data <- out$data %>%
  as_tibble() %>%
  slice_max(t) %>%
  distinct() %>%
  pivot_longer(!t, names_to = "ID", values_to = "count") %>%
  separate(ID, 
           into = c("state", "patch"), 
           sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(state = factor(state, levels = c("S", "E", "I", "R", "N")),
         persist = case_when(state=="I" & count > 0 ~ T, 
                             state=="I" & count == 0 ~ F)) %>%
  drop_na() %>%
  select(patch, count, persist)
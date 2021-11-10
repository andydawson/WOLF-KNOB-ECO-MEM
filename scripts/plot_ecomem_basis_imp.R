library(ggplot2)
library(rstan)
library(reshape2)
library(mgcv)

#######################################################################################
## spherical decay function
#######################################################################################

spherical <- function(x, phi) {
  
  # if (x <= 0) {
  #   reject("x negative spherical function undefined; found x = ", x)
  # }
  if (x > phi) { g=0 } else {
    g = (1 - 3.0/2.0 * x/phi + 1.0/2.0 * (x/phi)^3) 
  }
  
  return(g)
}


#######################################################################################
## load and extract posterior samples
#######################################################################################
# fit = readRDS('scripts-stan/output/fit_ecomem_basis_imp.RDS')
# 
# suff = 'imp'

fit = readRDS(paste0(path_output, '/fit_ecomem_basis_imp_', suffix, '.RDS'))

plot(fit)

names(fit)
post = extract(fit, inc_warmup=TRUE)

names(post)


#######################################################################################
## load data
#######################################################################################

dat = readRDS(paste0(path_output, '/data_ecomem_basis_imp_', suffix, '.RDS'))

n_basis = dat$n_basis
N_sites = dat$N_sites
lag = dat$lag

#######################################################################################
## plot output
#######################################################################################

if (!dir.exists(path_figures)){
  dir.create(path_figures)
}

dir.create(paste0(path_figures, '/', suffix))
path_figures = paste0(path_figures, '/', suffix)
  
## continuous memory
w_num = as.character(seq(0,n_basis))
w_num[which(nchar(w_num)==1)] = paste0('0', w_num[which(nchar(w_num)==1)])

w = post$w
w_quants = apply(w, 2, quantile, c(0.05, 0.5, 0.95))
w_dat = data.frame(par=paste0('w', w_num), t(w_quants))
colnames(w_dat) = c('par', 'q5', 'q50', 'q95')

ggplot(data=w_dat) +
  geom_point(aes(x=par, y=q50)) +
  geom_linerange(aes(x=par, ymin=q5, ymax=q95)) +
  theme_bw() #+
  #scale_y_discrete(limits=rev)

vals = seq(0, lag)
w_dat$vals = vals

ggplot(data=w_dat) + 
  geom_line(aes(x=vals, y=q50),) +
  geom_ribbon(aes(x=vals, ymin=q5, ymax=q95), fill='dodgerblue', alpha=0.5) +
  theme_bw() +
  theme(text = element_text(size=16)) +
  ylab("Continuous \n Antecedent Weight") +
  scale_x_continuous(name="Lag", breaks=seq(0, lag))
ggsave(paste0(path_figures, '/cmem_antecedent-weight-_', suffix, '.png'))

# fit GAM model
m <- gam(q50 ~ s(vals, k=4), data = w_dat, method='REML') 
mlo <- gam(q5 ~ s(vals, k=4), data = w_dat, method='REML') 
mhi <- gam(q95 ~ s(vals, k=4), data = w_dat, method='REML') 

# define finer grid of predictor values
lag_fine <- seq(0, lag, by = 0.1) 

# apply predict() function to the fitted GAM model
# using the finer grid of x values
p <- predict(m, newdata = data.frame(vals = lag_fine), se = TRUE) 
plo <- predict(mlo, newdata = data.frame(vals = lag_fine), se = TRUE) 
phi <- predict(mhi, newdata = data.frame(vals = lag_fine), se = TRUE) 


# plot the estimated mean values of y (fit) at given x values
# over the finer grid of x values;
# superimpose approximate 95% confidence band for the true 
# mean values of y at given x values in the finer grid
w_dat_fine <- data.frame(x = lag_fine, 
                fit = p$fit,
                lwr = plo$fit, 
                upr = phi$fit)


ggplot(data=w_dat_fine) + 
  geom_line(aes(x=x, y=fit),) +
  geom_ribbon(aes(x=x, ymin=lwr, ymax=upr), fill='dodgerblue', alpha=0.5) +
  theme_bw() +
  theme(text = element_text(size=16)) +
  ylab("Continuous \n Antecedent Weight") +
  scale_x_continuous(name="Lag", breaks=seq(0, lag))
ggsave(paste0(path_figures, '/cmem_antecedent-weight-smooth_', suffix, '.png'))

## gamma
gamma_idx = which(substr(names(post), 1, 5) == "gamma")

gamma_dim = rep(NA, length(gamma_idx))
gamma_names = c()
for (i in 1:length(gamma_idx)){
  gamma_dim[i] = ncol(data.frame(post[gamma_idx[i]][[1]] ))
  gamma_name = names(post[gamma_idx[i]])
  if (gamma_dim[i] == 1) {
    gamma_name_rep = gamma_name
  } else {
    gamma_name_rep = paste0(gamma_name, '_', seq(1, gamma_dim[i]))
    #gamma_name_rep = rep(gamma_name, gamma_dim[i])
    
  }
  gamma_names = c(gamma_names, gamma_name_rep)
}

gamma = matrix(unlist(post[gamma_idx]), ncol = sum(gamma_dim), byrow = FALSE)
colnames(gamma) = gamma_names

gamma_melt = melt(gamma)
colnames(gamma_melt) = c('iter', 'variable', 'value')
# gamma_melt$type = substr(gamma_melt$Var2, 1, 6)
# colnames(gamma_melt) = c('iter', 'variable', 'value', 'type')

ggplot(data=gamma_melt) + 
  geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
  geom_vline(xintercept=0, lty=2) +
  facet_wrap(~variable, scales='free_x') +
  xlab('gamma') +
  ylab('frequency') +
  theme_bw() #+ 
  #coord_flip()
ggsave(paste0(path_figures, '/gamma_hist_', suffix, '.png'))

ggplot(data=gamma_melt) + 
  geom_line(aes(x = iter, y = value)) + 
  facet_grid(variable~., scales = "free_y") +
  xlab('iter') +
  ylab('value') +
  theme_bw()
ggsave(paste0(path_figures, '/gamma_trace_', suffix, '.png'))


# #quantile(gamma1, c(0.10, 0.5, 0.90))
# gamma.quants = apply(gamma, 2, quantile, c(0.05, 0.5, 0.95))
# gamma.quants
# # coef.dat = data.frame(par=paste0('gamma', seq(0,ncol(gamma)-1)), t(gamma.quants))
# coef.dat = data.frame(par=colnames(gamma.quants), t(gamma.quants))



## beta

beta_idx = which(substr(names(post), 1, 4) == "beta")
beta_dim = ncol(data.frame(post[beta_idx][[1]]))

beta_names = paste0('beta', seq(0, beta_dim-1))

beta = matrix(unlist(post[beta_idx]), ncol = beta_dim, byrow = FALSE)
colnames(beta) = beta_names

beta_melt = melt(beta)
#beta_melt$type = substr(beta_melt$Var2, 1, 6)
colnames(beta_melt) = c('iter', 'variable', 'value')#, 'type')

ggplot(data=beta_melt) + 
  geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
  geom_vline(xintercept=0, lty=2) +
  facet_wrap(~variable, scales='free_x') +
  xlab('beta') +
  ylab('frequency') +
  theme_bw() #+ 
#coord_flip()
ggsave(paste0(path_figures, '/beta_hist_', suffix, '.png'))

ggplot(data=beta_melt) + 
  geom_line(aes(x = iter, y = value)) + 
  facet_grid(variable~., scales = "free_y") +
  xlab('iter') +
  ylab('value') +
  theme_bw()
ggsave(paste0(path_figures, '/beta_trace_', suffix, '.png'))



## both gamma and beta
both_melt = rbind(gamma_melt, beta_melt)

level_order = c('gamma0', 'gamma1', 'beta0', 'beta1')

ggplot(data=both_melt,aes(x=value,y=factor(variable, level=rev(level_order)))) + 
  stat_summary(fun.data=function(x) {
    res <- quantile(x,probs=c(0.025,0.5,0.975))
    names(res)<-c("ymin","y","ymax")
    res}) + theme_bw() + 
  geom_vline(xintercept=0, colour="red", linetype="dashed") +
  xlab("value") + 
  ylab("variable")
ggsave(paste0(path_figures, '/gamma_beta_linerange_', suffix, '.png'))

## sigma

sigma_idx = which(substr(names(post), 1, 5) == "sigma")

sigma_dim = rep(NA, length(sigma_idx))
sigma_names = c()
for (i in 1:length(sigma_idx)){
  sigma_dim[i] = ncol(data.frame(post[sigma_idx[i]][[1]] ))
  sigma_name = names(post[sigma_idx[i]])
  if (sigma_dim[i] == 1) {
    sigma_name_rep = sigma_name
  } else {
    sigma_name_rep = paste0(sigma_name, '_', seq(1, sigma_dim[i]))
    #gamma_name_rep = rep(gamma_name, gamma_dim[i])
    
  }
  sigma_names = c(sigma_names, sigma_name_rep)
}

sigma = matrix(unlist(post[sigma_idx]), ncol = sum(sigma_dim), byrow = FALSE)
colnames(sigma) = sigma_names

sigma_melt = melt(sigma)
#beta_melt$type = substr(beta_melt$Var2, 1, 6)
colnames(sigma_melt) = c('iter', 'variable', 'value')#, 'type')

ggplot(data=sigma_melt) + 
  geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
  #geom_vline(xintercept=0, lty=2) +
  facet_grid(variable~.) + #, scales='free_x') +
  xlab('sigma') +
  ylab('frequency') +
  theme_bw() #+ 
#coord_flip()
ggsave(paste0(path_figures, '/sigma_hist_', suffix, '.png'))

ggplot(data=sigma_melt) + 
  geom_line(aes(x = iter, y = value)) + 
  facet_grid(variable~., scales = "free_y") +
  xlab('iter') +
  ylab('value') +
  theme_bw()
ggsave(paste0(path_figures, '/sigma_trace_', suffix, '.png'))

## u individual effects

u_idx = which(substr(names(post), 1, 1) == "u")

u_dim = rep(NA, length(u_idx))
u_names = c()
for (i in 1:length(u_idx)){
  u_dim[i] = ncol(data.frame(post[u_idx[i]][[1]] ))
  u_name = names(post[u_idx[i]])
  if (u_dim[i] == 1) {
    u_name_rep = u_name
  } else {
    u_name_rep = paste0(u_name, '_', seq(1, u_dim[i]))
    #gamma_name_rep = rep(gamma_name, gamma_dim[i])
    
  }
  u_names = c(u_names, u_name_rep)
}

u = matrix(unlist(post[u_idx]), ncol = sum(u_dim), byrow = FALSE)
colnames(u) = u_names

u_means = colMeans(u)
hist(u_means)

u_melt = melt(u)
colnames(u_melt) = c('iter', 'variable', 'value')

N_site_u = length(u_means)/3

# u_melt$site = NA
# u_melt$site[1:N_site_u] = 1
# u_melt$site[(N_site_u+1):(2*N_site_u)] = 2
# u_melt$site[(2*N_site_u+1):(3*N_site_u)] = 3
# 
# ggplot(data=u_melt) + 
#   geom_histogram(aes(x=value, y = (..count..)/sum(..count..))) + 
#   #geom_vline(xintercept=0, lty=2) +
#   facet_grid(site~.) + #, scales='free_x') +
#   xlab('u') +
#   ylab('frequency') +
#   theme_bw() #+ 
# #coord_flip()
# ggsave(paste0(path_figures, '/sigma_hist_', suffix, '.png'))

u_mean_dat = data.frame(u_means = u_means, site = NA)
u_mean_dat$site[1:N_site_u] = 1
u_mean_dat$site[(N_site_u+1):(2*N_site_u)] = 2
u_mean_dat$site[(2*N_site_u+1):(3*N_site_u)] = 3

ggplot(data=u_mean_dat) +
  geom_histogram(aes(x=u_means, y = (..count..)/sum(..count..))) +
  #geom_vline(xintercept=0, lty=2) +
  facet_grid(site~.) + #, scales='free_x') +
  xlab('u mean') +
  ylab('frequency') +
  theme_bw() #+
#coord_flip()
ggsave(paste0(path_figures, '/u_hist_', suffix, '.png'))

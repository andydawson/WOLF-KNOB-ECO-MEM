
# dataset info
mem_var          = 'swe.yel'
covars           = c("stream.yel","tmax.lsum")
lag              = 6
sigma            = 0.001

# stan
N_warmup         = 500
N_iter           = 1000
model_name       = 'scripts/ecomem_basis_imp_logy_0dmem.stan'
include_inits    = 0
init_file        = 'data/inits/x.RDS'

# output
serial           = '01'
suffix           = paste0(format(Sys.time(), "%Y-%m-%d-%H%M"), "-", serial)
path_output      = 'output'
path_figures     = 'output'

# go!
source("scripts/fit_ecomem_basis_imp_ndmem.R")
source("scripts/plot_ecomem_basis_imp.R")

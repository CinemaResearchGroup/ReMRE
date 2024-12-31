##########################
### Cluster Path Setup ###
##########################

INPUT_ADDRESS <- "/home/yxy1423/DATA"
OUTPUT_ADDRESS <- INPUT_ADDRESS

test_no <- 131

rep_main_max <- 6e4
burn_in_main <- rep_main_max / 2
thinning_main <- 1

### events <- c("MI", "HF", "SK")
events <- c("MI")
### cohort <- c("CKD", "OLD", "FULL")
cohort <- "CKD"
### initial <- "random" or "given"
initial <- "given"

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
			"mvtnorm", 
			### mvtnorm::dmvnorm() and mvtnorm::rmvnorm(), compute the density and generate random vector of multivariate normal distribution
			"copula",
			### copula::cCopula(), generate multivariate random variables from a copula dependence
			"foreach",
			### foreach::%dopar%, implement loop for parallel computing
			"doParallel",
			### doParallel::registerDoParallel(), register clusters for parallel computing
			"survival", 
			### survival::coxph(), fits a Cox proportional hazards regression model
			"frailtypack",
			### frailtypack::frailtyPenal(), fit a shared gamma or log-normal frailty model using a semiparametric Penalized Likelihood estimation or parametric estimation on the hazard function
			"bayestestR"
			### bayestestR::ci(), compute confidence/credible/compatibility intervals (CI) or support intervals (SI) for Bayesian and frequentist models.
)

package.check <- lapply(X = REQUIRED_PACKAGES,
					  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    require(x, character.only = TRUE)
  }
)

#############################
### Loading the functions ###
#############################

GENERATE_SIMULATION_DATASET <- function(x_mat, omega_mat, J_max, term_c, beta00_vec, beta0k_mat, beta0_vec, betak_mat, phi_vec, logtheta, copula_gene, hzd_bs_dist){

	### INPUT
	### x_mat: n by p data matrix 
	### omega_mat: n by K
	### J_max: the common maximum index for all recurrent events
	### term_c: 1 by 1 parameter controls the termination rate
	### beta00_vec: dim by 1
	### beta0k_mat: dim by K
	### beta0_vec: p by 1
	### betak_mat: p by K
	### phi_vec: K by 1
	### logtheta: 1 by 1
	### copula_gene: 1 by 1, Clayton, Frank, Gumbel
	### hzd_bs_dist: 1 by 1, the selected distribution
	
	n <- nrow(omega_mat)
	K <- ncol(omega_mat)
	
	### S0(t) = exp(- H(t)), S1(t) = exp(- a * H(t))
	### then, invS1(t) = invS0(t ^ (1 / a))
	if(hzd_bs_dist == "Exp"){
		### note: 1 lambda = exp(beta)
		invS0 <- function(x, beta00_vec, covar_vec) - log(x) / (exp(beta00_vec[1]) * covar_vec)
		invSk_list <- list()
		for(k in 1 : K){
			invSk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) - log(x) / (exp(beta0k_mat[1, kp]) * covar_vec)
		}
	} else if(hzd_bs_dist == "Weibull"){
		### note: 1 nu, 2 lambda = exp(beta)
		invS0 <- function(x, beta00_vec, covar_vec) (- log(x) / covar_vec) ^ (1 / exp(beta00_vec[1])) / exp(beta00_vec[2])
		invSk_list <- list()
		for(k in 1 : K){
			invSk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) (- log(x) / covar_vec) ^ (1 / exp(beta0k_mat[1, kp])) / exp(beta0k_mat[2, kp])
		}
	} else if(hzd_bs_dist == "Gompertz"){
		### note: 1 a in R, 2 b = exp(beta)
		invS0 <- function(x, beta00_vec) (1 / beta00_vec[1]) * log(- (beta00_vec[1] / exp(beta00_vec[2])) * log(x) / covar_vec + 1)
		invSk_list <- list()
		for(k in 1 : K){
			invSk_list[[k]] <- function(x, beta0k_mat, kp) (1 / beta0k_mat[1, kp]) * log(- (beta0k_mat[1, kp] / exp(beta0k_mat[2, kp])) * log(x) / covar_vec + 1)
		}
	} else if(hzd_bs_dist == "Log_Logistic"){
		### note: 1 kappa, 2 eta = exp(beta)
		invS0 <- function(x, beta00_vec, covar_vec) exp(beta00_vec[2]) * (1 / (x ^ (1 / covar_vec)) - 1) ^ (1 / exp(beta00_vec[1]))
		invSk_list <- list()
		for(k in 1 : K){
			invSk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) exp(beta0k_mat[2, kp]) * (1 / (x ^ (1 / covar_vec)) - 1) ^ (1 / exp(beta0k_mat[1, kp]))
		}
	}
	
	### Step 1: generate terminal time D_vec (n by 1) and recurrent gap time R_list (n by J_max by K)
	### set.seed(2024)
	u_vec <- runif(n)
	### n by 1
	v_list <- list()
	### n by J_max by K
	R_list <- list()
	### n by J_max by K
	
	### S0(t) = exp(- H(t)), S1(t) = exp(- a * H(t))
	### then, invS1(t) = invS0(t ^ (1 / a))
	covar_vec <- exp(drop(x_mat %*% beta0_vec) + drop(omega_mat %*% phi_vec))
	### n by 1
	D_vec <- invS0(x = u_vec, beta00_vec = beta00_vec, covar_vec = covar_vec)
	### n by 1	
	
	if(copula_gene == "Clayton"){
		cop_object <- claytonCopula(exp(logtheta), dim = (K + 1))
	} else if (copula_gene == "Frank"){
		cop_object <- frankCopula(exp(logtheta), dim = (K + 1))
	} else if (copula_gene == "Gumbel"){
		cop_object <- gumbelCopula(exp(logtheta), dim = (K + 1))
	}
	
	for(i in 1 : n){
		### set.seed(2024)
		v_mat_temp <- cCopula(cbind(u_vec[i], matrix(runif(J_max * K), nrow = J_max, ncol = K)), copula = cop_object, inverse = TRUE)[, - 1, drop = FALSE]
		### add "drop = FALSE" to remain the matrix structure
		### J_max by K
		v_list[[i]] <- v_mat_temp
		### J_max by K
		R_mat_temp <- matrix(0, nrow = J_max, ncol = K)
		### J_max by K
		for(k in 1 : K){
			### S0(t) = exp(- H(t)), S1(t) = exp(- a * H(t))
			### then, invS1(t) = invS0(t ^ (1 / a))
			covar_temp <- exp(drop(t(betak_mat[, k]) %*% x_mat[i, ]) + omega_mat[i, k])
			### 1 by 1
			R_mat_temp[, k] <- invSk_list[[k]](x = v_mat_temp[, k], beta0k_mat = beta0k_mat, kp = k, covar_vec = covar_temp)
			### J_max by 1
		}
		R_list[[i]] <- R_mat_temp
	}
		
	### Step 2: generate censoring time C_vec (n by 1)
	
	C_vec <- runif(n = n, min = 0, max = term_c)
	### n by 1
	
	### Step 3: generate follow-up time Y_vec (n by 1) and its indicator Delta_vec (n by 1)
	
	Y_vec <- pmin(D_vec, C_vec)
	### n by 1
	Delta_vec <- ifelse(D_vec <= C_vec, 1, 0)
	### n by 1
	
	### Step 4: generate observed gap time Y_list (n by J_max by K) and its indicator Delta_list (n by J_max by K)
	
	T_list <- C_list <- Y_list <- Delta_list <- list()
	### (n by (J_max + 1) by K) for T_list
	### (n by J_max by K)
	m_mat <- matrix(0, nrow = n, ncol = K)
	### n by K
	
	for(i in 1 : n){
		
		T_mat_temp <- matrix(0, nrow = J_max + 1, ncol = K)
		### (J_max + 1) by K starting from j = 0
		C_mat_temp <- Y_mat_temp <- Delta_mat_temp <- matrix(0, nrow = J_max, ncol = K)
		### J_max by K
		
		for(j in 1 : J_max){
						
			T_mat_temp[j + 1, ] <- apply(R_list[[i]][1 : j, , drop = FALSE], 2, sum)
			### K by 1
			### j = 1, T_mat_temp[j + 1, ] = second one
			### add "drop = F" to keep matrix structure
			C_mat_temp[j, ] <- pmax(rep(0, K), Y_vec[i] - T_mat_temp[j , ])
			### K by 1
			### j = 1, T_mat_temp[j, ] = first one = all 0s
			Y_mat_temp[j, ] <- pmin(R_list[[i]][j, ], C_mat_temp[j, ])
			### K by 1
			Delta_mat_temp[j, ] <- ifelse(R_list[[i]][j, ] <= C_mat_temp[j, ], 1, 0)
			### K by 1
			
			for(k in 1 : K){
				### add this loop for correcting m_mat counts
				if(Y_vec[i] < T_mat_temp[2, k]){
					m_mat[i, k] <- m_mat[i, k]
					### ni.k = 0
				} else {
					m_mat[i, k] <- m_mat[i, k] + ifelse(Y_vec[i] > T_mat_temp[j + 1, k], 1, 0)
					### j = 1, T_mat_temp[j + 1, ] = second one
					### elements of m_mat as large as J_max
				}
			}
		}
		
		T_list[[i]] <- T_mat_temp
		### (J_max + 1) by K
		C_list[[i]] <- C_mat_temp
		### J_max by K
		Y_list[[i]] <- Y_mat_temp
		### J_max by K
		Delta_list[[i]] <- Delta_mat_temp
		### J_max by K
	}
	M_vec <- apply(m_mat, 1, max) + 1
	### n by 1
	
	### Step 5: replace J_max with M_vec for observed gap time Y_list (n by J_max by K) and its indicator Delta_list (n by J_max by K)
	
	Y_list_short <- Delta_list_short <- list()
	### n by M_vec by K
	
	for(i in 1 : n){
	
		Y_list_short[[i]] <- Y_list[[i]][1 : M_vec[i], , drop = FALSE]
		Delta_list_short[[i]] <- Delta_list[[i]][1 : M_vec[i], , drop = FALSE]
		### add "drop = FALSE" to remain the matrix structure
	}
	
	return(list(
			Y_vec = Y_vec, 
		Delta_vec = Delta_vec, 
		   Y_list = Y_list_short, 
	   Delta_list = Delta_list_short, 
			M_vec = M_vec,
			m_mat = m_mat
	))
}

CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS <- function(u_mat, deriv_index_mat, theta){
	
	### INPUT
	### u_mat: J by M
	### 	   J for the number of obs
	###        M for the dimension of copulas
	### deriv_index_mat: J by M binary
	### theta: 1 by 1 
	### Note: double check the density of Copula function and the
	### partial derivative of Copula function evaluated at 0 or 1
	### see Page 225 of Umberto (2004) 
	
	M <- ncol(u_mat) 
	### 1 by 1, in our case: K + 1
	L_vec <- apply(deriv_index_mat, 1, sum)
	### J by 1, in our case: k for each obs
	
	u_mat_temp <- u_mat ^ (- theta)
	### produces Inf when a number is close to 0
	### produces 1 when a number is close to 1
	u_mat_temp[sapply(u_mat_temp, is.infinite)] <- .Machine$double.xmax * 1e-5
	### u_mat_temp[u_mat_temp == 1] <- .Machine$double.xmax * 1e-5
	
	log_1 <- L_vec * log(theta) + lgamma(1 / theta + L_vec) - lgamma(1 / theta)
	### J by 1
	### log_2 <- (- theta - 1) * apply(u_mat * deriv_index_mat, 1, function(x) return(sum(log(x)[!is.infinite(log(x))])))
	### J by 1
	log_2 <- (- theta - 1) * apply(u_mat_temp ^ (- 1 / theta) * deriv_index_mat, 1, function(x) return(sum(log(x)[!is.infinite(log(x))])))
	### when deriv_index_mat = 1, we replace the corresponding 0 in
	### u_mat to a sufficiently small number
	log_3 <- (- 1 / theta - L_vec) * log(apply(u_mat_temp, 1, sum) - (M - 1))
	### J by 1
	### NOTE: when all delta = 0, by definition, C(0,..,0) = C
	### however, in this formula, plugin delta = 0, log_1 = 0
	### log_2 = 0, and log_3 = log(C), consistently
	
	return(log_1 + log_2 + log_3)
}

CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL <- function(d_vec, x_mat, beta00_vec, beta0_vec, phi_vec, omega_mat, hzd_bs_dist){

	### INPUT
	### d_vec: n by 1 for terminal time
	### x_mat: n by p covariate matrix
	### beta00_vec: dim by 1 baseline hazard parameters
	### beta0_vec: p by 1 regression vector
	### phi_vec: K by 1 frailty regression vector
	### omega_mat: n by K frailty matrix
	### hzd_bs_dist: 1 by 1, the selected distribution
	
	### h0(t) = h00(t), h1(t) = h00(t) * a
	### then, h1(t) = h0(t) * a
	### logS0(t) = - H(t), logS1(t) = - a * H(t)
	### then, logS1(t) = logS0 * a
	
	if(hzd_bs_dist == "Exp"){
		### note: 1 lambda = exp(beta)
		h0 <- function(x, beta00_vec, covar_vec) exp(beta00_vec[1]) * covar_vec + x - x
		logS0 <- function(x, beta00_vec, covar_vec) - exp(beta00_vec[1]) * x * covar_vec
	} else if(hzd_bs_dist == "Weibull"){
		### note: 1 nu, 2 lambda = exp(beta)
		h0 <- function(x, beta00_vec, covar_vec) (exp(beta00_vec[2]) * exp(beta00_vec[1])) * (exp(beta00_vec[2]) * x) ^ exp(beta00_vec[1]) * covar_vec
		logS0 <- function(x, beta00_vec, covar_vec) - (exp(beta00_vec[2]) * x) ^ exp(beta00_vec[1]) * covar_vec
	} else if(hzd_bs_dist == "Gompertz"){
		### note: 1 a in R, 2 b = exp(beta)
		h0 <- function(x, beta00_vec, covar_vec) exp(beta00_vec[2]) * exp(beta00_vec[1] * x) * covar_vec
		logS0 <- function(x, beta00_vec, covar_vec) - (exp(beta00_vec[2]) / beta00_vec[1]) * (exp(beta00_vec[1] * x) - 1) * covar_vec
	} else if(hzd_bs_dist == "Log_Logistic"){
		### note: 1 kappa, 2 eta = exp(beta)
		h0 <- function(x, beta00_vec, covar_vec) (exp(beta00_vec[1]) / exp(beta00_vec[2])) * (x / exp(beta00_vec[2])) ^ (exp(beta00_vec[1]) - 1) / (1 + (x / exp(beta00_vec[2])) ^ exp(beta00_vec[1])) * covar_vec
		logS0 <- function(x, beta00_vec, covar_vec) - log(1 + (x / exp(beta00_vec[2])) ^ exp(beta00_vec[1])) * covar_vec		
	}
	
	covar_vec <- exp(drop(x_mat %*% beta0_vec + omega_mat %*% phi_vec))
	### n by 1

	h0_vec <- h0(x = d_vec, beta00_vec = beta00_vec, covar_vec = covar_vec)
	### n by 1
	h0_vec <- pmax(h0_vec, .Machine$double.xmin * 1e5)
	### h0_vec must be > 0 <-> log(h0_vec) != NaN
	
	logS0_vec <- logS0(x = d_vec, beta00_vec = beta00_vec, covar_vec = covar_vec)
	### n by 1
	S0_vec <- pmax(exp(logS0_vec), .Machine$double.xmin * 1e5)
	### S0_vec must be > 0 
	logS0_vec <- log(S0_vec)
	
	f0_vec <- h0_vec * S0_vec
	f0_vec <- pmax(f0_vec, .Machine$double.xmin * 1e5)
	### f0_vec must be > 0 
	
	return(list(
	 logf0_vec = log(f0_vec), 
		S0_vec = S0_vec, 
	 logS0_vec = logS0_vec))
	### S0 = 0 -> Inf in likelihood
}

CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT <- function(r_list, M_vec, x_mat, beta0k_mat, betak_mat, omega_mat, hzd_bs_dist){

	### INPUT
	### r_list: n by M_vec by K for recurrent time
	### M_vec: n by 1
	### x_mat: n by p covariate matrix
	### beta0k_mat: dim by K baseline hazard parameters
	### betak_mat: p by K regression vector
	### omega_mat: n by K frailty matrix
	### hzd_bs_dist: 1 by 1, the selected distribution
	
	### hk(t) = h0k(t), h1(t) = h0k(t) * a
	### then, h1(t) = hk(t) * a
	### logSk(t) = - H(t), logS1(t) = - a * H(t)
	### then, logS1(t) = logSk * a

	n <- nrow(omega_mat)
	K <- ncol(omega_mat)

	if(hzd_bs_dist == "Exp"){
		### note: 1 lambda = exp(beta)
		hk_list <- logSk_list <- list()
		for(k in 1 : K){
			hk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) exp(beta0k_mat[1, kp]) * covar_vec + x - x
			logSk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) - exp(beta0k_mat[1, kp]) * x * covar_vec
		}
	} else if(hzd_bs_dist == "Weibull"){
		### note: 1 nu, 2 lambda = exp(beta)
		hk_list <- logSk_list <- list()
		for(k in 1 : K){
			hk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) - (exp(beta0k_mat[2, kp]) * x) ^ exp(beta0k_mat[1, kp]) * covar_vec
			logSk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) - (exp(beta0k_mat[2, kp]) * x) ^ exp(beta0k_mat[1, kp]) * covar_vec
		}
	} else if(hzd_bs_dist == "Gompertz"){
		### note: 1 a in R, 2 b = exp(beta)
		hk_list <- logSk_list <- list()
		for(k in 1 : K){
			hk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) exp(beta0k_mat[2, kp]) * exp(beta0k_mat[1, kp] * x) * covar_vec
			logSk_list[[k]] <- function(x, beta0k_mat, kp) - (exp(beta0k_mat[2, kp]) / beta0k_mat[1, kp]) * (exp(beta0k_mat[1, kp] * x) - 1) * covar_vec
		}
	} else if(hzd_bs_dist == "Log_Logistic"){
		### note: 1 kappa, 2 eta = exp(beta)
		hk_list <- logSk_list <- list()
		for(k in 1 : K){
			hk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) (exp(beta0k_mat[1, kp]) / exp(beta0k_mat[2, kp])) * (x / exp(beta0k_mat[2, kp])) ^ (exp(beta0k_mat[1, kp]) - 1) / (1 + (x / exp(beta0k_mat[2, kp])) ^ exp(beta0k_mat[1, kp])) * covar_vec
			logSk_list[[k]] <- function(x, beta0k_mat, kp, covar_vec) - log(1 + (x / exp(beta0k_mat[2, kp])) ^ exp(beta0k_mat[1, kp])) * covar_vec			
		}
	}

	const_mat <- exp(x_mat %*% betak_mat + omega_mat)
	### n by K
	logfk_list <- Sk_list <- list()
	### n by M_vec by K
	
	for(i in 1 : n){
		logfk_mat_temp <- matrix(0, nrow = M_vec[i], ncol = K)
		### M_vec by K
		Sk_mat_temp <- matrix(0, nrow = M_vec[i], ncol = K)
		### M_vec by K
		for(k in 1 : K){
			logSk_vec_temp <- logSk_list[[k]](x = r_list[[i]][, k], beta0k_mat = beta0k_mat, kp = k, covar_vec = const_mat[i, k])
			### M_vec by 1
			hk_vec_temp <- hk_list[[k]](x = r_list[[i]][, k], beta0k_mat = beta0k_mat, kp = k, covar_vec = const_mat[i, k])
			### M_vec by 1
			Sk_mat_temp[, k] <- exp(logSk_vec_temp)
			### M_vec by 1
			logfk_mat_temp[, k] <- log(hk_vec_temp) + logSk_vec_temp
			### M_vec by 1
		}
		fk_mat_temp <- exp(logfk_mat_temp)
		### M_vec by K
		fk_mat_temp <- pmax(fk_mat_temp, .Machine$double.xmin * 1e5)
		### M_vec by K
		logfk_list[[i]] <- log(fk_mat_temp)
		### M_vec by K
		Sk_list[[i]] <- pmax(Sk_mat_temp, .Machine$double.xmin * 1e5)
		### M_vec by K
	}
	
	
	
	return(list(
			logfk_list = logfk_list, 
			   Sk_list = Sk_list))
}

CALCULATE_LOG_PRODUCT_COPULAS <- function(S0_vec, Sk_list, Delta_list, theta, CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS){
		
	### INPUT
	### S0_vec: n by 1 for terminal time
	### Sk_list: n by M_vec by K for recurrent time
	### Delta_list: n by M_vec by K for recurrent time
	### theta: 1 by 1
	
	subscript_index1_list <- do.call(mapply, c("cbind", list(1, Delta_list)))
	### n by M_vec by (K + 1) 
	subscript_index0_list <- do.call(mapply, c("cbind", list(0, Delta_list)))
	### n by M_vec by (K + 1)
	u_list <- do.call(mapply, c("cbind", list(S0_vec, Sk_list)))
	### n by M_vec by (K + 1)
	
	log_prod1_list <- mapply(
		     FUN = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS, 
	       u_mat = u_list, 
 deriv_index_mat = subscript_index1_list, 
           theta = theta)
	### n by M_vec by 1
	
	log_prod0_list <- mapply(
	         FUN = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS, 
		   u_mat = u_list, 
 deriv_index_mat = subscript_index0_list, 
           theta = theta)
	### n by M_vec by 1
		
	return(list(log_prod0 = unlist(lapply(log_prod0_list, sum)), 
				log_prod1 = unlist(lapply(log_prod1_list, sum))))
	### n by 1
}

VECTORIZED_LOG_LIKELIHOOD_FUNCTION <- function(x_mat, M_vec, Y_vec, Delta_vec, Y_list, Delta_list, omega_mat, theta, beta00_vec, beta0_vec, beta0k_mat, betak_mat, phi_vec, hzd_bs_dist, CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL, CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT, CALCULATE_LOG_PRODUCT_COPULAS, CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS){

	### INPUT
	### x_mat: n by p covariate vector matrix
	### M_vec: n by 1 
	### Y_vec: n by 1 for terminal time
	### Delta_vec: n by 1 for terminal time
	### Y_list: n by Mi by K for recurrent time
	### Delta_list: n by Mi by K for recurrent time
	### omega_mat: n by K frailty matrix
	### theta: 1 by 1
	### beta00_vec: dim by 1
	### beta0_vec: p by 1
	### beta0k_mat: dim by K
	### betak_mat: p by K
	### phi_vec: K by 1
	### hzd_bs_dist: 1 by 1, the selected distribution
		
	### Step 1 compute the functions for the terminal event

	res_0 <- CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL(
				  d_vec = Y_vec, 
				  x_mat = x_mat,
			 beta00_vec = beta00_vec,				  
			  beta0_vec = beta0_vec, 
				phi_vec = phi_vec, 
			  omega_mat = omega_mat, 
		    hzd_bs_dist = hzd_bs_dist)
	logf0_vec <- res_0$logf0_vec
	### n by 1
	logS0_vec <- res_0$logS0_vec
	### n by 1
	S0_vec <- res_0$S0_vec
	### n by 1
	### use for copula computation
		
	### Step 2 compute the functions for the recurrent events
	
	res_k <- CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT(
				 r_list = Y_list, 
				  M_vec = M_vec, 
				  x_mat = x_mat, 
			 beta0k_mat = beta0k_mat,
			  betak_mat = betak_mat, 
			  omega_mat = omega_mat, 
		    hzd_bs_dist = hzd_bs_dist)
	logfk_list <- res_k$logfk_list
	### n by M_vec by K
	Sk_list <- res_k$Sk_list
	### n by M_vec by K
		
	### Step 3 compute the copulas
	
	res_c <- CALCULATE_LOG_PRODUCT_COPULAS(
		   S0_vec = S0_vec, 
          Sk_list = Sk_list, 
	   Delta_list = Delta_list, 
	        theta = theta, 
	CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS)
	log_prod0 <- res_c$log_prod0
	### n by 1
	log_prod1 <- res_c$log_prod1
	### n by 1
		
	log_1_vec <- unlist(lapply(Map("*", Delta_list, logfk_list), sum))
	### n by 1
	### Map("*", Delta_list, logfk_list) -> n by M_vec by K
	
	log_2_vec <- Delta_vec * (logf0_vec + log_prod1)
	### n by 1
	
	log_3_vec <- (1 - Delta_vec) * (- (M_vec - 1) * logS0_vec + log_prod0)
	
	###return(sum(log_1_vec + log_2_vec + log_3_vec))
	return(log_1_vec + log_2_vec + log_3_vec)
}

SAMPLING_MAIN_BAYES <- function(x_mat, Y_vec, Delta_vec, Y_list, Delta_list, M_vec, omega_mat_ini, Theta_vec_ini, rep_main_max, size_istep_mat, size_pstep_vec, index_cond_1, index_cond_2, hzd_bs_dist, CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL, CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT, CALCULATE_LOG_PRODUCT_COPULAS, CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS, VECTORIZED_LOG_LIKELIHOOD_FUNCTION){

	### INPUT
	### x_mat: n by p
	### Y_vec: n by 1
	### Delta_vec: n by 1
	### Y_list: n by J_max by K
	### Delta_list: n by J_max by K
	### M_vec: n by 1
	### omega_mat_ini: n by K
	### Theta_vec_ini: G by 1
	### rep_main_max: 1 by 1
	### size_istep_mat: n by K
	### size_pstep_vec: G by 1
	### index_cond_1: dim by 1
	### index_cond_2: dim by 1
	### hzd_bs_dist: 1 by 1, the selected distribution 
	
	n <- dim(omega_mat_ini)[1]
	K <- dim(omega_mat_ini)[2]
	p <- dim(x_mat)[2]	
	G <- length(Theta_vec_ini)
	eps <- 1e-3
	
	if(hzd_bs_dist == "Exp"){
		dim_hzd <- 1
	} else if (hzd_bs_dist == "Weibull"){
		dim_hzd <- 2
	} else if (hzd_bs_dist == "Gompertz"){
		dim_hzd <- 2
	} else if (hzd_bs_dist == "Log-Logistic"){
		dim_hzd <- 2
	}
	
	omegaikq_mat_updt <- omega_mat_ini
	### n by K, initial values
	### do not save the sample of omega_mat
	accept_rate_mat <- matrix(0, nrow = n, ncol = K)
	### n by K, for omega_mat
	
	Thetagq_vec_updt <- Theta_vec_ini
	### G by 1, initial values
	Theta_est_mat <- matrix(0, nrow = rep_main_max, ncol = G)
	### rep_main_max by G
	accept_rate_vec <- rep(0, G)
	### G by 1, for Theta_vec
	
	pb <- txtProgressBar(min = 1, max = rep_main_max, style = 3)
	rep_main <- 1
		
	while(rep_main <= rep_main_max){

		setTxtProgressBar(pb, rep_main)
		### show the progress bar

		sigma2_vec <- Thetagq_vec_updt[(dim_hzd + dim_hzd * K + p + p * K + K + 1 + 1) : (dim_hzd + dim_hzd * K + p + p * K + K + 1 + K)]
		### K by 1
		if(K == 1){
			Sigma_mat <- sqrt(sigma2_vec) %*% t(sqrt(sigma2_vec))
		} else {
			rho <- Thetagq_vec_updt[dim_hzd + dim_hzd * K + p + p * K + K + 1 + K + 1]
			### 1 by 1
			Sigma_mat <- sqrt(sigma2_vec) %*% t(sqrt(sigma2_vec)) * rho
			### K by K
			diag(Sigma_mat) <- sigma2_vec
		}

		for(k in 1 : K){
		
			omegaikq_mat_temp <- omegaikq_mat_prev <- omegaikq_mat_updt
			### n by K
		
			omegaikq_mat_temp[, k] <- omegaikq_mat_prev[, k] + size_istep_mat[, k] * rnorm(n)
			### n by 1
			### only k-th columns are different
			
			logARikq_vec <- VECTORIZED_LOG_LIKELIHOOD_FUNCTION(
				x_mat = x_mat, 
				M_vec = M_vec, 
				Y_vec = Y_vec, 
			Delta_vec = Delta_vec, 
			   Y_list = Y_list, 
		   Delta_list = Delta_list, 
			omega_mat = omegaikq_mat_temp,
		   beta00_vec = Thetagq_vec_updt[1 : dim_hzd],
		   beta0k_mat = matrix(Thetagq_vec_updt[(dim_hzd + 1) : (dim_hzd + dim_hzd * K)], dim_hzd, K, byrow = FALSE),
			beta0_vec = Thetagq_vec_updt[(dim_hzd + dim_hzd * K + 1) : (dim_hzd + dim_hzd * K + p)],
			betak_mat = matrix(Thetagq_vec_updt[(dim_hzd + dim_hzd * K + p + 1) : (dim_hzd + dim_hzd * K + p + p * K)], p, K, byrow = FALSE),
			  phi_vec = Thetagq_vec_updt[(dim_hzd + dim_hzd * K + p + p * K + 1) : (dim_hzd + dim_hzd * K + p + p * K + K)], 
				theta = exp(Thetagq_vec_updt[dim_hzd + dim_hzd * K + p + p * K + K + 1]),
		  hzd_bs_dist = hzd_bs_dist, 
							CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL = CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL, 
						   CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT = CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT, 
									 CALCULATE_LOG_PRODUCT_COPULAS = CALCULATE_LOG_PRODUCT_COPULAS, 
			CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS) - 
			VECTORIZED_LOG_LIKELIHOOD_FUNCTION(
				x_mat = x_mat, 
				M_vec = M_vec, 
				Y_vec = Y_vec, 
			Delta_vec = Delta_vec, 
			   Y_list = Y_list, 
		   Delta_list = Delta_list, 
			omega_mat = omegaikq_mat_prev,
		   beta00_vec = Thetagq_vec_updt[1 : dim_hzd],
		   beta0k_mat = matrix(Thetagq_vec_updt[(dim_hzd + 1) : (dim_hzd + dim_hzd * K)], dim_hzd, K, byrow = FALSE),
			beta0_vec = Thetagq_vec_updt[(dim_hzd + dim_hzd * K + 1) : (dim_hzd + dim_hzd * K + p)],
			betak_mat = matrix(Thetagq_vec_updt[(dim_hzd + dim_hzd * K + p + 1) : (dim_hzd + dim_hzd * K + p + p * K)], p, K, byrow = FALSE),
			  phi_vec = Thetagq_vec_updt[(dim_hzd + dim_hzd * K + p + p * K + 1) : (dim_hzd + dim_hzd * K + p + p * K + K)], 
				theta = exp(Thetagq_vec_updt[dim_hzd + dim_hzd * K + p + p * K + K + 1]),
		  hzd_bs_dist = hzd_bs_dist,
							CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL = CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL, 
						   CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT = CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT, 
									 CALCULATE_LOG_PRODUCT_COPULAS = CALCULATE_LOG_PRODUCT_COPULAS, 
			CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS) +
			log(apply(omegaikq_mat_temp, 1, dmvnorm, mean = rep(0, K), sigma = Sigma_mat)) - 
			log(apply(omegaikq_mat_prev, 1, dmvnorm, mean = rep(0, K), sigma = Sigma_mat))
			### n by 1
			
			compare_vec <- log(runif(n)) <= logARikq_vec
			### n by 1
			
			accept_rate_mat[, k] <- accept_rate_mat[, k] + 1 * compare_vec
			
			omegaikq_mat_updt[, k] <- ifelse(compare_vec, omegaikq_mat_temp[, k], omegaikq_mat_prev[, k])
			### n by 1
		}
		### omegaikq_mat_updt
		### n by K
		
		for(g in 1 : G){
		
			Thetagq_vec_temp <- Thetagq_vec_prev <- Thetagq_vec_updt
			### G by 1
			
			if(g %in% index_cond_1){
			### sigmak2, K by 1, para in (0, inf)
				repeat{
					Thetagq_vec_temp[g] <- Thetagq_vec_prev[g] + size_pstep_vec[g] * rnorm(1)
					### G by 1
					### only g-th elements are different
					if(Thetagq_vec_temp[g] > 0) break
				}
				logARgq_cor <- pnorm(Thetagq_vec_prev[g], log.p = TRUE) - pnorm(Thetagq_vec_temp[g], log.p = TRUE)
			} else if (g %in% index_cond_2){
			### rho, 1 by 1, para in (- 1, 1)
				repeat{
					Thetagq_vec_temp[g] <- Thetagq_vec_prev[g] + size_pstep_vec[g] * rnorm(1)
					### G by 1
					### only g-th elements are different
					if(Thetagq_vec_temp[g] > - 1 & Thetagq_vec_temp[g] < 1) break
				}
				logARgq_cor <- log(pnorm(1 - Thetagq_vec_prev[g]) - pnorm(- 1 - Thetagq_vec_prev[g])) - log(pnorm(1 - Thetagq_vec_temp[g]) - pnorm(- 1 - Thetagq_vec_temp[g]))
			} else {
				Thetagq_vec_temp[g] <- Thetagq_vec_prev[g] + size_pstep_vec[g] * rnorm(1)
				### G by 1
				### only g-th elements are different
				logARgq_cor <- 0
			}
			
			theta_prev <- exp(Thetagq_vec_prev[dim_hzd + dim_hzd * K + p + p * K + K + 1])
			sigma2_vec_prev <- Thetagq_vec_prev[(dim_hzd + dim_hzd * K + p + p * K + K + 1 + 1) : (dim_hzd + dim_hzd * K + p + p * K + K + 1 + K)]
			if(K == 1){
				Sigma_mat_prev <- as.matrix(sigma2_vec_prev)
			} else {
				rho_prev <- Thetagq_vec_prev[dim_hzd + dim_hzd * K + p + p * K + K + 1 + K + 1]
				sigma_vec_prev <- sqrt(sigma2_vec_prev)
				Sigma_mat_prev <- sigma_vec_prev %*% t(sigma_vec_prev) * rho_prev
				### K by K
				diag(Sigma_mat_prev) <- sigma2_vec_prev
			}

			theta_temp <- exp(Thetagq_vec_temp[dim_hzd + dim_hzd * K + p + p * K + K + 1])
			sigma2_vec_temp <- Thetagq_vec_temp[(dim_hzd + dim_hzd * K + p + p * K + K + 1 + 1) : (dim_hzd + dim_hzd * K + p + p * K + K + 1 + K)]
			if(K == 1){
				Sigma_mat_temp <- as.matrix(sigma2_vec_temp)
			} else {
				rho_temp <- Thetagq_vec_temp[dim_hzd + dim_hzd * K + p + p * K + K + 1 + K + 1]
				sigma_vec_temp <- sqrt(sigma2_vec_temp)
				Sigma_mat_temp <- sigma_vec_temp %*% t(sigma_vec_temp) * rho_temp
				### K by K
				diag(Sigma_mat_temp) <- sigma2_vec_temp
			}
			
			logARgq <- sum(VECTORIZED_LOG_LIKELIHOOD_FUNCTION(
			   	 x_mat = x_mat, 
				 M_vec = M_vec, 
				 Y_vec = Y_vec, 
			 Delta_vec = Delta_vec, 
	    		Y_list = Y_list, 
		    Delta_list = Delta_list, 
			 omega_mat = omegaikq_mat_updt,
			beta00_vec = Thetagq_vec_temp[1 : dim_hzd],
		    beta0k_mat = matrix(Thetagq_vec_temp[(dim_hzd + 1) : (dim_hzd + dim_hzd * K)], dim_hzd, K, byrow = FALSE),
			 beta0_vec = Thetagq_vec_temp[(dim_hzd + dim_hzd * K + 1) : (dim_hzd + dim_hzd * K + p)],
			 betak_mat = matrix(Thetagq_vec_temp[(dim_hzd + dim_hzd * K + p + 1) : (dim_hzd + dim_hzd * K + p + p * K)], p, K, byrow = FALSE),
			   phi_vec = Thetagq_vec_temp[(dim_hzd + dim_hzd * K + p + p * K + 1) : (dim_hzd + dim_hzd * K + p + p * K + K)], 
			     theta = exp(Thetagq_vec_temp[dim_hzd + dim_hzd * K + p + p * K + K + 1]),
		   hzd_bs_dist = hzd_bs_dist,
							CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL = CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL, 
						   CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT = CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT, 
									 CALCULATE_LOG_PRODUCT_COPULAS = CALCULATE_LOG_PRODUCT_COPULAS, 
			CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS)) - sum(VECTORIZED_LOG_LIKELIHOOD_FUNCTION(
				 x_mat = x_mat, 
				 M_vec = M_vec, 
			     Y_vec = Y_vec, 
			 Delta_vec = Delta_vec, 
			    Y_list = Y_list, 
		    Delta_list = Delta_list, 
			 omega_mat = omegaikq_mat_updt, 
			beta00_vec = Thetagq_vec_prev[1 : dim_hzd],
		    beta0k_mat = matrix(Thetagq_vec_prev[(dim_hzd + 1) : (dim_hzd + dim_hzd * K)], dim_hzd, K, byrow = FALSE),
			 beta0_vec = Thetagq_vec_prev[(dim_hzd + dim_hzd * K + 1) : (dim_hzd + dim_hzd * K + p)],
			 betak_mat = matrix(Thetagq_vec_prev[(dim_hzd + dim_hzd * K + p + 1) : (dim_hzd + dim_hzd * K + p + p * K)], p, K, byrow = FALSE),
			   phi_vec = Thetagq_vec_prev[(dim_hzd + dim_hzd * K + p + p * K + 1) : (dim_hzd + dim_hzd * K + p + p * K + K)], 
			     theta = exp(Thetagq_vec_prev[dim_hzd + dim_hzd * K + p + p * K + K + 1]),
		   hzd_bs_dist = hzd_bs_dist,
							CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL = CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL, 
						   CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT = CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT, 
									 CALCULATE_LOG_PRODUCT_COPULAS = CALCULATE_LOG_PRODUCT_COPULAS, 
	        CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS)) +
			sum(apply(omegaikq_mat_updt, 1, dmvnorm, mean = rep(0, K), sigma = Sigma_mat_temp, log = TRUE)) - 
			sum(apply(omegaikq_mat_updt, 1, dmvnorm, mean = rep(0, K), sigma = Sigma_mat_prev, log = TRUE)) + 
			### dgamma(theta_temp, shape = eps, rate = eps, log = TRUE) -
			### dgamma(theta_prev, shape = eps, rate = eps, log = TRUE) +
			sum(eps * log(eps) - lgamma(eps) - (eps + 1) * log(sigma2_vec_temp) - (eps / sigma2_vec_temp)) -
			sum(eps * log(eps) - lgamma(eps) - (eps + 1) * log(sigma2_vec_prev) - (eps / sigma2_vec_prev)) 
			
			if(K == 1){
				logARgq <- logARgq 
			} else {
				logARgq <- logARgq + dunif(rho_temp, min = - 1, max = 1, log = TRUE) - dunif(rho_prev, min = - 1, max = 1, log = TRUE)
			}
			### 1 by 1
			
			logARgq <- logARgq + logARgq_cor
			### 1 by 1
			
			compare_scale <- log(runif(1)) <= logARgq
			### 1 by 1
			
			accept_rate_vec[g] <- accept_rate_vec[g] + 1 * compare_scale
			### G by 1
			
			Thetagq_vec_updt[g] <- ifelse(compare_scale, Thetagq_vec_temp[g], Thetagq_vec_prev[g])
			### G by 1
		}
		Theta_est_mat[rep_main, ] <- Thetagq_vec_updt
		### G by 1
		rep_main <- rep_main + 1
	}
	
	return(list(
		omega_est_mat = omegaikq_mat_updt, 
	  accept_rate_mat = accept_rate_mat / rep_main_max,
		Theta_est_mat = Theta_est_mat,
	  accept_rate_vec = accept_rate_vec / rep_main_max
		   )
	)
}

EXTRACT_EVENT_FROM_LIST <- function(events, cohort, workingdat){
	
	### x_mat_temp: n by num of covariates
	### Y_vec_temp: n by 1
	### Delta_vec_temp: n by 1
	### Y_list_temp: n by Mi by K
	### Delta_list_temp: n by Mi by K
	### M_mat_temp: n by K
	
	x_mat_temp <- working_data$x_mat
	### 13720 by 10
	Y_vec_temp <- working_data$Y_vec
	### 13720 by 1
	Delta_vec_temp <- working_data$Delta_vec
	### 13720 by 1
	Y_list_temp <- working_data$Y_list 
	### 13720 by Mi by 3
	Delta_list_temp <- working_data$Delta_list
	### 13720 by Mi by 3
	M_mat_temp <- working_data$M_mat
	### 13720 by 3
	
	#################################
	### Step 1 extract the cohort ###
	#################################
	
	if(cohort == "FULL"){
		row_cohort_vec <- seq(nrow(x_mat_temp))
	} else if(cohort == "OLD"){
		row_cohort_vec <- which(x_mat_temp$AGE > 60)
	} else if(cohort == "CKD"){
		row_cohort_vec <- which(x_mat_temp$eGFR < 60)
	}
	x_mat <- x_mat_temp[row_cohort_vec, ]
	### n by 10
	x_mat[, c("SBP", "BMI", "AGE", "eGFR")] <- scale(x_mat[, c("SBP", "BMI", "AGE", "eGFR")])
	### n by 10		
	Y_vec <- Y_vec_temp[row_cohort_vec]
	### n by 1
	Delta_vec <- Delta_vec_temp[row_cohort_vec]
	### n by 1
	Y_list_temp <- Y_list_temp[row_cohort_vec]
	### n by Mi by 3
	Delta_list_temp <- Delta_list_temp[row_cohort_vec]
	### n by Mi by 3
	M_mat_temp <- M_mat_temp[row_cohort_vec, ]
	### n by 3
	
	#################################
	### Step 2 extract the events ###
	#################################
	
	### split form: our approach
	### long form: coxph() or frailtyPenal() approach for K = 1 only
	
	K <- length(events)
	Y_list <- Delta_list <- list()
	long_dat <- c()
	
	### case of K = 3
	if(K == 3){
		Y_list <- Y_list_temp
		Delta_list <- Delta_list_temp
		M_vec <- apply(M_mat_temp, 1, max)
		### order: MI, HF, SK
		
		long_dat_temp <- NULL
	} 
	
	### case of K = 2
	
	if(setequal(events, c("MI", "HF"))){
		M_vec <- apply(M_mat_temp[, c(1, 2)], 1, max)
		### order: MI, HF
	} else if (setequal(events, c("MI", "SK"))) {
		M_vec <- apply(M_mat_temp[, c(1, 3)], 1, max)
		### order: MI, SK
	} else if (setequal(events, c("HF", "SK"))) {
		M_vec <- apply(M_mat_temp[, c(2, 3)], 1, max)
		### order: HF, SK
	} else if (setequal(events, c("MI"))) {
		M_vec <- M_mat_temp[, c(1)]
		### order: MI
	} else if (setequal(events, c("HF"))) {
		M_vec <- M_mat_temp[, c(2)]
		### order: HF
	} else if (setequal(events, c("SK"))) {
		M_vec <- M_mat_temp[, c(3)]
		### order: SK
	}
	
	for(rep_row in 1 : length(Y_list_temp)){
		if(setequal(events, c("MI", "HF"))){
			Y_list[[rep_row]] <- Y_list_temp[[rep_row]][1 : M_vec[rep_row], c(1, 2), drop = FALSE]
			Delta_list[[rep_row]] <- Delta_list_temp[[rep_row]][1 : M_vec[rep_row], c(1, 2), drop = FALSE]
			### order: MI, HF
			
			long_dat <- NULL
			
		} else if (setequal(events, c("MI", "SK"))) {
			Y_list[[rep_row]] <- Y_list_temp[[rep_row]][1 : M_vec[rep_row], c(1, 3), drop = FALSE]
			Delta_list[[rep_row]] <- Delta_list_temp[[rep_row]][1 : M_vec[rep_row], c(1, 3), drop = FALSE]
			### order: MI, SK
			
			long_dat <- NULL
			
		} else if (setequal(events, c("HF", "SK"))) {
			Y_list[[rep_row]] <- Y_list_temp[[rep_row]][1 : M_vec[rep_row], c(2, 3), drop = FALSE]
			Delta_list[[rep_row]] <- Delta_list_temp[[rep_row]][1 : M_vec[rep_row], c(2, 3), drop = FALSE]
			### order: HF, SK
			
			long_dat <- NULL
			
		} else if (setequal(events, c("MI"))) {
			Y_list[[rep_row]] <- Y_list_temp[[rep_row]][1 : M_vec[rep_row], c(1), drop = FALSE]
			Delta_list[[rep_row]] <- Delta_list_temp[[rep_row]][1 : M_vec[rep_row], c(1), drop = FALSE]
			### order: MI
			
			long_dat_temp <- x_mat[rep_row, ]
			### 1 by 10
			long_dat_temp <- long_dat_temp[rep(seq_len(nrow(long_dat_temp)), each = M_vec[rep_row]), ]
			### Mi by 10
			long_dat_temp$GAP_TIME <- drop(Y_list[[rep_row]])
			### Mi by 1
			long_dat_temp$EVENT <- drop(Delta_list[[rep_row]])
			### Mi by 1
			long_dat_temp$DEATH <- c(rep(0, M_vec[rep_row] - 1), Delta_vec[rep_row])
			### Mi by 1
			long_dat_temp$VISIT_NUM <- seq(M_vec[rep_row])
			### Mi by 1
			long_dat <- rbind(long_dat, long_dat_temp)
			
		} else if (setequal(events, c("HF"))) {
			Y_list[[rep_row]] <- Y_list_temp[[rep_row]][1 : M_vec[rep_row], c(2), drop = FALSE]
			Delta_list[[rep_row]] <- Delta_list_temp[[rep_row]][1 : M_vec[rep_row], c(2), drop = FALSE]
			### order: HF
			
			long_dat_temp <- x_mat[rep_row, ]
			### 1 by 10
			long_dat_temp <- long_dat_temp[rep(seq_len(nrow(long_dat_temp)), each = M_vec[rep_row]), ]
			### Mi by 10
			long_dat_temp$GAP_TIME <- drop(Y_list[[rep_row]])
			### Mi by 1
			long_dat_temp$EVENT <- drop(Delta_list[[rep_row]])
			### Mi by 1
			long_dat_temp$DEATH <- c(rep(0, M_vec[rep_row] - 1), Delta_vec[rep_row])
			### Mi by 1
			long_dat_temp$VISIT_NUM <- seq(M_vec[rep_row])
			### Mi by 1
			long_dat <- rbind(long_dat, long_dat_temp)
			
		} else if (setequal(events, c("SK"))) {
			Y_list[[rep_row]] <- Y_list_temp[[rep_row]][1 : M_vec[rep_row], c(3), drop = FALSE]
			Delta_list[[rep_row]] <- Delta_list_temp[[rep_row]][1 : M_vec[rep_row], c(3), drop = FALSE]
			### order: SK
			
			long_dat_temp <- x_mat[rep_row, ]
			### 1 by 10
			long_dat_temp <- long_dat_temp[rep(seq_len(nrow(long_dat_temp)), each = M_vec[rep_row]), ]
			### Mi by 10
			long_dat_temp$GAP_TIME <- drop(Y_list[[rep_row]])
			### Mi by 1
			long_dat_temp$EVENT <- drop(Delta_list[[rep_row]])
			### Mi by 1
			long_dat_temp$DEATH <- c(rep(0, M_vec[rep_row] - 1), Delta_vec[rep_row])
			### Mi by 1
			long_dat_temp$VISIT_NUM <- seq(M_vec[rep_row])
			### Mi by 1
			long_dat <- rbind(long_dat, long_dat_temp)
		}
	}
	return(list(
		 x_mat = x_mat, 
		 Y_vec = Y_vec,
	 Delta_vec = Delta_vec, 
	    Y_list = Y_list, 
	Delta_list = Delta_list, 
		 M_vec = M_vec, 
	  long_dat = long_dat
	))
}

##############################
### Simulation Study Setup ###
##############################

DATA_FILE_NAME <- paste(INPUT_ADDRESS, "ARIC_RECURRENT_EVENT.RDS", sep = "/")
working_data <- readRDS(DATA_FILE_NAME)
res_dat <- EXTRACT_EVENT_FROM_LIST(events = events, cohort = cohort, workingdat = working_data)

x_mat_temp <- res_dat$x_mat
### n by 10
x_mat <- unname(as.matrix(subset(x_mat_temp, select = c(GEN, RACE, SBP, COM, DIA, SMK, BMI, AGE))))
### n by 8, exclude ID_C and eGFR
Y_vec <- res_dat$Y_vec
Delta_vec <- res_dat$Delta_vec
Y_list <- res_dat$Y_list
Delta_list <- res_dat$Delta_list
M_vec <- res_dat$M_vec
long_dat <- data.frame(res_dat$long_dat)

### length(M_vec[M_vec > 1]) / length(M_vec)
### FULL: 13744, MI: 28.5%, HF: 6%, SK: 7%, three: 34%
### OLD: 2565, MI: 40%, HF: 8%, SK: 9%, three: 47%
### CKD: 185, MI: 48%, HF: 5%, SK: 9%, three: 53%

K <- length(events)
n <- dim(x_mat)[1]
p <- dim(x_mat)[2]

### copula_gene <- "Clayton"
### "Clayton", "Frank", "Gumbel"
hzd_bs_dist <- "Exp"
### "Exp", "Weibull", "Gompertz", "Log-Logistic"

if(K == 1){
	### beta0_vec , p by 1
	### betak_mat , p by K
	### beta00    , 1 by 1
	### beta0k_vec, K by 1
	### hzd_bs_para_vec, (K + 1) by 1
	### phi_vec   , K by 1
	### sigma2_vec , K by 1
	### Theta_vec_truth <- c(beta00, beta0k_vec, beta0_vec, as.vector(betak_mat), phi_vec, logtheta, sigma2_vec)
	### 1 + K + p + pK + K + 1 + K
	
	G <- 1 + K + p + p * K + K + 1 + K
	
	size_pstep_vec <- c(
		3e-1, 									### beta00
		rep(3e-1, K), 							### beta0k_vec
		rep(2e-1, p), 							### beta0_vec
		rep(2.5e-1, p),						    ### betak_mat
		rep(2.5e-1, K),							### phi_vec
		2.5e-1, 								### logtheta
		rep(8e-2, K)							### sigma2_vec
	)
	
	index_cond_1 <- (1 + K + p + p * K + K + 1 + 1) : 
					(1 + K + p + p * K + K + 1 + K)
	### sigmak2, K by 1, para in (0, inf)
	
	index_cond_2 <- NULL
	### 1 by 1, para in (- 1, 1)
} else {
	### beta0_vec , p by 1
	### betak_mat , p by K
	### beta00    , 1 by 1
	### beta0k_vec, K by 1
	### hzd_bs_para_vec, (K + 1) by 1
	### phi_vec   , K by 1
	### sigma2_vec , K by 1
	### Theta_vec_truth <- c(beta00, beta0k_vec, beta0_vec, as.vector(betak_mat), phi_vec, logtheta, sigma2_vec, rho)
	### 1 + K + p + pK + K + 1 + K + 1
	
	G <- 1 + K + p + p * K + K + 1 + K + 1
	
	size_pstep_vec <- c(
		8e-1, 			  ### beta00
		rep(1e+0, K),     ### beta0k_vec
		rep(5e-2, p),     ### beta0_vec
		rep(5e-1, p * K), ### betak_mat
		rep(5e-2, K),	  ### phi_vec
		1e-2, 			  ### logtheta
		rep(5e-1, K),	  ### sigma2_vec
		1e-1 			  ### rho
	)
	
	index_cond_1 <- (1 + K + p + p * K + K + 1 + 1) : 
					(1 + K + p + p * K + K + 1 + K)
	### sigmak2, K by 1, para in (0, inf)
	
	index_cond_2 <- 1 + K + p + p * K + K + 1 + K + 1
	### rho, 1 by 1, para in (- 1, 1)
}
### use smaller size for a greater acceptation rate
### use greater size for a smaller acceptation rate
size_istep_mat <- matrix(rep(rep(8e-1, n), K), n, K, byrow = FALSE)
size_pstep_vec <- size_pstep_vec

if(initial == "random"){
	Theta_vec_ini <- rep(0.5, G)
	### G by 1
	
	dim_hzd <- 1
	sigma2_vec <- Theta_vec_ini[(dim_hzd + dim_hzd * K + p + p * K + K + 1 + 1) : (dim_hzd + dim_hzd * K + p + p * K + K + 1 + K)]
	### K by 1
	if(K == 1){
		Sigma_mat <- sqrt(sigma2_vec) %*% t(sqrt(sigma2_vec))
	} else {
		rho <- Theta_vec_ini[dim_hzd + dim_hzd * K + p + p * K + K + 1 + K + 1]
		### 1 by 1
		Sigma_mat <- sqrt(sigma2_vec) %*% t(sqrt(sigma2_vec)) * rho
		### K by K
		diag(Sigma_mat) <- sigma2_vec
	}
	
	omega_mat_ini <- rmvnorm(n = n, mean = rep(0, K), sigma = Sigma_mat)
	### n by K, remaining p.d.
} else if(initial == "given"){
	
	if(setequal(events, c("MI", "HF", "SK"))){
		
		Theta_vec_ini <- c(
			log(0.0002719967), 	 ### h00
			log(0.0001460534),### h0k
			log(7.313616E-06),### h0k
			log(3.314442E-05),### h0k
			0.7501423, -0.6245529, 0.5928256, 
			1.2350077,  1.1719920, 1.1620574,
			0.0148206,  0.6840326, 					 ### terminal
			0.7486732, -0.2668187, 0.5262305,
			1.0383464,  0.6842282, 0.6008014,
			0.0445354,  0.2238983, 				     ### MI
			1.146586,   1.053520,   1.249231,
			0,  	   -0.755979,  -0.720915,
			-0.397364, -0.316332, 				     ### HF
			1.146586,   1.053520,   1.249231,
			0,  	   -0.755979,  -0.720915,
			-0.397364, -0.316332, 				     ### SK
			0.849849, 					 			 ### phi
			-0.0464516, 					 		 ### phi
			0.635858, 								 ### phi
			0.5,								     ### logtheta
			1.19273, 					 			 ### sigma2
			3.9975, 					 			 ### sigma2
			1.90441, 					 			 ### sigma2
			0.5										 ### rho
		)
		
		Sigma_mat <- sqrt(c(1.19273, 3.9975, 1.90441)) %*% t(sqrt(c(1.19273, 3.9975, 1.90441))) * 0.5
		### K by K
		diag(Sigma_mat) <- c(1.19273, 3.9975, 1.90441)
		omega_mat_ini <- rmvnorm(n, rep(0, K), Sigma_mat)
		### n by K
		
	} else if (setequal(events, c("MI"))){
		Theta_vec_ini <- c(
			log(0.0002719967), 	 					 ### h00
			log(0.0001460534),						 ### h0k
			0.7501423, -0.6245529, 0.5928256, 
			1.2350077,  1.1719920, 1.1620574,
			0.0148206,  0.6840326, 					 ### terminal
			0.7486732, -0.2668187, 0.5262305,
			1.0383464,  0.6842282, 0.6008014,
			0.0445354,  0.2238983, 				     ### recurrent
			0.849849, 					 			 ### phi
			0.5,								     ### logtheta
			1.19273 					 			 ### sigma2
		)
		### G by 1
		omega_mat_ini <- as.matrix(rnorm(n * K, mean = 0, sd = sqrt(1.19273)))
		### n by K
	} else if (setequal(events, c("HF"))){
		Theta_vec_ini <- c(
			log(0.000267391), 	 					 ### h00
			log(7.313616E-06),						 ### h0k
			0.49552456, -0.49957036, 0.42583515, 
			0,           0.87552589, 0.70283481,
			0.00988959,  0.50106165, 				 ### terminal
			1.146586,      1.053520,   1.249231,
			0,  	      -0.755979,  -0.720915,
			-0.397364,    -0.316332, 				 ### recurrent
			-0.0464516, 					 		 ### phi
			0.5,								     ### logtheta
			3.9975 					 			     ### sigma2
		)
		### G by 1
		omega_mat_ini <- as.matrix(rnorm(n * K, mean = 0, sd = sqrt(3.9975)))
		### n by K
	} else if (setequal(events, c("SK"))){
		Theta_vec_ini <- c(
			log(0.0003441284), 	 					 ### h00
			log(3.314442E-05),						 ### h0k
			0.7228319, -0.6241645, 0.5814563, 
			1.2165999,  1.1421369, 1.0328898,
			-0.0617494, 0.6172256, 				     ### terminal
			0.0220007,  0.1051154, 0.5705223,
			-0.1774374, 1.4779496, 0.6266484,
			-0.4552599, 0.2596003, 				     ### recurrent
			0.635858, 					 		     ### phi
			0.5,								     ### logtheta
			1.90441 					 			 ### sigma2
		)
		### G by 1
		omega_mat_ini <- as.matrix(rnorm(n * K, mean = 0, sd = sqrt(1.90441)))
		### n by K
	}
}

time_bg <- Sys.time()

res_main <- SAMPLING_MAIN_BAYES( 
					   x_mat = x_mat, 
					   Y_vec = Y_vec, 
				   Delta_vec = Delta_vec, 
					  Y_list = Y_list, 
				  Delta_list = Delta_list, 
					   M_vec = M_vec, 
			   omega_mat_ini = omega_mat_ini, 
			   Theta_vec_ini = Theta_vec_ini, 
				rep_main_max = rep_main_max,			   
			  size_istep_mat = size_istep_mat, 
			  size_pstep_vec = size_pstep_vec, 
				index_cond_1 = index_cond_1, 
				index_cond_2 = index_cond_2, 
				 hzd_bs_dist = hzd_bs_dist,
						    CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL = CALCULATE_SURVIVAL_VECTOR_FOR_TERMINAL, 
						   CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT = CALCULATE_SURVIVAL_VECTOR_FOR_RECURRENT, 
									 CALCULATE_LOG_PRODUCT_COPULAS = CALCULATE_LOG_PRODUCT_COPULAS, 
			CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS = CACULATE_MULTI_DIM_DERIVATIVE_CLAYTON_SURVIVAL_COPULAS, 
								VECTORIZED_LOG_LIKELIHOOD_FUNCTION = VECTORIZED_LOG_LIKELIHOOD_FUNCTION)

time_ed <- Sys.time()

save.image(paste(INPUT_ADDRESS, paste("TEST", test_no, "Event", paste(events, collapse = "_"), "Cohort", paste(cohort, collapse = "_"), "Chain", rep_main_max, "ARIC_ReMRE.RData", sep = "_"), sep = "/"))

burn_in_main <- rep_main_max * (4 / 5)

Theta_est_mat <- res_main$Theta_est_mat
### rep_main_max by G
Theta_est_cl_mat <- Theta_est_mat[seq(from = burn_in_main + 1, to = nrow(Theta_est_mat), by = thinning_main), ]
### rep_main_max / 2 by G

Theta_est_mean <- apply(Theta_est_cl_mat, 2, mean)
### G by 1

time_ed - time_bg
### 1. view the comp. time

Theta_est_mean
### 2. view the changes

res_main$accept_rate_vec
### 3. view the acceptation rates for parameters

res_main$accept_rate_mat
### 4. view the acceptation rates for omega

t(apply(Theta_est_cl_mat, 2, function(x) {ci_temp <- ci(x, ci = 0.95, method = "HDI"); return(c(ci_temp$CI_low, ci_temp$CI_high))}))
### 5. view the Bayesian highest density intervals for parameters

#######################################
### apply the PWP model to the data ###
#######################################

res_PWP <- coxph(Surv(time = rep(0, nrow(long_dat)), time2 = GAP_TIME, event = EVENT) ~ GEN + RACE + SBP + COM + DIA + SMK + BMI + AGE + cluster(ID_C) + strata(VISIT_NUM), data = long_dat)
summary(res_PWP)

#################################
### apply the JFM to the data ###
#################################

res_FRAILTY <- frailtyPenal(
		formula = Surv(GAP_TIME, EVENT) ~ cluster(ID_C) + GEN + RACE + SBP + COM + DIA + SMK + BMI + AGE + terminal(DEATH), 
		formula.terminalEvent = ~ GEN + RACE + SBP + COM + DIA + SMK + BMI + AGE, 
		   data = long_dat, 
	   RandDist = "LogN", 
		n.knots = 8, 
		  kappa = c(9.55e+9, 1.41e+12),
	recurrentAG = FALSE
	)
summary(res_FRAILTY)











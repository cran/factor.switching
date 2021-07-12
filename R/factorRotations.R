#library('lpSolve')
#library('HDInterval')
#library('coda')
#library('bayesSurv')

rsp_exact <- function( lambda_mcmc, maxIter = 100, threshold = 1e-6, verbose=TRUE, rotate=TRUE, printIter=1000 ){
	lambda_mcmc <- as.matrix(lambda_mcmc)
	#with colnames: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q
	cnames <- colnames(lambda_mcmc)
	d <- dim(lambda_mcmc)[2]
	if(sum(grepl('LambdaV', cnames)) < d){
		stop('Column names of lambda_mcmc should be formatted as: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q, where p and q denote the number of variables and factors, respectively.')
	}
	q <- as.numeric(strsplit(cnames[d],split='_')[[1]][2])
	p <- d/q
	mcmcIterations <- dim(lambda_mcmc)[1]
	threshold <- threshold*mcmcIterations*p*q
	checkNames <- paste0(rep(paste0('LambdaV',1:p,'_'),each=q), 1:q)
	if(sum(checkNames==cnames) != d){stop('Column names of lambda_mcmc are wrong, please check input.')}

	lambda_mcmc_varimax <- lambda_mcmc
	if(rotate){
		if(q > 1){
			for(iter in 1:mcmcIterations){
				lambda <- matrix(lambda_mcmc[iter, ], nrow = p, byrow=T)
				v <- varimax(lambda, normalize=F)
				rotated_lambda <- v$loadings
				class(rotated_lambda) <- 'matrix'
				lambda_mcmc_varimax[iter, ] <- c(t(rotated_lambda))
			}
		}
	}

	all_c <- array(data = 1, dim = c(2^q, q))
	l <- 1
	if(q>1){
		for(i in 1:(q-1)){
			pos <- combn(1:q, i)
			j <- dim(pos)[2]
			for(k in 1:j){
				l <- l + 1
				all_c[l, pos[,k]] <- rep(-1, i)
			}
		}
	}	
	all_c[l+1, ] <- rep(-1, q)

	lambda_hat <-  matrix(colMeans(lambda_mcmc_varimax), ncol = q, nrow = p, byrow = TRUE)
	lambda_hat_zero <- lambda_hat
	c_vectors <- matrix( rep(1, q), ncol = q, nrow = mcmcIterations, byrow = TRUE)
	v_vectors <- matrix( 1:q, ncol = q, nrow = mcmcIterations, byrow = TRUE)
	st <- 1:q
	dim_all_c <- 2^q
	dim_all_v <- factorial(q)

	perm <- matrix( 1:q, ncol = q, nrow = dim_all_c, byrow = TRUE)
	costs <- numeric(dim_all_c)
	f <- numeric(maxIter)
	lambda_hat_values <- array(data = NA, dim = c(p, q, maxIter))
	totalIterations <- 0
	criterion = TRUE
	cost.matrix <- matrix(numeric(q * q), nrow = q, ncol = q)
	t1 <- numeric(maxIter)
	start_time <- Sys.time()
	while((criterion == TRUE)&(totalIterations < maxIter)){
		totalIterations <- totalIterations + 1
		cat(paste0('* iteration: ', totalIterations),'\n')
		lambda_hat_new <- 0*lambda_hat
		objective <- 0

		for(iter in 1:mcmcIterations){
			lambda <- matrix(lambda_mcmc_varimax[iter,], ncol = q, nrow = p, byrow = TRUE)
			for(i in 1:dim_all_c){
				c_vec <- all_c[i, ]
				lambda_switch <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda_hat
				for(j in 1:q){
					temp <- (lambda - lambda_switch[,j])^2
					cost.matrix[j, ] <- colSums(temp)
				}
				matr <- lp.assign(cost.matrix)$solution
				for (j in 1:q) {
					perm[i, j] <- st[matr[, j] > 0]
				}
				perm[i, ] <- order(perm[i, ])
				costs[i] <- sum(cost.matrix * matr)
			}
			minIndex <- order(costs)[1]
			v_vectors[iter, ] <- perm[minIndex, ]
			c_vectors[iter, ] <- all_c[minIndex, ]
			switchedMatrix <- matrix(c_vectors[iter, ], ncol = q, nrow = p, byrow=T) * lambda[ ,v_vectors[iter, ]]
			objective <- objective + costs[minIndex]
			lambda_hat_new <- lambda_hat_new + switchedMatrix
			#sum((switchedMatrix - lambda_hat)^2) should be equal to costs[minIndex]
			if((iter %% printIter == 0)&&(verbose==TRUE)){
				cat(paste0('          mcmc draw = ', iter, ':  sum f = ', round(objective,3)),'\r')}
		}

		f[totalIterations] <- objective
		if( totalIterations > 1){
			if( f[totalIterations - 1] - f[totalIterations] < threshold){criterion = FALSE}
		}
		if(verbose==TRUE){
			cat(paste0('   -----  objective function = ', round(objective,3), ' -----'),'\n')
			cat('\n')
		}
		lambda_hat_new <- lambda_hat_new/mcmcIterations
		lambda_hat <- lambda_hat_new
		lambda_hat_values[,,totalIterations] <- lambda_hat
		end_time <- Sys.time()
		t1[totalIterations] <- as.numeric(difftime(end_time, start_time, units='min'))

	}
	####compute f(0)
	c_vec <- rep(1,q)
	v_vec <- 1:q 
	f_zero <- 0
	for(i in 1:mcmcIterations){
		lambda <- matrix(lambda_mcmc_varimax[iter,], ncol = q, nrow = p, byrow = TRUE)
		switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda[ ,v_vec]
		f_zero <- f_zero + sum((switchedMatrix - lambda_hat)^2)
	}
	t_exact <- c(0, t1[1:totalIterations])
	f_exact <- c(f_zero,f[1:totalIterations])
	objective_function <- data.frame(time=t_exact, value=f_exact)
	#
	lambda_reordered_mcmc <- lambda_mcmc
	for(i in 1:mcmcIterations){
		lambda <- matrix(lambda_mcmc_varimax[i, ], ncol = q, nrow = p, byrow = TRUE)
		switchedMatrix <- matrix(c_vectors[i, ], ncol = q, nrow = p, byrow=T) * lambda[ ,v_vectors[i, ]] 
		lambda_reordered_mcmc[i, ] <- c(t(switchedMatrix))
	}
	lambda_reordered_mcmc <- as.mcmc(lambda_reordered_mcmc)
	result <- vector('list',length=5)
	result[[1]] <- lambda_reordered_mcmc
	result[[2]] <- c_vectors
	result[[3]] <- v_vectors
	result[[4]] <- lambda_hat
	result[[5]] <- objective_function
	names(result) <- c('lambda_reordered_mcmc', 'sign_vectors', 'permute_vectors', 'lambda_hat','objective_function')
	class(result) <- c('list', 'rsp')
	return(result)
	
}

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########

rsp_partial_sa <- function( lambda_mcmc, maxIter = 1000, threshold = 1e-6, verbose=TRUE, sa_loops, rotate=TRUE, increaseIter=FALSE, temperatureSchedule=NULL, printIter=1000 ){
	lambda_mcmc <- as.matrix(lambda_mcmc)
	#with colnames: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q
	cnames <- colnames(lambda_mcmc)
	d <- dim(lambda_mcmc)[2]
	if(sum(grepl('LambdaV', cnames)) < d){
		stop('Column names of lambda_mcmc should be formatted as: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q, where p and q denote the number of variables and factors, respectively.')
	}
	q <- as.numeric(strsplit(cnames[d],split='_')[[1]][2])
	p <- d/q
	mcmcIterations <- dim(lambda_mcmc)[1]
	threshold <- threshold*mcmcIterations*p*q
	checkNames <- paste0(rep(paste0('LambdaV',1:p,'_'),each=q), 1:q)
	if(sum(checkNames==cnames) != d){stop('Something is wrong, please check input.')}
	if(is.null(temperatureSchedule)){
		temperatureSchedule <- function(sim_ann_iter){1/log(sim_ann_iter+1)}
	}


	objective_function <- function(lambda, lambda_hat, c_vec, v_vec, q, p){
		switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda[ ,v_vec]
		f <- sum((switchedMatrix - lambda_hat)^2)
		result <- vector('list', length = 2)    
		result[[1]] <- f
		result[[2]] <- switchedMatrix
		names(result) <- c('value', 'matrix')
		return(result)
	}


	lambda_mcmc_varimax <- lambda_mcmc
	if(rotate){
		if(q > 1){
			for(iter in 1:mcmcIterations){
				lambda <- matrix(lambda_mcmc[iter, ], nrow = p, byrow=T)
				v <- varimax(lambda, normalize=F)
				rotated_lambda <- v$loadings
				class(rotated_lambda) <- 'matrix'
				lambda_mcmc_varimax[iter, ] <- c(t(rotated_lambda))
			}
		}
	}

	lambda_hat <-  matrix(colMeans(lambda_mcmc_varimax), ncol = q, nrow = p, byrow = TRUE)
	lambda_hat_zero <- lambda_hat

	st <- 1:q
	c_vectors <- matrix( rep(1, q), ncol = q, nrow = mcmcIterations, byrow = TRUE)
	v_vectors <- matrix( 1:q, ncol = q, nrow = mcmcIterations, byrow = TRUE)
	dim_all_c <- 2^q
	dim_all_v <- factorial(q)
	perm <- 1:q
	sim_an_iters <- sa_loops
	f <- numeric(maxIter)
	lambda_hat_values <- array(data = NA, dim = c(p, q, maxIter))
	totalIterations <- 0
	criterion = TRUE
	cost.matrix <- matrix(numeric(q * q), nrow = q, ncol = q)
	t1 <- numeric(maxIter)
	start_time <- Sys.time()
	while((criterion == TRUE)&(totalIterations < maxIter)){
		totalIterations <- totalIterations + 1
		cat(paste0('* iteration: ', totalIterations),'\n')
		lambda_hat_new <- 0*lambda_hat
		objective <- 0

		for(iter in 1:mcmcIterations){
			lambda <- matrix(lambda_mcmc_varimax[iter,], ncol = q, nrow = p, byrow = TRUE)
			f_values <- numeric(sim_an_iters)
		        c_vec <- c_vectors[iter, ]
		        v_vec <- v_vectors[iter, ]
			obj <- objective_function(lambda = lambda, lambda_hat = lambda_hat,  
					c_vec = c_vec, v_vec = v_vec, q = q, p = p)
			f_old <- obj$value
			f_values[1] <- f_old
			alpha <- 0.2
			temperature <- 1000
			accept_rate <- 1

			for(sim_ann_iter in 2:sim_an_iters){
	#			temperature <- temperature*alpha
				temperature <- temperatureSchedule(sim_ann_iter)
				if (runif(1) < 0.9){
				#	random switching of one sign
					c_proposed <- c_vec
					rIndex <- sample(q,1)
					c_proposed[rIndex] <- -c_vec[rIndex]
				}else{
				#	random proposal
					c_proposed <- sample(c(-1,1),q,replace=TRUE)
				}
				# 	computing v_proposed by solving the transportation problem with c_proposed
				lambda_switch <- matrix(c_proposed, ncol = q, nrow = p, byrow=T) * lambda_hat
				for(j in 1:q){
					temp <- (lambda - lambda_switch[,j])^2
					cost.matrix[j, ] <- colSums(temp)
				}
				matr <- lp.assign(cost.matrix)$solution
				for (j in 1:q) {
					perm[j] <- st[matr[, j] > 0]
				}
				v_proposed <- order(perm)
				cost_proposed <- sum(cost.matrix * matr)
				obj <- cost_proposed
				f_new <- obj

				if(f_new < f_old){
					c_vec <- c_proposed
					v_vec <- v_proposed
					f_old <- f_new
					accept_rate <- accept_rate + 1
					#switchedMatrix <- matrix(c_proposed, ncol = q, nrow = p, byrow=T) * lambda[ ,v_proposed]
				}else{
					u <- runif(1)
					ar <-  - (f_new - f_old)/temperature
					if(log(u) < ar){
						c_vec <- c_proposed
						v_vec <- v_proposed
						f_old <- f_new
						accept_rate <- accept_rate + 1
						#switchedMatrix <- matrix(c_proposed, ncol = q, nrow = p, byrow=T) * lambda[ ,v_proposed]
					} 
				}
				f_values[sim_ann_iter] <- f_old
			}
			#plot(f_values)
			v_vectors[iter, ] <- v_vec
			c_vectors[iter, ] <- c_vec
			switchedMatrix <- switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda[ ,v_vec]
			objective <- objective + f_old
			lambda_hat_new <- lambda_hat_new + switchedMatrix
			#sum((switchedMatrix - lambda_hat)^2) should be equal to costs[minIndex]
			if((iter %% printIter == 0)&&(verbose==TRUE)){
				cat(paste0('          mcmc draw = ', iter, ':  sum f = ', round(objective,3)),'\r')}
		}

		f[totalIterations] <- objective
		if( totalIterations > 10){
			f_diff <- f[totalIterations - 1] - f[totalIterations]
			if((f_diff < 0)&&(increaseIter)){
				sim_an_iters <- sim_an_iters + sa_loops
				cat(paste0('  WARNING: objective function increased.'),'\n')
				cat(paste0('           Simulated Annealing loops increased to ', sim_an_iters,'.'),'\n')
				f_diff <- threshold+1
			}
			if( abs(f_diff) < threshold){criterion = FALSE}		}
		if(verbose==TRUE){
			cat(paste0('   -----  objective function = ', round(objective,3), ' -----'),'\n')
			cat('\n')
		}
		lambda_hat_new <- lambda_hat_new/mcmcIterations
		lambda_hat <- lambda_hat_new
		lambda_hat_values[,,totalIterations] <- lambda_hat
		end_time <- Sys.time()
		t1[totalIterations] <- as.numeric(difftime(end_time, start_time, units='min'))

	}
	####compute f(0)
	c_vec <- rep(1,q)
	v_vec <- 1:q 
	f_zero <- 0
	for(i in 1:mcmcIterations){
		lambda <- matrix(lambda_mcmc_varimax[iter,], ncol = q, nrow = p, byrow = TRUE)
		switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda[ ,v_vec]
		f_zero <- f_zero + sum((switchedMatrix - lambda_hat)^2)
	}
	t_exact <- c(0, t1[1:totalIterations])
	f_exact <- c(f_zero,f[1:totalIterations])
	objective_function <- data.frame(time=t_exact, value=f_exact)
	#
	lambda_reordered_mcmc <- lambda_mcmc
	for(i in 1:mcmcIterations){
		lambda <- matrix(lambda_mcmc_varimax[i, ], ncol = q, nrow = p, byrow = TRUE)
		switchedMatrix <- matrix(c_vectors[i, ], ncol = q, nrow = p, byrow=T) * lambda[ ,v_vectors[i, ]] 
		lambda_reordered_mcmc[i, ] <- c(t(switchedMatrix))
	}
	lambda_reordered_mcmc <- as.mcmc(lambda_reordered_mcmc)
	result <- vector('list',length=5)
	result[[1]] <- lambda_reordered_mcmc
	result[[2]] <- c_vectors
	result[[3]] <- v_vectors
	result[[4]] <- lambda_hat
	result[[5]] <- objective_function
	names(result) <- c('lambda_reordered_mcmc', 'sign_vectors', 'permute_vectors', 'lambda_hat','objective_function')
	class(result) <- c('list', 'rsp')
	return(result)
}


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########

rsp_full_sa <- function( lambda_mcmc, maxIter = 1000, threshold = 1e-6, verbose=TRUE, sa_loops, rotate=TRUE, increaseIter=FALSE, temperatureSchedule=NULL, printIter=1000 ){
	lambda_mcmc <- as.matrix(lambda_mcmc)
	#with colnames: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q
	cnames <- colnames(lambda_mcmc)
	d <- dim(lambda_mcmc)[2]
	if(sum(grepl('LambdaV', cnames)) < d){
		stop('Column names of lambda_mcmc should be formatted as: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q, where p and q denote the number of variables and factors, respectively.')
	}
	q <- as.numeric(strsplit(cnames[d],split='_')[[1]][2])
	p <- d/q
	mcmcIterations <- dim(lambda_mcmc)[1]
	threshold <- threshold*mcmcIterations*p*q
	checkNames <- paste0(rep(paste0('LambdaV',1:p,'_'),each=q), 1:q)
	if(sum(checkNames==cnames) != d){stop('Something is wrong, please check input.')}
	if(is.null(temperatureSchedule)){
		temperatureSchedule <- function(sim_ann_iter){1/log(sim_ann_iter+1)}
	}

	objective_function <- function(lambda, lambda_hat, c_vec, v_vec, q, p){
		switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda[ ,v_vec]
		f <- sum((switchedMatrix - lambda_hat)^2)
		result <- vector('list', length = 2)    
		result[[1]] <- f
		result[[2]] <- switchedMatrix
		names(result) <- c('value', 'matrix')
		return(result)
	}


	lambda_mcmc_varimax <- lambda_mcmc
	if(rotate){
		if(q > 1){
			for(iter in 1:mcmcIterations){
				lambda <- matrix(lambda_mcmc[iter, ], nrow = p, byrow=T)
				v <- varimax(lambda, normalize=F)
				rotated_lambda <- v$loadings
				class(rotated_lambda) <- 'matrix'
				lambda_mcmc_varimax[iter, ] <- c(t(rotated_lambda))
			}
		}
	}

	lambda_hat <-  matrix(colMeans(lambda_mcmc_varimax), ncol = q, nrow = p, byrow = TRUE)
	lambda_hat_zero <- lambda_hat

	st <- 1:q
	c_vectors <- matrix( rep(1, q), ncol = q, nrow = mcmcIterations, byrow = TRUE)
	v_vectors <- matrix( 1:q, ncol = q, nrow = mcmcIterations, byrow = TRUE)
	dim_all_c <- 2^q
	dim_all_v <- factorial(q)
	perm <- 1:q
	sim_an_iters <- sa_loops
	f <- numeric(maxIter)
	lambda_hat_values <- array(data = NA, dim = c(p, q, maxIter))
	totalIterations <- 0
	criterion = TRUE
	cost.matrix <- matrix(numeric(q * q), nrow = q, ncol = q)
	t1 <- numeric(maxIter)
	start_time <- Sys.time()
	while((criterion == TRUE)&(totalIterations < maxIter)){
		totalIterations <- totalIterations + 1
		cat(paste0('* iteration: ', totalIterations),'\n')
		lambda_hat_new <- 0*lambda_hat
		objective <- 0

		for(iter in 1:mcmcIterations){
			lambda <- matrix(lambda_mcmc_varimax[iter,], ncol = q, nrow = p, byrow = TRUE)
			f_values <- numeric(sim_an_iters)
		        c_vec <- c_vectors[iter, ]
		        v_vec <- v_vectors[iter, ]
			obj <- objective_function(lambda = lambda, lambda_hat = lambda_hat, 
				c_vec = c_vec, v_vec = v_vec, q = q, p = p)
			f_old <- obj$value
			f_values[1] <- f_old
			alpha <- 0.2
			temperature <- 1000
			accept_rate <- 1

			for(sim_ann_iter in 2:sim_an_iters){
	#			temperature <- temperature*alpha
				temperature <- temperatureSchedule(sim_ann_iter)
				uR <- runif(1)
				if (uR < 0.8){
				#	both c and v
					c_proposed <- c_vec
					rIndex <- sample(q,1)
					c_proposed[rIndex] <- -c_vec[rIndex]
					random_pair <- sample(q,2,replace=F)
					v_proposed <- v_vec
					v_proposed[random_pair] <- v_vec[random_pair[c(2,1)]] 
				}else{
					if(uR < 0.9){
				#	only c
						c_proposed <- c_vec
						rIndex <- sample(q,1)
						c_proposed[rIndex] <- -c_vec[rIndex]
						v_proposed <- v_vec
					}else{
				#	only v
						c_proposed <- c_vec
						random_pair <- sample(q,2,replace=F)
						v_proposed <- v_vec
						v_proposed[random_pair] <- v_vec[random_pair[c(2,1)]] 
					}
				}
				f_new <- objective_function(lambda = lambda, lambda_hat = lambda_hat, 
						c_vec = c_proposed, v_vec = v_proposed, q = q, p = p)$value
				if(f_new < f_old){
					c_vec <- c_proposed
					v_vec <- v_proposed
					f_old <- f_new
					accept_rate <- accept_rate + 1
					#switchedMatrix <- matrix(c_proposed, ncol = q, nrow = p, byrow=T) * lambda[ ,v_proposed]
				}else{
					u <- runif(1)
					ar <-  - (f_new - f_old)/temperature
					if(log(u) < ar){
						c_vec <- c_proposed
						v_vec <- v_proposed
						f_old <- f_new
						accept_rate <- accept_rate + 1
						#switchedMatrix <- matrix(c_proposed, ncol = q, nrow = p, byrow=T) * lambda[ ,v_proposed]
					} 
				}
				f_values[sim_ann_iter] <- f_old
			}
			#plot(f_values)
			v_vectors[iter, ] <- v_vec
			c_vectors[iter, ] <- c_vec
			switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda[ ,v_vec]
			objective <- objective + f_old
			lambda_hat_new <- lambda_hat_new + switchedMatrix
			#sum((switchedMatrix - lambda_hat)^2) should be equal to costs[minIndex]
			if((iter %% printIter == 0)&&(verbose==TRUE)){
				cat(paste0('          mcmc draw = ', iter, ':  sum f = ', round(objective,3)),'\r')}
		}

		f[totalIterations] <- objective
		if( totalIterations > 10){
			f_diff <- f[totalIterations - 1] - f[totalIterations]
			if((f_diff < 0)&&(increaseIter)){
				sim_an_iters <- sim_an_iters + sa_loops
				cat(paste0('  WARNING: objective function increased.'),'\n')
				cat(paste0('           Simulated Annealing loops increased to ', sim_an_iters,'.'),'\n')
				f_diff <- threshold+1
			}
			if( abs(f_diff) < threshold){criterion = FALSE}
		}
		if(verbose==TRUE){
			cat(paste0('   -----  objective function = ', round(objective,3), ' -----'),'\n')
			cat('\n')
		}
		lambda_hat_new <- lambda_hat_new/mcmcIterations
		lambda_hat <- lambda_hat_new
		lambda_hat_values[,,totalIterations] <- lambda_hat
		end_time <- Sys.time()
		t1[totalIterations] <- as.numeric(difftime(end_time, start_time, units='min'))

	}
	####compute f(0)
	c_vec <- rep(1,q)
	v_vec <- 1:q 
	f_zero <- 0
	for(i in 1:mcmcIterations){
		lambda <- matrix(lambda_mcmc_varimax[iter,], ncol = q, nrow = p, byrow = TRUE)
		switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow=T) * lambda[ ,v_vec]
		f_zero <- f_zero + sum((switchedMatrix - lambda_hat)^2)
	}
	t_exact <- c(0, t1[1:totalIterations])
	f_exact <- c(f_zero,f[1:totalIterations])
	objective_function <- data.frame(time=t_exact, value=f_exact)
	#
	lambda_reordered_mcmc <- lambda_mcmc
	for(i in 1:mcmcIterations){
		lambda <- matrix(lambda_mcmc_varimax[i, ], ncol = q, nrow = p, byrow = TRUE)
		switchedMatrix <- matrix(c_vectors[i, ], ncol = q, nrow = p, byrow=T) * lambda[ ,v_vectors[i, ]] 
		lambda_reordered_mcmc[i, ] <- c(t(switchedMatrix))
	}
	lambda_reordered_mcmc <- as.mcmc(lambda_reordered_mcmc)
	result <- vector('list',length=5)
	result[[1]] <- lambda_reordered_mcmc
	result[[2]] <- c_vectors
	result[[3]] <- v_vectors
	result[[4]] <- lambda_hat
	result[[5]] <- objective_function
	names(result) <- c('lambda_reordered_mcmc', 'sign_vectors', 'permute_vectors', 'lambda_hat','objective_function')
	class(result) <- c('list', 'rsp')
	return(result)
}


procrustes_switching <- function (lambda_mcmc, maxIter = 100, threshold = 1e-06, verbose = TRUE, 
    rotate = FALSE, printIter = 1000) 
{
    lambda_mcmc <- as.matrix(lambda_mcmc)
    cnames <- colnames(lambda_mcmc)
    d <- dim(lambda_mcmc)[2]
    if (sum(grepl("LambdaV", cnames)) < d) {
        stop("Column names of lambda_mcmc should be formatted as: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q, where p and q denote the number of variables and factors, respectively.")
    }
    q <- as.numeric(strsplit(cnames[d], split = "_")[[1]][2])
    p <- d/q
    mcmcIterations <- dim(lambda_mcmc)[1]
    threshold <- threshold * mcmcIterations * p * q
    checkNames <- paste0(rep(paste0("LambdaV", 1:p, "_"), each = q), 
        1:q)
    if (sum(checkNames == cnames) != d) {
        stop("Column names of lambda_mcmc are wrong, please check input.")
    }
    lambda_mcmc_varimax <- lambda_mcmc
    if (rotate) {
        if (q > 1) {
            for (iter in 1:mcmcIterations) {
                lambda <- matrix(lambda_mcmc[iter, ], nrow = p, 
                  byrow = T)
                v <- varimax(lambda, normalize = F)
                rotated_lambda <- v$loadings
                class(rotated_lambda) <- "matrix"
                lambda_mcmc_varimax[iter, ] <- c(t(rotated_lambda))
            }
        }
    }
    lambda_hat <- matrix(colMeans(lambda_mcmc_varimax), ncol = q, 
        nrow = p, byrow = TRUE)
    lambda_hat_zero <- lambda_hat
    f <- numeric(maxIter)
    lambda_hat_values <- array(data = NA, dim = c(p, q, maxIter))
    totalIterations <- 0
    criterion = TRUE
    t1 <- numeric(maxIter)
    lambda_reordered_mcmc <- lambda_mcmc
    start_time <- Sys.time()
    while ((criterion == TRUE) & (totalIterations < maxIter)) {
        totalIterations <- totalIterations + 1
        cat(paste0("* iteration: ", totalIterations), "\n")
        lambda_hat_new <- 0 * lambda_hat
        objective <- 0
        for (iter in 1:mcmcIterations) {
            lambda <- matrix(lambda_mcmc_varimax[iter, ], ncol = q, 
                nrow = p, byrow = TRUE)
            switchedMatrix <- procrustes(lambda, lambda_hat)$X.new
            objective <- objective + sum((switchedMatrix - lambda_hat)^2)
            lambda_hat_new <- lambda_hat_new + switchedMatrix
            lambda_reordered_mcmc[iter, ] <- c(t(switchedMatrix))
            if ((iter%%printIter == 0) && (verbose == TRUE)) {
                cat(paste0("          mcmc draw = ", iter, ":  sum f = ", 
                  round(objective, 3)), "\r")
            }
        }
        f[totalIterations] <- objective
        if (totalIterations > 1) {
            if (f[totalIterations - 1] - f[totalIterations] < 
                threshold) {
                criterion = FALSE
            }
        }
        if (verbose == TRUE) {
            cat(paste0("   -----  objective function = ", round(objective, 
                3), " -----"), "\n")
            cat("\n")
        }
        lambda_hat_new <- lambda_hat_new/mcmcIterations
        lambda_hat <- lambda_hat_new
        lambda_hat_values[, , totalIterations] <- lambda_hat
        end_time <- Sys.time()
        t1[totalIterations] <- as.numeric(difftime(end_time, 
            start_time, units = "min"))
    }
    c_vec <- rep(1, q)
    v_vec <- 1:q
    f_zero <- 0
   for (i in 1:mcmcIterations) {
        lambda <- matrix(lambda_mcmc_varimax[iter, ], ncol = q, 
            nrow = p, byrow = TRUE)
        switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow = T) * 
            lambda[, v_vec]
        f_zero <- f_zero + sum((switchedMatrix - lambda_hat)^2)
    }
    t_exact <- c(0, t1[1:totalIterations])
    f_exact <- c(f_zero, f[1:totalIterations])
    objective_function <- data.frame(time = t_exact, value = f_exact)
#	final rotation of all mcmc according to a varimax rotation based on lambda_hat
        if (q > 1) {
		overall_rotation <- varimax(lambda_hat, normalize = F)$rotmat 
		lambda_hat <- lambda_hat %*% overall_rotation
		for (iter in 1:mcmcIterations) {
			lambda <- matrix(lambda_reordered_mcmc[iter, ], nrow = p, 
			  byrow = T)
			rotated_lambda <- lambda %*% overall_rotation
			class(rotated_lambda) <- "matrix"
			lambda_reordered_mcmc[iter, ] <- c(t(rotated_lambda))
		}
        }




    lambda_reordered_mcmc <- as.mcmc(lambda_reordered_mcmc)
    result <- vector("list", length = 5)
    result[[1]] <- lambda_reordered_mcmc
    result[[2]] <- NULL
    result[[3]] <- NULL
    result[[4]] <- lambda_hat
    result[[5]] <- objective_function
    names(result) <- c("lambda_reordered_mcmc", "sign_vectors", 
        "permute_vectors", "lambda_hat", "objective_function")
    class(result) <- c("list", "rsp")
    return(result)
}




weighted_procrustes_switching <- function (lambda_mcmc, maxIter = 100, threshold = 1e-06, verbose = TRUE, 
    weight = TRUE, printIter = 1000) 
{
    lambda_mcmc <- as.matrix(lambda_mcmc)
    cnames <- colnames(lambda_mcmc)
    d <- dim(lambda_mcmc)[2]
    if (sum(grepl("LambdaV", cnames)) < d) {
        stop("Column names of lambda_mcmc should be formatted as: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q, where p and q denote the number of variables and factors, respectively.")
    }
    q <- as.numeric(strsplit(cnames[d], split = "_")[[1]][2])
    p <- d/q
    mcmcIterations <- dim(lambda_mcmc)[1]
    threshold <- threshold * mcmcIterations * p * q
    checkNames <- paste0(rep(paste0("LambdaV", 1:p, "_"), each = q), 
        1:q)
    if (sum(checkNames == cnames) != d) {
        stop("Column names of lambda_mcmc are wrong, please check input.")
    }
    lambda_mcmc_varimax <- lambda_mcmc
    lambda_hat <- matrix(colMeans(lambda_mcmc_varimax), ncol = q, 
        nrow = p, byrow = TRUE)
    lambda_hat_zero <- lambda_hat
    f <- numeric(maxIter)
    lambda_hat_values <- array(data = NA, dim = c(p, q, maxIter))
    totalIterations <- 0
    criterion = TRUE
    t1 <- numeric(maxIter)
    lambda_reordered_mcmc <- lambda_mcmc
    start_time <- Sys.time()
    w <- matrix(0, p, p)
    while ((criterion == TRUE) & (totalIterations < maxIter)) {
        totalIterations <- totalIterations + 1
        cat(paste0("* iteration: ", totalIterations), "\n")

	
	if (weight) {
		z <- array(0, dim=c(q, q, p))
		for (iter in 1:mcmcIterations) {
			lambda <- matrix(lambda_mcmc[iter, ], nrow = p, byrow = T)
			for (r in 1:p){
				u <- lambda[r,] - lambda_hat[r,]
				z[,,r] <- z[,,r] + u %*% t(u)
			}
		}
		z <- z/mcmcIterations
		for(r in 1:p){
			w[r, r] <- (det(z[,,r]))^{-1/q}
		}
	}



        weighted_lambda_hat <- w %*% lambda_hat    
        lambda_hat_new <- 0 * lambda_hat
        objective <- 0
        for (iter in 1:mcmcIterations) {
            lambda <- matrix(lambda_mcmc_varimax[iter, ], ncol = q, 
                nrow = p, byrow = TRUE)
            switchedMatrix <- procrustes(lambda, weighted_lambda_hat)$X.new
            objective <- objective + sum((switchedMatrix - lambda_hat)^2)
            lambda_hat_new <- lambda_hat_new + switchedMatrix
            lambda_reordered_mcmc[iter, ] <- c(t(switchedMatrix))
            if ((iter%%printIter == 0) && (verbose == TRUE)) {
                cat(paste0("          mcmc draw = ", iter, ":  sum f = ", 
                  round(objective, 3)), "\r")
            }
        }
        f[totalIterations] <- objective
        if (totalIterations > 1) {
            if (f[totalIterations - 1] - f[totalIterations] < 
                threshold) {
                criterion = FALSE
            }
        }
        if (verbose == TRUE) {
            cat(paste0("   -----  objective function = ", round(objective, 
                3), " -----"), "\n")
            cat("\n")
        }
        lambda_hat_new <- lambda_hat_new/mcmcIterations
        lambda_hat <- lambda_hat_new
        lambda_hat_values[, , totalIterations] <- lambda_hat
        end_time <- Sys.time()
        t1[totalIterations] <- as.numeric(difftime(end_time, 
            start_time, units = "min"))
    }
    c_vec <- rep(1, q)
    v_vec <- 1:q
    f_zero <- 0
   for (i in 1:mcmcIterations) {
        lambda <- matrix(lambda_mcmc_varimax[iter, ], ncol = q, 
            nrow = p, byrow = TRUE)
        switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow = T) * 
            lambda[, v_vec]
        f_zero <- f_zero + sum((switchedMatrix - lambda_hat)^2)
    }
    t_exact <- c(0, t1[1:totalIterations])
    f_exact <- c(f_zero, f[1:totalIterations])
    objective_function <- data.frame(time = t_exact, value = f_exact)
#	final rotation of all mcmc according to a varimax rotation based on lambda_hat
        if (q > 1) {
		overall_rotation <- varimax(lambda_hat, normalize = F)$rotmat 
		lambda_hat <- lambda_hat %*% overall_rotation
		for (iter in 1:mcmcIterations) {
			lambda <- matrix(lambda_reordered_mcmc[iter, ], nrow = p, 
			  byrow = T)
			rotated_lambda <- lambda %*% overall_rotation
			class(rotated_lambda) <- "matrix"
			lambda_reordered_mcmc[iter, ] <- c(t(rotated_lambda))
		}
        }




    lambda_reordered_mcmc <- as.mcmc(lambda_reordered_mcmc)
    result <- vector("list", length = 5)
    result[[1]] <- lambda_reordered_mcmc
    result[[2]] <- NULL
    result[[3]] <- NULL
    result[[4]] <- lambda_hat
    result[[5]] <- objective_function
    names(result) <- c("lambda_reordered_mcmc", "sign_vectors", 
        "permute_vectors", "lambda_hat", "objective_function")
    class(result) <- c("list", "rsp")
    return(result)
}




#' @export
plot.rsp <- function(x, prob = 0.99, myCol = c('red','blue'), mfrow=NULL, subSet = NULL, simCR=TRUE,HDI=TRUE, ...){
	if( 'rsp' %in% class(x) ){
		oldpar <- par(no.readonly = TRUE)
		on.exit(par(oldpar)) 
		q <- dim(x$lambda_hat)[2]
		p <- dim(x$lambda_hat)[1]
		if(is.null(subSet)){subSet=1:q}else{q <- length(subSet)}
		if(is.null(mfrow)){mfrow=c(1,q)}
		errb<-apply(x$lambda_reordered_mcmc, 2, function(y){hdi(y, credMass = prob)})
		cr <- credible.region(x$lambda_reordered_mcmc,probs=prob)[as.character(prob)][[1]]
		par(mfcol = c(1,q), mar = numeric(4), oma = c(4, 4, 4, 4), mgp = c(2, .6, 0))
		epsilon = 0.2
		par(mfrow=mfrow)
		jindex<-0
		for(j in subSet){
			jindex<-jindex+1
		        plot(x$lambda_hat[,j],
		                ylim = c(-1.25,1.25),
		                type = 'n',pch=16,axes=FALSE)
		        if(jindex%%mfrow[2]==1){
		                axis(2L,cex.axis = 1.3)
		        }
		        if(jindex < q/mfrow[1]+1){
				axis(3L,cex.axis = 1.3)
		        }
		        if(jindex > (mfrow[1]-1)*mfrow[2]){
				axis(1L,cex.axis = 1.3)
		        }
		        if(jindex %% mfrow[2] == 0){
				axis(4L,cex.axis = 1.3)
		        }
		        box()
		        abline(h=0, lwd = 1.5, col = 'gray')
		        for(i in 1:p){
				if(simCR==TRUE){
				        points(c(i,i),cr[,paste0('LambdaV',i,'_',j)],type='l',col=myCol[1],lwd=1.3)
				        points(c(i - epsilon,i+epsilon), rep(cr[,paste0('LambdaV',i,'_',j)][1],2),
				                type='l',col=myCol[1],lwd=1.3)
				        points(c(i - epsilon,i+epsilon), rep(cr[,paste0('LambdaV',i,'_',j)][2],2),
				                type='l',col=myCol[1],lwd=1.3)
				}
				if(HDI==TRUE){
				        points(c(i,i), errb[,paste0('LambdaV',i,'_',j)],type='l',lty=3,col=myCol[2],lwd=1.3)
				        points(c(i - epsilon,i+epsilon), rep(errb[,paste0('LambdaV',i,'_',j)][1],2),
				                type='l',col=myCol[2],lwd=1.3)
				        points(c(i - epsilon,i+epsilon), rep(errb[,paste0('LambdaV',i,'_',j)][2],2),
				                type='l',col=myCol[2],lwd=1.3)
				}
		        }
		        points(x$lambda_hat[,j],pch=16,col='black')
		        legend('bottomright',lty=0,as.expression(bquote(italic(j)*' = '*.(j))),bty='n',cex=1.2)
		}
	#       mtext(bquote(italic(r)), side = 1, outer = TRUE, line = 2.2, cex = 1.5)
		mtext(bquote(ring(lambda)[italic(rj)]), side = 2, outer = TRUE, line = 1.8, cex = 1.0)
		mtext(bquote(italic(r)), side = 1, outer = TRUE, line = 1.8, cex = 1.0)
	}else{
		cat(paste("    The input is not in class `rsp`"),'\n')
	}
}

#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########

switch_and_permute <- function(lambda_mcmc, switch_vectors, permute_vectors){
	mcmcIterations <- dim(lambda_mcmc)[1]
	reordered_mcmc <- lambda_mcmc
	q <- dim(switch_vectors)[2]
	p <- length(lambda_mcmc[1,])/q
	for(i in 1:mcmcIterations){
	        lambda <- matrix(lambda_mcmc[i, ], ncol = q, nrow = p, byrow = TRUE)
	        switchedMatrix <- matrix(switch_vectors[i, ], ncol = q, nrow = p, byrow=T) * lambda[ ,permute_vectors[i, ]] 
	        reordered_mcmc[i, ] <- c(t(switchedMatrix))
	}
	return(reordered_mcmc)
}


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########

compareMultipleChains <- function(rspObjectList, scheme='exact', sa_loops=10, maxIter=1000, threshold=1e-6){
	nChains <- length(rspObjectList)

	cnames <- colnames(rspObjectList[[1]]$lambda_reordered_mcmc)
	d <- dim(rspObjectList[[1]]$lambda_reordered_mcmc)[2]
	q <- as.numeric(strsplit(cnames[d],split='_')[[1]][2])
	p <- d/q

	lambda_hat_values <- matrix(nrow=nChains,ncol=d)
	colnames(lambda_hat_values) <-  colnames(rspObjectList[[1]]$lambda_reordered_mcmc)
	for(i in 1:nChains){
		lambda_hat_values[i, ] <- c(t(rspObjectList[[i]]$lambda_hat))
	}
	if(scheme=='exact'){
		tankard <- rsp_exact( lambda_mcmc = lambda_hat_values, 
			maxIter = 100, threshold = threshold, verbose=TRUE, rotate=FALSE )
	}else{
		if(scheme=='partial_sa'){
			tankard <- rsp_partial_sa( lambda_mcmc = lambda_hat_values, maxIter = maxIter, 
				threshold = threshold, verbose=TRUE, sa_loops = sa_loops, rotate=FALSE )
		}else{
			tankard <- rsp_full_sa( lambda_mcmc = lambda_hat_values, maxIter = maxIter, 
				threshold = threshold, verbose=TRUE, sa_loops = sa_loops,rotate=FALSE )
		}
	}
	reorderedChains <- vector('list',length=nChains)
	for(i in 1:nChains){
		mcmcIterations <- dim(rspObjectList[[i]]$lambda_reordered_mcmc)[1]
		v_vectors <- matrix(tankard$permute_vectors[i,],nrow=mcmcIterations,ncol=q,byrow=TRUE)
		c_vectors <- matrix(tankard$sign_vectors[i,],nrow=mcmcIterations,ncol=q,byrow=TRUE)
		reorderedChains[[i]] <- as.mcmc(switch_and_permute(lambda_mcmc = rspObjectList[[i]]$lambda_reordered_mcmc, 
						switch_vectors = c_vectors, 
						permute_vectors = v_vectors))
	}
	reorderedChains <- as.mcmc.list(reorderedChains)
	return(reorderedChains)
}

credible.region <- function (sample, probs = c(0.9, 0.975)) {
    if (!length(dim(sample))) 
        stop("Incorrect sample parameter supplied")
    if (sum(probs < 0.5)) 
        stop("probs must be each at least 0.50")
    if (sum(probs >= 1)) 
        stop("probs must be strictly lower than 1")
    nn <- dim(sample)[1]
    p <- dim(sample)[2]
    k <- floor(nn * probs)
    ranks <- apply(sample, 2, rank, ties.method = "first")
    S.left <- nn + 1 - apply(ranks, 1, min)
    S.right <- apply(ranks, 1, max)
    SS <- apply(cbind(S.left, S.right), 1, max)
    SS <- SS[order(SS)]
    jstar <- SS[k]
    result <- list()
    for (kk in 1:length(jstar)) {
        up <- sample[ranks == jstar[kk]]
        low <- sample[ranks == nn + 1 - jstar[kk]]
        result[[kk]] <- rbind(low, up)
        rownames(result[[kk]]) <- c("Lower", "Upper")
        colnames(result[[kk]]) <- colnames(sample)
    }
    names(result) <- paste(probs)
    return(result)
}


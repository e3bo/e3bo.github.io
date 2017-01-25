SIR.onestep <- function (x, params) { #function to calculate one step of stochastic SIR
    X <- x[2] #local variable for susceptibles
    Y <- x[3] #local variable for infecteds
    Z <- x[4] #local variable for recovereds
    N <- X + Y + Z #total population size (subject to demographic change)
    with( #use with as in deterministic model to simplify code
        as.list(params), {
            rates <- c(mu * N, beta * X * Y / N, mu * X, mu * Y, gamma * Y, mu * Z)
            changes <- matrix(c( 1, 0, 0,
                                -1, 1, 0,
                                -1, 0, 0,
                                0,-1, 1,
                                0,-1, 0,
                                0, 0,-1), ncol = 3, byrow = TRUE)
            tau <- rexp(n = 1,rate = sum(rates)) # exponential waiting time
            U <- runif(1) #uniform random deviate
            m <- min(which(cumsum(rates) >= U * sum(rates)))
            x <- x[2:4] + changes[m,]
            return(out <- c(tau, x))
        }
    )
}

SIR.model <- function (x, params) { #function to simulate stochastic SIR
    output <- list() #set up list to store results
    output[[1]] <- x #first record of output is initial condition
    i <- 1
    while (x[3] > 0) { #iterate until epidemic ends
        output[[i + 1]] <- x <- SIR.onestep(x, params)
        i <- i + 1
    }
    output <- do.call(rbind, output)
    colnames(output) <- c("time", "X", "Y", "Z") #name variables
    output #return output
}

get_analytic_prob_fun <- function(R){
    if (R == 1L){
        function(n) {
            1 / (sqrt(n * pi))
        }
    } else {
        function(n, trunc = 500) {
            cum <- 0
            for (x in seq(n, trunc)){
                lterm <- (x - 1) * log(R) - (2 * x - 1) * log(R + 1)
                lterm <- lterm + lfactorial(2 * x - 2) - lfactorial(x)
                lterm <- lterm - lfactorial(x - 1) 
                cum <- cum + exp(lterm)
            }
            cum
        }
    }
}

make_outbreak_plots <- function(R = 0.5, nsims = 10){
    if (R < 0 | nsims < 1 ){
        stop("Invalid parameters")
    }
    if (R > 1){
        warning("Major outbreaks possible, simulations could take a while ...")
    }
    pop.size <- 1e5 #total population size
    Y0 <- 1 #initial number infected
    X0 <- pop.size - Y0 #initial number susceptible (~98% of population)
    xstart <- c(time = 0, X = X0, Y = Y0, Z = pop.size - X0 - Y0) #initial conditions
    params <- list(mu = 0.0000, beta = R, gamma = 1) #parameters
    data <- vector(mode = 'list', length = nsims) #initialize list to store the output
    for (k in 1:nsims) { #simulate nsims times
        data[[k]] <- as.data.frame(SIR.model(xstart, params))
        data[[k]]$cum.time <- cumsum(data[[k]]$time)
    }
    max.time <- max(sapply(data, function(d) max(d$cum.time)))
    max.y <- max(sapply(data, function(d) max(d$Y))) * 1.08
    ## ^^^ find max infected in all runs and increase by 8% for plot
    par(mfrow=c(2, 1)) # make a plot with 2 panels
    plot(Y ~ cum.time, data = data[[1]], xlab='\nTime', ylab='No. infected', col=1,
         xlim=c(0, max.time), ylim=c(0, max.y), type = 's')
    for (k in 1:nsims) { #add multiple epidemics to plot
        lines(Y ~ cum.time, data = data[[k]], col = k, type = 's')
    }
    # now plot the cumulative distribution of outbreak sizes
    nsteps <- sapply(data, nrow) - 1
    ob_sizes <- (nsteps + 1) / 2
    F <- ecdf(ob_sizes)
    knts <- knots(F)
    to <- knts[length(knts) - 1]
    if(to == 1){
        to <- 1.01
    }
    curve(1 - F(x), from=1, to = to, log="xy", n = 1001,
          xlab = "\nX", ylab = "Pr(outbreak size > X)")
    analytic_Fc <- Vectorize(get_analytic_prob_fun(R = R))
    curve(analytic_Fc(x), from=1, to = to, log="xy", n = 1001,
          add = TRUE, col = 2)    
}
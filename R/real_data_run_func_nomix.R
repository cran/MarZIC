real_data_run_func_nomix <-
  function(yi_vec,
           obs_m_vec,
           xi_vec,
           li_vec,
           confound_mat,
           x4_inter,
           x5_inter,
           k) {
    ini_value <- ini_bound_nomix(yi_vec, obs_m_vec, xi_vec, k)
    theta1 <- ini_value[[1]]
    lb_est <- ini_value[[2]]
    ub_est <- ini_value[[3]]

    if (length(confound_mat) > 0) {
      num_confound <- ncol(confound_mat)
      theta1 <- c(theta1, rep(0, num_confound * 3))
      lb_est <-
        c(lb_est, rep(-Inf, num_confound), rep(-10, num_confound * 2))
      ub_est <-
        c(ub_est, rep(Inf, num_confound), rep(10, num_confound * 2))
      conf_ind <- TRUE
    } else {
      num_confound <- 1
      theta1 <- c(theta1, rep(0, num_confound * 3))
      lb_est <-
        c(lb_est, rep(0, num_confound), rep(0, num_confound * 2))
      ub_est <-
        c(ub_est, rep(0, num_confound), rep(0, num_confound * 2))
      confound_mat <- matrix(0, nrow = length(yi_vec), ncol = 1)
      conf_ind <- FALSE
    }
    
    if (!x4_inter) {
      theta1[5]<-0
    }
    
    if (!x5_inter) {
      theta1[6]<-0
    }
    
    t1 <- Sys.time()
    theta0 <- 0
    
    while(sum(abs(theta1 - theta0)) > 1e-3) {
      theta0 <- theta1
      est1 <- solnl(
        theta0,
        objfun = function(par) {
          # cat(par,"\n")
          # par[6]<-0
          Q_theta_cpp_nomix(
            par,
            para_vec_0_ori = theta0,
            yi_vec,
            m_star_vec = obs_m_vec,
            x_i_vec = xi_vec,
            l_i_vec = li_vec,
            confound_mat = confound_mat,
            x4_inter = x4_inter,
            x5_inter = x5_inter
          )
        },
        lb = lb_est,
        ub = ub_est
      )
      
      theta1 <- as.numeric(est1$par)
    }
    
    est1$par <- theta1
    
    hess_mat <- hessian(function(x) {
      # x[6]<-0
      Q_theta_cpp_nomix(
        x,
        as.numeric(est1$par),
        yi_vec,
        m_star_vec = obs_m_vec,
        x_i_vec = xi_vec,
        l_i_vec = li_vec,
        confound_mat = confound_mat,
        x4_inter = x4_inter,
        x5_inter = x5_inter
      )
    }, as.numeric(est1$par))
    
    Jac_mat <- jacobian(function(y) {
      # y[6]<-0
      grad(
        function(x) {
          # x[6]<-0
          Q_theta_cpp_nomix(
            x,
            y,
            yi_vec,
            m_star_vec = obs_m_vec,
            x_i_vec = xi_vec,
            l_i_vec = li_vec,
            confound_mat = confound_mat,
            x4_inter = x4_inter,
            x5_inter = x5_inter
          )
        },
        as.numeric(est1$par)
      )
    }, as.numeric(est1$par))
    
    hess_est <- hess_mat + Jac_mat
    
    if (!x4_inter) {
      diag(hess_est)[5]<-1
    }
    if (!x5_inter) {
      diag(hess_est)[6]<-1
    }
    


    t2 <- Sys.time()

    # est1$par_true<-par_true
    est1$time <- difftime(t2, t1, units = "hours")
    est1$hess_est <- hess_est


    est1$mediation_effect <-
      mediation_effect_cal(
        para_vec = as.numeric(est1$par),
        x_1 = 0,
        x_2 = 1,
        confound_mat = confound_mat,
        x4_inter = x4_inter,
        x5_inter = x5_inter
      )

    mediation_var <- function(x) {
      return(mediation_effect_cal(
        x,
        x_1 = 0,
        x_2 = 1,
        confound_mat = confound_mat,
        x4_inter = x4_inter,
        x5_inter = x5_inter
      ))
    }

    Med_jac <- jacobian(mediation_var, as.numeric(est1$par))



    est1$Med_jac <- Med_jac
    
    true_fn <- li_total_raw(
      est1$par,
      yi_vec = yi_vec,
      obs_m_vec = obs_m_vec,
      xi_vec = xi_vec,
      li_vec = li_vec,
      confound_mat = confound_mat
    )
    est1$fn <- true_fn
    
    BIC_est <- 2 * est1$fn + log(length(yi_vec)) * length(theta1)
    AIC_est <- 2 * est1$fn + 2 * length(theta1)

    est1$BIC_est <- BIC_est
    est1$AIC_est <- AIC_est

    if (conf_ind == FALSE) {
      col_exclude <-
        c((length(est1$par) - num_confound * 3 + 1):length(est1$par))
      NIE_sd <-
        try(sqrt(diag(
          est1$Med_jac[, -col_exclude] %*%
            solve(est1$hess_est[-col_exclude, -col_exclude]) %*%
            t(est1$Med_jac[, -col_exclude])
        )), TRUE)


      est1$NIE_sd <- NIE_sd
      est1$par_sd <-
        sqrt(diag(solve(est1$hess_est[-col_exclude, -col_exclude])))
    } else {
      NIE_sd <-
        try(sqrt(diag(est1$Med_jac %*% solve(est1$hess_est) %*% t(est1$Med_jac))), TRUE)


      est1$NIE_sd <- NIE_sd
      est1$par_sd <- sqrt(diag(solve(est1$hess_est)))
    }

    # cat("k ", k, " theta ", est1$par, "\n")

    return(est1)
  }

hypers_q <- function(phi_ref = 0.2, D_ref = 20, type = "exponential", alpha_phi = 1.4, sd_log_b = 1, delay_seq = NULL) {
    # phi ~ beta(alpha, beta)
    beta_phi <- alpha_phi * (1 - phi_ref) / phi_ref + (2 - 1 / phi_ref) / 3
    # b ~ lognormal(mean_log_b, sd_log_b)
    if (type == "exponential") {
      b_ref <-  - log(0.05 / (1-phi_ref)) / D_ref
    } else if (type == "gompertz") {
      b_ref <- - log(log(0.95) / log(phi_ref)) / D_ref
    }
    mean_log_b <- log(b_ref)
    # logit(phi) ~ normal(mean_logit_phi, sd_logit_phi)
    mean_logit_phi <- digamma(alpha_phi) - digamma(beta_phi)
    sd_logit_phi <- sqrt(trigamma(alpha_phi) + trigamma(beta_phi))

    if (!is.null(delay_seq)) {
        qfun <- ifelse(type == "exponential",
            function(x, y) 1 - (1 - x) * exp(- delay_seq * y),
            function(x, y) exp(log(x) * exp(- delay_seq * y))
        )

        n <- 4000
        phi <- rbeta(n, alpha_phi, beta_phi)
        b <- rlnorm(n, mean_log_b, sd_log_b)
        B <- do.call(rbind, purrr::map2(phi, b, qfun))

        par(mfrow = c(3, 1))
        phi_seq <- seq(0, 1, length = 1000)
        plot(phi_seq, dbeta(phi_seq, alpha_phi, beta_phi), type = "l", lwd = 2)
        phi_median <- (alpha_phi - 1/3) / (alpha_phi + beta_phi - 2/3)
        abline(v = phi_median, col = 2, lwd = 3)
        b_seq <- seq(0, 3, length = 500)
        plot(b_seq, dlnorm(b_seq, mean_log_b, sd_log_b), type = "l", lwd = 2)
        b_median <- exp(mean_log_b)
        abline(v = b_median, col = 2, lwd = 3)
        matplot(delay_seq, t(B), type = "l", lty = 1, col = rgb(0,0,0, 0.2))
        lines(delay_seq, qfun(phi_median, b_median), ylim = c(0, 1), lwd = 3, col = 2)
    }

    list(alpha_phi = alpha_phi, beta_phi = beta_phi,
        mean_log_b = mean_log_b, sd_log_b = sd_log_b,
        mean_logit_phi = mean_logit_phi, sd_logit_phi = sd_logit_phi)
}

# 5-State Markov Cohort Model Engine
# AD DMT Cost-Effectiveness Analysis
# Sevinc Elif Sen

annual_to_monthly <- function(p) 1 - (1 - p)^(1/12)

run_markov <- function(params, drug = "lecanemab", price_discount = 0, include_cg = FALSE) {
  n_cycles <- params$horizon_years * params$cycles_per_year
  monthly_dr <- (1 + params$discount_rate)^(1/12) - 1

  tp <- list(
    mci_mild   = annual_to_monthly(params$tp_mci_mild),
    mild_mod   = annual_to_monthly(params$tp_mild_mod),
    mod_sev    = annual_to_monthly(params$tp_mod_sev),
    mci_death  = annual_to_monthly(params$tp_mci_death),
    mild_death = annual_to_monthly(params$tp_mild_death),
    mod_death  = annual_to_monthly(params$tp_mod_death),
    sev_death  = annual_to_monthly(params$tp_sev_death)
  )

  build_tm <- function(hr = 1) {
    m <- matrix(0, 5, 5)
    m[1,2] <- tp$mci_mild * hr;  m[1,5] <- tp$mci_death
    m[1,1] <- 1 - m[1,2] - m[1,5]
    m[2,3] <- tp$mild_mod * hr;  m[2,5] <- tp$mild_death
    m[2,2] <- 1 - m[2,3] - m[2,5]
    m[3,4] <- tp$mod_sev * hr;   m[3,5] <- tp$mod_death
    m[3,3] <- 1 - m[3,4] - m[3,5]
    m[4,5] <- tp$sev_death;      m[4,4] <- 1 - m[4,5]
    m[5,5] <- 1
    m
  }

  if (drug == "lecanemab") {
    hr <- params$lec_hr; tx_months <- params$lec_tx_months
    drug_monthly <- params$lec_drug_annual * (1 - price_discount/100) / 12
    inf_monthly <- params$lec_infusion * params$lec_n_infusions / tx_months
    diag <- params$lec_diag; mon_monthly <- params$lec_monitoring / tx_months
    aria_e <- params$lec_aria_e; aria_h <- params$lec_aria_h; aria_c <- params$lec_aria_cost
  } else {
    hr <- params$don_hr; tx_months <- params$don_tx_months
    drug_monthly <- params$don_drug_total * (1 - price_discount/100) / tx_months
    inf_monthly <- params$don_infusion * params$don_n_infusions / tx_months
    diag <- params$don_diag; mon_monthly <- params$don_monitoring / tx_months
    aria_e <- params$don_aria_e; aria_h <- params$don_aria_h; aria_c <- params$don_aria_cost
  }

  aria_monthly <- (aria_e + aria_h) * aria_c / tx_months
  tm_soc <- build_tm(1); tm_tx <- build_tm(hr)

  utilities <- c(params$u_mci, params$u_mild, params$u_mod, params$u_sev, params$u_death)
  hc_monthly <- c(params$c_mci, params$c_mild, params$c_mod, params$c_sev, 0) / 12
  cg_monthly <- c(params$cg_mci, params$cg_mild, params$cg_mod, params$cg_sev, 0) / 12

  run_arm <- function(is_tx) {
    trace <- matrix(0, n_cycles + 1, 5)
    trace[1, ] <- c(params$prop_mci, 1 - params$prop_mci, 0, 0, 0)
    total_cost <- 0; total_qaly <- 0; total_ly <- 0

    for (t in 1:n_cycles) {
      df <- 1 / (1 + monthly_dr)^(t - 0.5)
      trace[t+1, ] <- trace[t, ] %*% if (is_tx && t <= tx_months) tm_tx else tm_soc
      state_prop <- (trace[t, ] + trace[t+1, ]) / 2

      cycle_cost <- sum(state_prop * hc_monthly)
      if (include_cg) cycle_cost <- cycle_cost + sum(state_prop * cg_monthly)
      if (is_tx && t <= tx_months) cycle_cost <- cycle_cost + drug_monthly + inf_monthly + mon_monthly + aria_monthly
      if (is_tx && t == 1) cycle_cost <- cycle_cost + diag
      total_cost <- total_cost + cycle_cost * df

      u_adj <- utilities
      if (is_tx && t <= tx_months) {
        u_adj[1:4] <- u_adj[1:4] - (aria_e * params$u_aria_decrement + aria_h * params$u_aria_decrement * 0.5) / tx_months
      }
      total_qaly <- total_qaly + sum(state_prop * u_adj) / 12 * df
      total_ly <- total_ly + sum(state_prop[1:4]) / 12 * df
    }
    list(cost = total_cost, qaly = total_qaly, ly = total_ly, trace = trace)
  }

  soc <- run_arm(FALSE); tx <- run_arm(TRUE)
  inc_cost <- tx$cost - soc$cost; inc_qaly <- tx$qaly - soc$qaly; inc_ly <- tx$ly - soc$ly
  icer <- if (inc_qaly > 0) inc_cost / inc_qaly else Inf

  list(soc_cost = soc$cost, tx_cost = tx$cost, inc_cost = inc_cost,
       soc_qaly = soc$qaly, tx_qaly = tx$qaly, inc_qaly = inc_qaly,
       soc_ly = soc$ly, tx_ly = tx$ly, inc_ly = inc_ly,
       icer = icer, soc_trace = soc$trace, tx_trace = tx$trace)
}

##############################################################################
# Transferable Health Economic Evaluation Framework for Alzheimer's
# Disease-Modifying Therapies — Integrated CEA & BIA Platform
#
# Author: Sevinç Elif Şen, MSc
# Independent Health Economist & HEOR Researcher
#
# Model: 5-state Markov cohort model | HSE payer perspective | HIQA guidelines
# Treatments: Lecanemab (CLARITY-AD) & Donanemab (TRAILBLAZER-ALZ 2)
# Country: Ireland (parameterised for extension to additional systems)
#
# Components:
#   Module 1: Cost-Effectiveness Analysis (CEA) with PSA
#   Module 2: Budget Impact Analysis (BIA) with population cascade
#   Module 3: Threshold Price Analysis
#   Module 4: Head-to-Head Comparison
#   Module 5: Cross-Country Framework (architecture ready)
##############################################################################

if (!require(shiny)) install.packages("shiny")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(scales)) install.packages("scales")
if (!require(DT)) install.packages("DT")
if (!require(shinythemes)) install.packages("shinythemes")
if (!require(plotly)) install.packages("plotly")

library(shiny)
library(ggplot2)
library(scales)
library(DT)
library(shinythemes)
library(plotly)

# ============================================================================
# COLOUR PALETTE
# ============================================================================
col_primary   <- "#1B4F72"
col_secondary <- "#2E86C1"
col_accent    <- "#E74C3C"
col_success   <- "#27AE60"
col_warning   <- "#F39C12"
col_light     <- "#EBF5FB"
col_lec       <- "#2E86C1"
col_don       <- "#E74C3C"
col_soc       <- "#7F8C8D"

state_cols <- c("MCI-AD" = "#27AE60", "Mild AD" = "#F39C12",
                "Moderate AD" = "#E67E22", "Severe AD" = "#E74C3C",
                "Death" = "#7F8C8D")

# ============================================================================
# DEFAULT PARAMETERS — IRELAND
# ============================================================================
defaults <- list(
  country = "Ireland",
  wtp = 45000,
  discount_rate = 0.04,
  horizon_years = 30,
  cycles_per_year = 12,

  # Population (AD epidemiology — Ireland)
  total_pop = 5150000,
  prev_65plus = 0.165,     # proportion 65+
  ad_prevalence = 0.07,    # among 65+
  early_stage_pct = 0.35,  # MCI + mild among AD
  diagnosed_pct = 0.55,    # diagnosed among early-stage
  amyloid_pos_pct = 0.70,  # amyloid-positive among diagnosed early
  eligible_pct = 0.80,     # meeting all DMT eligibility criteria

  # Annual transition probabilities (natural history)
  tp_mci_mild = 0.15,
  tp_mild_mod = 0.20,
  tp_mod_sev  = 0.25,
  tp_mci_death = 0.02,
  tp_mild_death = 0.04,
  tp_mod_death = 0.08,
  tp_sev_death = 0.25,

  # Starting distribution
  prop_mci = 0.70,

  # Utilities (EQ-5D, Landeiro et al. 2020)
  u_mci  = 0.73,
  u_mild = 0.68,
  u_mod  = 0.54,
  u_sev  = 0.38,
  u_death = 0.00,
  u_aria_decrement = 0.05,

  # Health state costs (EUR 2024, HSE payer perspective)
  c_mci  = 4200,
  c_mild = 12800,
  c_mod  = 28500,
  c_sev  = 52000,

  # Caregiver costs (societal perspective, optional)
  cg_mci  = 2000,
  cg_mild = 8500,
  cg_mod  = 18000,
  cg_sev  = 12000,

  # --- Lecanemab ---
  lec_hr = 0.69,
  lec_tx_months = 18,
  lec_drug_annual = 24766,
  lec_infusion = 385,
  lec_n_infusions = 26,
  lec_diag = 3200,
  lec_monitoring = 2100,
  lec_aria_e = 0.126,
  lec_aria_h = 0.173,
  lec_aria_cost = 1500,


  # --- Donanemab ---
  don_hr = 0.73,
  don_tx_months = 12,
  don_drug_total = 24500,
  don_infusion = 385,
  don_n_infusions = 13,
  don_diag = 3200,
  don_monitoring = 2450,
  don_aria_e = 0.240,
  don_aria_h = 0.313,
  don_aria_cost = 1500
)

# ============================================================================
# MARKOV MODEL ENGINE
# ============================================================================
annual_to_monthly <- function(p) 1 - (1 - p)^(1/12)

run_markov <- function(params, drug = "lecanemab", price_discount = 0, include_cg = FALSE) {

  n_cycles <- params$horizon_years * params$cycles_per_year
  monthly_dr <- (1 + params$discount_rate)^(1/12) - 1

  # Monthly transition probabilities (natural history)
  tp <- list(
    mci_mild  = annual_to_monthly(params$tp_mci_mild),
    mild_mod  = annual_to_monthly(params$tp_mild_mod),
    mod_sev   = annual_to_monthly(params$tp_mod_sev),
    mci_death = annual_to_monthly(params$tp_mci_death),
    mild_death = annual_to_monthly(params$tp_mild_death),
    mod_death = annual_to_monthly(params$tp_mod_death),
    sev_death = annual_to_monthly(params$tp_sev_death)
  )

  # Build transition matrix
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

  # Drug-specific parameters
  if (drug == "lecanemab") {
    hr <- params$lec_hr
    tx_months <- params$lec_tx_months
    drug_monthly <- params$lec_drug_annual * (1 - price_discount/100) / 12
    inf_cost <- params$lec_infusion
    n_inf <- params$lec_n_infusions
    diag <- params$lec_diag
    mon <- params$lec_monitoring
    aria_e <- params$lec_aria_e
    aria_h <- params$lec_aria_h
    aria_c <- params$lec_aria_cost
  } else {
    hr <- params$don_hr
    tx_months <- params$don_tx_months
    drug_monthly <- params$don_drug_total * (1 - price_discount/100) / tx_months
    inf_cost <- params$don_infusion
    n_inf <- params$don_n_infusions
    diag <- params$don_diag
    mon <- params$don_monitoring
    aria_e <- params$don_aria_e
    aria_h <- params$don_aria_h
    aria_c <- params$don_aria_cost
  }

  inf_monthly <- inf_cost * n_inf / tx_months
  mon_monthly <- mon / tx_months
  aria_total <- (aria_e + aria_h) * aria_c
  aria_monthly <- aria_total / tx_months

  tm_soc <- build_tm(1)
  tm_tx  <- build_tm(hr)

  utilities <- c(params$u_mci, params$u_mild, params$u_mod, params$u_sev, params$u_death)
  hc_monthly <- c(params$c_mci, params$c_mild, params$c_mod, params$c_sev, 0) / 12
  cg_monthly <- c(params$cg_mci, params$cg_mild, params$cg_mod, params$cg_sev, 0) / 12

  run_arm <- function(is_tx) {
    trace <- matrix(0, n_cycles + 1, 5)
    trace[1, ] <- c(params$prop_mci, 1 - params$prop_mci, 0, 0, 0)

    total_cost <- 0; total_qaly <- 0; total_ly <- 0

    for (t in 1:n_cycles) {
      df <- 1 / (1 + monthly_dr)^(t - 0.5)

      if (is_tx && t <= tx_months) {
        trace[t+1, ] <- trace[t, ] %*% tm_tx
      } else {
        trace[t+1, ] <- trace[t, ] %*% tm_soc
      }

      # Half-cycle correction
      state_prop <- (trace[t, ] + trace[t+1, ]) / 2

      # Costs
      cycle_cost <- sum(state_prop * hc_monthly)
      if (include_cg) cycle_cost <- cycle_cost + sum(state_prop * cg_monthly)
      if (is_tx && t <= tx_months) {
        cycle_cost <- cycle_cost + drug_monthly + inf_monthly + mon_monthly + aria_monthly
      }
      if (is_tx && t == 1) cycle_cost <- cycle_cost + diag
      total_cost <- total_cost + cycle_cost * df

      # QALYs
      u_adj <- utilities
      if (is_tx && t <= tx_months) {
        u_adj[1:4] <- u_adj[1:4] - (aria_e * params$u_aria_decrement +
                                      aria_h * params$u_aria_decrement * 0.5) / tx_months
      }
      total_qaly <- total_qaly + sum(state_prop * u_adj) / 12 * df

      # Life years
      total_ly <- total_ly + sum(state_prop[1:4]) / 12 * df
    }

    list(cost = total_cost, qaly = total_qaly, ly = total_ly, trace = trace)
  }

  soc <- run_arm(FALSE)
  tx  <- run_arm(TRUE)

  inc_cost <- tx$cost - soc$cost
  inc_qaly <- tx$qaly - soc$qaly
  inc_ly   <- tx$ly - soc$ly
  icer <- if (inc_qaly > 0) inc_cost / inc_qaly else Inf

  list(
    soc_cost = soc$cost, tx_cost = tx$cost, inc_cost = inc_cost,
    soc_qaly = soc$qaly, tx_qaly = tx$qaly, inc_qaly = inc_qaly,
    soc_ly = soc$ly, tx_ly = tx$ly, inc_ly = inc_ly,
    icer = icer,
    soc_trace = soc$trace, tx_trace = tx$trace
  )
}

# ============================================================================
# PSA ENGINE
# ============================================================================
run_psa <- function(params, drug = "lecanemab", n_sim = 1000, price_discount = 0, include_cg = FALSE) {
  results <- data.frame(inc_cost = numeric(n_sim), inc_qaly = numeric(n_sim), icer = numeric(n_sim))

  for (i in 1:n_sim) {
    p <- params

    # Sample parameters
    p$u_mci  <- min(1, max(0, rnorm(1, params$u_mci, 0.05)))
    p$u_mild <- min(1, max(0, rnorm(1, params$u_mild, 0.05)))
    p$u_mod  <- min(1, max(0, rnorm(1, params$u_mod, 0.05)))
    p$u_sev  <- min(1, max(0, rnorm(1, params$u_sev, 0.05)))

    p$c_mci  <- max(0, rnorm(1, params$c_mci, params$c_mci * 0.15))
    p$c_mild <- max(0, rnorm(1, params$c_mild, params$c_mild * 0.15))
    p$c_mod  <- max(0, rnorm(1, params$c_mod, params$c_mod * 0.15))
    p$c_sev  <- max(0, rnorm(1, params$c_sev, params$c_sev * 0.15))

    p$tp_mci_mild  <- min(0.5, max(0.01, rnorm(1, params$tp_mci_mild, 0.03)))
    p$tp_mild_mod  <- min(0.5, max(0.01, rnorm(1, params$tp_mild_mod, 0.03)))
    p$tp_mod_sev   <- min(0.5, max(0.01, rnorm(1, params$tp_mod_sev, 0.03)))

    if (drug == "lecanemab") {
      p$lec_hr <- min(0.99, max(0.3, rnorm(1, params$lec_hr, 0.08)))
      p$lec_drug_annual <- max(0, rnorm(1, params$lec_drug_annual * (1 - price_discount/100),
                                         params$lec_drug_annual * 0.1))
    } else {
      p$don_hr <- min(0.99, max(0.3, rnorm(1, params$don_hr, 0.08)))
      p$don_drug_total <- max(0, rnorm(1, params$don_drug_total * (1 - price_discount/100),
                                        params$don_drug_total * 0.1))
    }

    r <- tryCatch(run_markov(p, drug, 0, include_cg), error = function(e) NULL)
    if (!is.null(r)) {
      results$inc_cost[i] <- r$inc_cost
      results$inc_qaly[i] <- r$inc_qaly
      results$icer[i]     <- r$icer
    }
  }
  results[results$icer < Inf & !is.na(results$icer), ]
}

# ============================================================================
# BIA ENGINE
# ============================================================================
run_bia <- function(params, uptake_scenario = "base", years = 5) {

  pop_65plus <- params$total_pop * params$prev_65plus
  ad_pop <- pop_65plus * params$ad_prevalence
  early_pop <- ad_pop * params$early_stage_pct
  diag_pop <- early_pop * params$diagnosed_pct
  amyloid_pop <- diag_pop * params$amyloid_pos_pct
  eligible_pop <- amyloid_pop * params$eligible_pct

  # Uptake trajectories by scenario
  uptake <- switch(uptake_scenario,
    "conservative" = c(0.02, 0.05, 0.08, 0.12, 0.15),
    "base"         = c(0.05, 0.10, 0.18, 0.25, 0.30),
    "optimistic"   = c(0.10, 0.20, 0.30, 0.40, 0.50)
  )

  # Market share: Year 1 lecanemab only, then donanemab enters
  lec_share <- c(1.0, 0.60, 0.50, 0.45, 0.40)
  don_share <- 1 - lec_share

  bia <- data.frame(
    Year = 1:years,
    Eligible = round(eligible_pop),
    Uptake_Pct = uptake * 100,
    Treated = round(eligible_pop * uptake),
    Lec_Patients = round(eligible_pop * uptake * lec_share[1:years]),
    Don_Patients = round(eligible_pop * uptake * don_share[1:years])
  )

  # Cost per patient per year
  lec_annual <- params$lec_drug_annual + params$lec_infusion * params$lec_n_infusions / (params$lec_tx_months/12) +
                params$lec_monitoring / (params$lec_tx_months/12)
  don_annual <- params$don_drug_total / (params$don_tx_months/12) + params$don_infusion * params$don_n_infusions / (params$don_tx_months/12) +
                params$don_monitoring / (params$don_tx_months/12)

  bia$Drug_Cost <- bia$Lec_Patients * lec_annual + bia$Don_Patients * don_annual
  bia$Diag_Cost <- bia$Treated * params$lec_diag  # same diagnostics for both
  bia$Monitor_Cost <- bia$Lec_Patients * params$lec_monitoring / (params$lec_tx_months/12) +
                       bia$Don_Patients * params$don_monitoring / (params$don_tx_months/12)
  bia$Infusion_Cost <- bia$Lec_Patients * params$lec_infusion * params$lec_n_infusions / (params$lec_tx_months/12) +
                        bia$Don_Patients * params$don_infusion * params$don_n_infusions / (params$don_tx_months/12)
  bia$ARIA_Cost <- bia$Lec_Patients * (params$lec_aria_e + params$lec_aria_h) * params$lec_aria_cost +
                    bia$Don_Patients * (params$don_aria_e + params$don_aria_h) * params$don_aria_cost
  bia$Total_Annual <- bia$Drug_Cost + bia$Diag_Cost + bia$Monitor_Cost + bia$Infusion_Cost + bia$ARIA_Cost
  bia$Cumulative <- cumsum(bia$Total_Annual)

  bia
}

# ============================================================================
# TORNADO / SENSITIVITY
# ============================================================================
run_tornado <- function(params, drug = "lecanemab", price_discount = 0) {
  base <- run_markov(params, drug, price_discount)$icer

  param_ranges <- list(
    "Drug Cost"              = if(drug=="lecanemab") c("lec_drug_annual", 0.7, 1.3)
                                else c("don_drug_total", 0.7, 1.3),
    "HR (Progression)"       = if(drug=="lecanemab") c("lec_hr", 0.5, 0.9)
                                else c("don_hr", 0.5, 0.9),
    "MCI\u2192Mild TP"       = c("tp_mci_mild", 0.7, 1.3),
    "Mild\u2192Moderate TP"  = c("tp_mild_mod", 0.7, 1.3),
    "Discount Rate"          = c("discount_rate", 0.75, 1.25),
    "Utility MCI"            = c("u_mci", 0.9, 1.1),
    "Utility Mild"           = c("u_mild", 0.85, 1.15),
    "Cost Severe AD"         = c("c_sev", 0.7, 1.3),
    "Cost Moderate AD"       = c("c_mod", 0.7, 1.3),
    "Starting Prop MCI"      = c("prop_mci", 0.8, 1.2)
  )

  tornado_data <- data.frame()
  for (nm in names(param_ranges)) {
    pr <- param_ranges[[nm]]
    pname <- pr[1]; lo_mult <- as.numeric(pr[2]); hi_mult <- as.numeric(pr[3])

    p_lo <- params; p_hi <- params
    p_lo[[pname]] <- params[[pname]] * lo_mult
    p_hi[[pname]] <- params[[pname]] * hi_mult

    icer_lo <- tryCatch(run_markov(p_lo, drug, price_discount)$icer, error = function(e) base)
    icer_hi <- tryCatch(run_markov(p_hi, drug, price_discount)$icer, error = function(e) base)

    tornado_data <- rbind(tornado_data, data.frame(
      Parameter = nm,
      Low = min(icer_lo, icer_hi),
      High = max(icer_lo, icer_hi),
      Base = base,
      stringsAsFactors = FALSE
    ))
  }

  tornado_data$Range <- tornado_data$High - tornado_data$Low
  tornado_data <- tornado_data[order(tornado_data$Range, decreasing = TRUE), ]
  tornado_data
}


# ============================================================================
# UI
# ============================================================================
ui <- fluidPage(
  theme = shinytheme("flatly"),

  tags$head(tags$style(HTML(sprintf("
    body { font-family: 'Segoe UI', Tahoma, sans-serif; }
    .navbar-default { background-color: %s !important; border-color: %s !important; }
    .navbar-default .navbar-brand { color: white !important; font-weight: 700; font-size: 16px; }
    .navbar-default .navbar-nav > li > a { color: #dce6f0 !important; }
    .navbar-default .navbar-nav > .active > a { background-color: %s !important; color: white !important; }
    .metric-box { background: white; border-left: 4px solid %s; padding: 15px; margin-bottom: 15px;
                  border-radius: 4px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
    .metric-box .value { font-size: 28px; font-weight: 700; color: %s; }
    .metric-box .label { font-size: 12px; color: #7f8c8d; text-transform: uppercase; }
    .sidebar-panel { background: #f8f9fa; border-radius: 6px; padding: 15px; }
    .info-text { font-size: 11px; color: #95a5a6; margin-top: 5px; }
    hr.section-break { border-top: 2px solid %s; margin: 20px 0; }
  ", col_primary, col_primary, col_secondary, col_secondary, col_primary, col_light)))),

  navbarPage(
    title = "AD-DMT Health Economic Evaluation Framework",
    id = "main_nav",

    # ── TAB 1: CEA ──
    tabPanel("Cost-Effectiveness",
      sidebarLayout(
        sidebarPanel(width = 3, class = "sidebar-panel",
          h4("Model Parameters"),
          selectInput("cea_drug", "Treatment",
                      choices = c("Lecanemab" = "lecanemab", "Donanemab" = "donanemab"),
                      selected = "lecanemab"),
          hr(),
          h5("Clinical"),
          conditionalPanel("input.cea_drug == 'lecanemab'",
            sliderInput("lec_hr", "Hazard Ratio", 0.3, 0.99, defaults$lec_hr, 0.01),
            sliderInput("lec_tx", "Treatment Duration (months)", 6, 36, defaults$lec_tx_months, 1)
          ),
          conditionalPanel("input.cea_drug == 'donanemab'",
            sliderInput("don_hr", "Hazard Ratio", 0.3, 0.99, defaults$don_hr, 0.01),
            sliderInput("don_tx", "Treatment Duration (months)", 6, 36, defaults$don_tx_months, 1)
          ),
          hr(),
          h5("Economic"),
          sliderInput("price_discount", "Price Discount (%)", 0, 80, 0, 5),
          sliderInput("discount_rate", "Discount Rate (%)", 0, 8, defaults$discount_rate * 100, 0.5),
          numericInput("wtp", "WTP Threshold (\u20AC/QALY)", defaults$wtp, min = 0, step = 5000),
          hr(),
          h5("Perspective"),
          checkboxInput("include_cg", "Include Caregiver Costs (Societal)", FALSE),
          hr(),
          p(class = "info-text", "CLARITY-AD: HR = 0.69 | TRAILBLAZER-ALZ 2: HR = 0.73"),
          p(class = "info-text", "HIQA WTP threshold: \u20AC45,000/QALY")
        ),

        mainPanel(width = 9,
          fluidRow(
            column(3, uiOutput("metric_icer")),
            column(3, uiOutput("metric_inc_cost")),
            column(3, uiOutput("metric_inc_qaly")),
            column(3, uiOutput("metric_threshold_ratio"))
          ),
          tabsetPanel(
            tabPanel("Results", br(), DT::dataTableOutput("cea_table"), br(), plotOutput("cost_breakdown", height = "350px")),
            tabPanel("Markov Trace", br(), plotOutput("markov_soc", height = "300px"), plotOutput("markov_tx", height = "300px")),
            tabPanel("Sensitivity", br(), plotOutput("tornado", height = "450px")),
            tabPanel("CE Plane (PSA)", br(),
                     fluidRow(column(4, numericInput("n_psa", "PSA Iterations", 1000, min = 100, max = 10000, step = 100)),
                              column(4, actionButton("run_psa", "Run PSA", class = "btn-primary", style = "margin-top:25px;"))),
                     plotOutput("ce_plane", height = "450px"),
                     plotOutput("ceac", height = "350px")),
            tabPanel("Price Analysis", br(), plotOutput("price_curve", height = "400px"), br(), DT::dataTableOutput("price_table"))
          )
        )
      )
    ),

    # ── TAB 2: HEAD-TO-HEAD ──
    tabPanel("Comparison",
      fluidRow(
        column(12, h3("Lecanemab vs. Donanemab — Head-to-Head Comparison"),
               p("Both treatments compared against standard of care at current Irish list prices."))
      ),
      fluidRow(
        column(6, h4("Lecanemab (CLARITY-AD)", style = sprintf("color:%s;", col_lec)),
               DT::dataTableOutput("comp_lec_table")),
        column(6, h4("Donanemab (TRAILBLAZER-ALZ 2)", style = sprintf("color:%s;", col_don)),
               DT::dataTableOutput("comp_don_table"))
      ),
      hr(),
      fluidRow(column(12, plotOutput("comp_chart", height = "400px")))
    ),

    # ── TAB 3: BIA ──
    tabPanel("Budget Impact",
      sidebarLayout(
        sidebarPanel(width = 3, class = "sidebar-panel",
          h4("BIA Parameters"),
          selectInput("bia_scenario", "Uptake Scenario",
                      choices = c("Conservative" = "conservative", "Base Case" = "base", "Optimistic" = "optimistic"),
                      selected = "base"),
          sliderInput("bia_years", "Time Horizon (years)", 1, 10, 5, 1),
          hr(),
          h5("Population Cascade"),
          numericInput("bia_pop", "Total Population", defaults$total_pop, min = 100000, step = 100000),
          sliderInput("bia_prev65", "% Aged 65+", 5, 30, defaults$prev_65plus * 100, 0.5),
          sliderInput("bia_ad_prev", "AD Prevalence (65+) %", 1, 15, defaults$ad_prevalence * 100, 0.5),
          sliderInput("bia_early", "Early-Stage %", 10, 60, defaults$early_stage_pct * 100, 5),
          sliderInput("bia_diag", "Diagnosed %", 20, 90, defaults$diagnosed_pct * 100, 5),
          sliderInput("bia_amyloid", "Amyloid Confirmed %", 30, 90, defaults$amyloid_pos_pct * 100, 5),
          sliderInput("bia_eligible", "Meeting Eligibility %", 40, 100, defaults$eligible_pct * 100, 5),
          hr(),
          p(class = "info-text", "Population data: CSO Ireland 2024")
        ),
        mainPanel(width = 9,
          fluidRow(
            column(4, uiOutput("bia_eligible_n")),
            column(4, uiOutput("bia_treated_yr5")),
            column(4, uiOutput("bia_cumulative"))
          ),
          tabsetPanel(
            tabPanel("Annual Impact", br(), plotOutput("bia_annual", height = "400px")),
            tabPanel("Cost Decomposition", br(), plotOutput("bia_decomp", height = "400px")),
            tabPanel("Population Cascade", br(), plotOutput("bia_cascade", height = "350px")),
            tabPanel("Data Table", br(), DT::dataTableOutput("bia_table"))
          )
        )
      )
    ),

    # ── TAB 4: ABOUT ──
    tabPanel("About",
      fluidRow(column(8, offset = 2,
        h3("Transferable Health Economic Evaluation Framework for AD DMTs"),
        p("This interactive platform accompanies the manuscript:"),
        p(tags$b("\u015Een, S.E. (2026). A Transferable Health Economic Evaluation Framework for Alzheimer\u2019s
          Disease-Modifying Therapies Across Small and Medium Healthcare Systems: A Multi-Country
          Decision-Analytic Modelling Study with Interactive R Shiny Platform.")),
        hr(),
        h4("Model Specifications"),
        tags$ul(
          tags$li("5-state Markov cohort model (MCI-AD, Mild AD, Moderate AD, Severe AD, Death)"),
          tags$li("Monthly cycles, lifetime horizon (30 years), half-cycle correction"),
          tags$li("HSE payer perspective (societal perspective optional)"),
          tags$li("4% annual discount rate per HIQA guidelines"),
          tags$li("WTP threshold: \u20AC45,000/QALY (HIQA)"),
          tags$li("Probabilistic sensitivity analysis: 10,000 Monte Carlo iterations"),
          tags$li("Budget impact analysis: ISPOR Principles of Good Practice")
        ),
        h4("Clinical Data Sources"),
        tags$ul(
          tags$li("Lecanemab: CLARITY-AD trial (van Dyck et al., NEJM 2023)"),
          tags$li("Donanemab: TRAILBLAZER-ALZ 2 trial (Sims et al., JAMA 2023)")
        ),
        h4("Cost Data Sources"),
        tags$ul(
          tags$li("HSE Fair Deal scheme tariffs"),
          tags$li("HIQA reference costs"),
          tags$li("Connolly et al. 2014 (dementia costs Ireland)"),
          tags$li("Landeiro et al. 2020 (EQ-5D utilities)")
        ),
        hr(),
        h4("Author"),
        p("Sevin\u00e7 Elif \u015Een, MSc"),
        p("Independent Health Economist & HEOR Researcher"),
        p("Expertise: Decision-analytic modelling | HTA | R Shiny | Health Economics"),
        hr(),
        p(class = "info-text", "\u00A9 2026 Sevin\u00e7 Elif \u015Een. All rights reserved."),
        p(class = "info-text", "Framework designed for extension to additional country contexts.")
      ))
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {

  # Reactive params
  current_params <- reactive({
    p <- defaults
    p$lec_hr <- input$lec_hr
    p$lec_tx_months <- input$lec_tx
    p$don_hr <- input$don_hr
    p$don_tx_months <- input$don_tx
    p$discount_rate <- input$discount_rate / 100
    p$wtp <- input$wtp
    # BIA params
    p$total_pop <- input$bia_pop
    p$prev_65plus <- input$bia_prev65 / 100
    p$ad_prevalence <- input$bia_ad_prev / 100
    p$early_stage_pct <- input$bia_early / 100
    p$diagnosed_pct <- input$bia_diag / 100
    p$amyloid_pos_pct <- input$bia_amyloid / 100
    p$eligible_pct <- input$bia_eligible / 100
    p
  })

  # ── CEA Results ──
  cea_results <- reactive({
    run_markov(current_params(), input$cea_drug, input$price_discount, input$include_cg)
  })

  # Metric boxes
  output$metric_icer <- renderUI({
    r <- cea_results()
    tags$div(class = "metric-box",
      tags$div(class = "value", sprintf("\u20AC%s", format(round(r$icer), big.mark = ","))),
      tags$div(class = "label", "ICER (\u20AC/QALY)")
    )
  })
  output$metric_inc_cost <- renderUI({
    r <- cea_results()
    tags$div(class = "metric-box",
      tags$div(class = "value", sprintf("\u20AC%s", format(round(r$inc_cost), big.mark = ","))),
      tags$div(class = "label", "Incremental Cost")
    )
  })
  output$metric_inc_qaly <- renderUI({
    r <- cea_results()
    tags$div(class = "metric-box",
      tags$div(class = "value", sprintf("%.3f", r$inc_qaly)),
      tags$div(class = "label", "Incremental QALYs")
    )
  })
  output$metric_threshold_ratio <- renderUI({
    r <- cea_results()
    ratio <- r$icer / input$wtp
    tags$div(class = "metric-box",
      tags$div(class = "value", style = if(ratio > 1) "color:#E74C3C;" else "color:#27AE60;",
               sprintf("%.1fx", ratio)),
      tags$div(class = "label", "Times WTP Threshold")
    )
  })

  # Results table
  output$cea_table <- DT::renderDataTable({
    r <- cea_results()
    drug_label <- if(input$cea_drug == "lecanemab") "Lecanemab + SoC" else "Donanemab + SoC"
    df <- data.frame(
      Outcome = c("Total Costs (\u20AC)", "Total QALYs", "Total Life-Years", "ICER (\u20AC/QALY)"),
      Treatment = c(sprintf("\u20AC%s", format(round(r$tx_cost), big.mark = ",")),
                     sprintf("%.2f", r$tx_qaly), sprintf("%.2f", r$tx_ly), ""),
      SoC = c(sprintf("\u20AC%s", format(round(r$soc_cost), big.mark = ",")),
               sprintf("%.2f", r$soc_qaly), sprintf("%.2f", r$soc_ly), ""),
      Increment = c(sprintf("\u20AC%s", format(round(r$inc_cost), big.mark = ",")),
                     sprintf("%.3f", r$inc_qaly), sprintf("%.3f", r$inc_ly),
                     sprintf("\u20AC%s", format(round(r$icer), big.mark = ",")))
    )
    colnames(df)[2] <- drug_label
    datatable(df, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })

  # Cost breakdown chart
  output$cost_breakdown <- renderPlot({
    r <- cea_results()
    drug_label <- if(input$cea_drug == "lecanemab") "Lecanemab + SoC" else "Donanemab + SoC"
    df <- data.frame(
      Arm = c(drug_label, "SoC"),
      Cost = c(r$tx_cost, r$soc_cost)
    )
    ggplot(df, aes(x = Arm, y = Cost, fill = Arm)) +
      geom_col(width = 0.5) +
      scale_fill_manual(values = c(col_secondary, col_soc)) +
      scale_y_continuous(labels = label_comma(prefix = "\u20AC")) +
      labs(title = "Total Lifetime Costs", y = "Cost (\u20AC)", x = "") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })

  # Markov traces
  output$markov_soc <- renderPlot({
    r <- cea_results()
    trace <- r$soc_trace
    years <- (0:nrow(trace)[-1]) / 12
    states <- c("MCI-AD", "Mild AD", "Moderate AD", "Severe AD", "Death")
    df <- data.frame()
    for (i in 1:5) {
      df <- rbind(df, data.frame(Year = (0:(nrow(trace)-1))/12, Proportion = trace[,i], State = states[i]))
    }
    df$State <- factor(df$State, levels = states)
    ggplot(df, aes(x = Year, y = Proportion, fill = State)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = state_cols) +
      labs(title = "Markov Trace \u2014 Standard of Care", x = "Years", y = "Proportion") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })

  output$markov_tx <- renderPlot({
    r <- cea_results()
    trace <- r$tx_trace
    drug_label <- if(input$cea_drug == "lecanemab") "Lecanemab" else "Donanemab"
    states <- c("MCI-AD", "Mild AD", "Moderate AD", "Severe AD", "Death")
    df <- data.frame()
    for (i in 1:5) {
      df <- rbind(df, data.frame(Year = (0:(nrow(trace)-1))/12, Proportion = trace[,i], State = states[i]))
    }
    df$State <- factor(df$State, levels = states)
    ggplot(df, aes(x = Year, y = Proportion, fill = State)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = state_cols) +
      labs(title = paste("Markov Trace \u2014", drug_label, "+ SoC"), x = "Years", y = "Proportion") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })

  # Tornado
  output$tornado <- renderPlot({
    td <- run_tornado(current_params(), input$cea_drug, input$price_discount)
    td <- head(td, 10)
    td$Parameter <- factor(td$Parameter, levels = rev(td$Parameter))
    ggplot(td) +
      geom_segment(aes(x = Low, xend = High, y = Parameter, yend = Parameter),
                   linewidth = 8, color = col_secondary, alpha = 0.7) +
      geom_vline(xintercept = td$Base[1], linetype = "dashed", color = col_accent, linewidth = 1) +
      scale_x_continuous(labels = label_comma(prefix = "\u20AC")) +
      labs(title = "One-Way Sensitivity Analysis (Tornado Diagram)",
           subtitle = sprintf("Base ICER: \u20AC%s/QALY", format(round(td$Base[1]), big.mark = ",")),
           x = "ICER (\u20AC/QALY)", y = "") +
      theme_minimal(base_size = 13)
  })

  # PSA
  psa_data <- eventReactive(input$run_psa, {
    withProgress(message = "Running PSA...", {
      run_psa(current_params(), input$cea_drug, input$n_psa, input$price_discount, input$include_cg)
    })
  })

  output$ce_plane <- renderPlot({
    df <- psa_data()
    if (nrow(df) == 0) return(NULL)
    ggplot(df, aes(x = inc_qaly, y = inc_cost)) +
      geom_point(alpha = 0.3, color = col_secondary, size = 1.5) +
      geom_abline(slope = input$wtp, intercept = 0, linetype = "dashed", color = col_accent, linewidth = 1) +
      geom_hline(yintercept = 0, color = "grey40") +
      geom_vline(xintercept = 0, color = "grey40") +
      scale_y_continuous(labels = label_comma(prefix = "\u20AC")) +
      labs(title = "Cost-Effectiveness Plane",
           subtitle = sprintf("WTP threshold: \u20AC%s/QALY | %d iterations",
                             format(input$wtp, big.mark = ","), nrow(df)),
           x = "Incremental QALYs", y = "Incremental Cost (\u20AC)") +
      theme_minimal(base_size = 13)
  })

  output$ceac <- renderPlot({
    df <- psa_data()
    if (nrow(df) == 0) return(NULL)
    wtp_range <- seq(0, 200000, by = 5000)
    prob_ce <- sapply(wtp_range, function(w) mean(df$inc_cost / df$inc_qaly <= w, na.rm = TRUE))
    ceac_df <- data.frame(WTP = wtp_range, Probability = prob_ce)
    ggplot(ceac_df, aes(x = WTP, y = Probability)) +
      geom_line(color = col_secondary, linewidth = 1.2) +
      geom_vline(xintercept = input$wtp, linetype = "dashed", color = col_accent) +
      scale_x_continuous(labels = label_comma(prefix = "\u20AC")) +
      scale_y_continuous(labels = percent, limits = c(0, 1)) +
      labs(title = "Cost-Effectiveness Acceptability Curve",
           x = "Willingness-to-Pay (\u20AC/QALY)", y = "Probability Cost-Effective") +
      theme_minimal(base_size = 13)
  })

  # Price analysis
  output$price_curve <- renderPlot({
    discounts <- seq(0, 80, by = 5)
    icers <- sapply(discounts, function(d) run_markov(current_params(), input$cea_drug, d, input$include_cg)$icer)
    df <- data.frame(Discount = discounts, ICER = icers)
    ggplot(df, aes(x = Discount, y = ICER)) +
      geom_line(color = col_secondary, linewidth = 1.2) +
      geom_point(color = col_secondary, size = 2) +
      geom_hline(yintercept = input$wtp, linetype = "dashed", color = col_accent, linewidth = 1) +
      annotate("text", x = 70, y = input$wtp + 15000, label = sprintf("WTP = \u20AC%s",
               format(input$wtp, big.mark = ",")), color = col_accent, size = 4) +
      scale_y_continuous(labels = label_comma(prefix = "\u20AC")) +
      labs(title = "ICER vs. Drug Price Discount",
           x = "Price Discount (%)", y = "ICER (\u20AC/QALY)") +
      theme_minimal(base_size = 13)
  })

  output$price_table <- DT::renderDataTable({
    discounts <- seq(0, 80, by = 10)
    results <- lapply(discounts, function(d) {
      r <- run_markov(current_params(), input$cea_drug, d, input$include_cg)
      data.frame(
        `Discount (%)` = d,
        `ICER (\u20AC/QALY)` = sprintf("\u20AC%s", format(round(r$icer), big.mark = ",")),
        `Inc Cost (\u20AC)` = sprintf("\u20AC%s", format(round(r$inc_cost), big.mark = ",")),
        `Cost-Effective?` = if(r$icer <= input$wtp) "\u2705 Yes" else "\u274C No",
        check.names = FALSE
      )
    })
    datatable(do.call(rbind, results), options = list(dom = 't', paging = FALSE), rownames = FALSE)
  })

  # ── COMPARISON TAB ──
  comp_lec <- reactive({ run_markov(current_params(), "lecanemab", 0, input$include_cg) })
  comp_don <- reactive({ run_markov(current_params(), "donanemab", 0, input$include_cg) })

  output$comp_lec_table <- DT::renderDataTable({
    r <- comp_lec()
    df <- data.frame(
      Outcome = c("ICER (\u20AC/QALY)", "Inc. Cost (\u20AC)", "Inc. QALYs", "Inc. Life-Years", "Times WTP"),
      Value = c(sprintf("\u20AC%s", format(round(r$icer), big.mark = ",")),
                sprintf("\u20AC%s", format(round(r$inc_cost), big.mark = ",")),
                sprintf("%.3f", r$inc_qaly), sprintf("%.3f", r$inc_ly),
                sprintf("%.1fx", r$icer / input$wtp))
    )
    datatable(df, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })

  output$comp_don_table <- DT::renderDataTable({
    r <- comp_don()
    df <- data.frame(
      Outcome = c("ICER (\u20AC/QALY)", "Inc. Cost (\u20AC)", "Inc. QALYs", "Inc. Life-Years", "Times WTP"),
      Value = c(sprintf("\u20AC%s", format(round(r$icer), big.mark = ",")),
                sprintf("\u20AC%s", format(round(r$inc_cost), big.mark = ",")),
                sprintf("%.3f", r$inc_qaly), sprintf("%.3f", r$inc_ly),
                sprintf("%.1fx", r$icer / input$wtp))
    )
    datatable(df, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE)
  })

  output$comp_chart <- renderPlot({
    rl <- comp_lec(); rd <- comp_don()
    df <- data.frame(
      Drug = rep(c("Lecanemab", "Donanemab"), each = 2),
      Metric = rep(c("ICER (\u20AC/QALY)", "Incremental QALYs"), 2),
      Value = c(rl$icer, rl$inc_qaly * 100000, rd$icer, rd$inc_qaly * 100000)
    )
    # Side-by-side ICER comparison
    icer_df <- data.frame(
      Drug = c("Lecanemab", "Donanemab"),
      ICER = c(rl$icer, rd$icer)
    )
    ggplot(icer_df, aes(x = Drug, y = ICER, fill = Drug)) +
      geom_col(width = 0.5) +
      geom_hline(yintercept = input$wtp, linetype = "dashed", color = col_accent, linewidth = 1) +
      annotate("text", x = 1.5, y = input$wtp + 15000,
               label = sprintf("HIQA WTP: \u20AC%s", format(input$wtp, big.mark = ",")),
               color = col_accent, size = 4) +
      scale_fill_manual(values = c("Lecanemab" = col_lec, "Donanemab" = col_don)) +
      scale_y_continuous(labels = label_comma(prefix = "\u20AC")) +
      labs(title = "ICER Comparison: Lecanemab vs. Donanemab",
           subtitle = "Both vs. Standard of Care at list prices",
           y = "ICER (\u20AC/QALY)", x = "") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })

  # ── BIA TAB ──
  bia_data <- reactive({
    p <- current_params()
    run_bia(p, input$bia_scenario, input$bia_years)
  })

  output$bia_eligible_n <- renderUI({
    b <- bia_data()
    tags$div(class = "metric-box",
      tags$div(class = "value", format(b$Eligible[1], big.mark = ",")),
      tags$div(class = "label", "Eligible Patients")
    )
  })
  output$bia_treated_yr5 <- renderUI({
    b <- bia_data()
    tags$div(class = "metric-box",
      tags$div(class = "value", format(tail(b$Treated, 1), big.mark = ",")),
      tags$div(class = "label", paste0("Treated (Year ", nrow(b), ")"))
    )
  })
  output$bia_cumulative <- renderUI({
    b <- bia_data()
    tags$div(class = "metric-box",
      tags$div(class = "value", sprintf("\u20AC%sM", format(round(tail(b$Cumulative, 1) / 1e6, 1), big.mark = ","))),
      tags$div(class = "label", "Cumulative Budget Impact")
    )
  })

  output$bia_annual <- renderPlot({
    b <- bia_data()
    ggplot(b, aes(x = Year, y = Total_Annual / 1e6)) +
      geom_col(fill = col_secondary, width = 0.6) +
      geom_line(aes(y = Cumulative / 1e6), color = col_accent, linewidth = 1.2) +
      geom_point(aes(y = Cumulative / 1e6), color = col_accent, size = 3) +
      scale_y_continuous(labels = label_comma(suffix = "M", prefix = "\u20AC")) +
      labs(title = sprintf("Budget Impact \u2014 %s Scenario",
                          tools::toTitleCase(input$bia_scenario)),
           x = "Year", y = "Cost (\u20AC Millions)",
           caption = "Bars = annual | Line = cumulative") +
      theme_minimal(base_size = 14)
  })

  output$bia_decomp <- renderPlot({
    b <- bia_data()
    decomp <- data.frame(
      Year = rep(b$Year, 4),
      Category = rep(c("Drug Acquisition", "Diagnostics", "Monitoring & Infusion", "ARIA Management"),
                     each = nrow(b)),
      Cost = c(b$Drug_Cost, b$Diag_Cost, b$Monitor_Cost + b$Infusion_Cost, b$ARIA_Cost) / 1e6
    )
    decomp$Category <- factor(decomp$Category,
                              levels = c("Drug Acquisition", "Diagnostics", "Monitoring & Infusion", "ARIA Management"))
    ggplot(decomp, aes(x = Year, y = Cost, fill = Category)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c(col_primary, col_secondary, col_warning, col_accent)) +
      scale_y_continuous(labels = label_comma(suffix = "M", prefix = "\u20AC")) +
      labs(title = "Budget Impact by Cost Category", x = "Year", y = "Cost (\u20AC Millions)") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })

  output$bia_cascade <- renderPlot({
    p <- current_params()
    cascade <- data.frame(
      Stage = c("Total Population", "Aged 65+", "With AD",
                "Early-Stage", "Diagnosed", "Amyloid+", "Eligible"),
      N = c(p$total_pop,
            p$total_pop * p$prev_65plus,
            p$total_pop * p$prev_65plus * p$ad_prevalence,
            p$total_pop * p$prev_65plus * p$ad_prevalence * p$early_stage_pct,
            p$total_pop * p$prev_65plus * p$ad_prevalence * p$early_stage_pct * p$diagnosed_pct,
            p$total_pop * p$prev_65plus * p$ad_prevalence * p$early_stage_pct * p$diagnosed_pct * p$amyloid_pos_pct,
            p$total_pop * p$prev_65plus * p$ad_prevalence * p$early_stage_pct * p$diagnosed_pct * p$amyloid_pos_pct * p$eligible_pct)
    )
    cascade$Stage <- factor(cascade$Stage, levels = rev(cascade$Stage))
    ggplot(cascade, aes(x = Stage, y = N)) +
      geom_col(fill = col_secondary, width = 0.6) +
      geom_text(aes(label = format(round(N), big.mark = ",")), hjust = -0.1, size = 3.5) +
      scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.3))) +
      coord_flip() +
      labs(title = "Population Cascade \u2014 DMT Eligibility", x = "", y = "Number of Patients") +
      theme_minimal(base_size = 13)
  })

  output$bia_table <- DT::renderDataTable({
    b <- bia_data()
    b$Drug_Cost <- sprintf("\u20AC%s", format(round(b$Drug_Cost), big.mark = ","))
    b$Diag_Cost <- sprintf("\u20AC%s", format(round(b$Diag_Cost), big.mark = ","))
    b$Monitor_Cost <- sprintf("\u20AC%s", format(round(b$Monitor_Cost), big.mark = ","))
    b$Infusion_Cost <- sprintf("\u20AC%s", format(round(b$Infusion_Cost), big.mark = ","))
    b$ARIA_Cost <- sprintf("\u20AC%s", format(round(b$ARIA_Cost), big.mark = ","))
    b$Total_Annual <- sprintf("\u20AC%s", format(round(b$Total_Annual), big.mark = ","))
    b$Cumulative <- sprintf("\u20AC%s", format(round(b$Cumulative), big.mark = ","))
    datatable(b, options = list(dom = 't', paging = FALSE), rownames = FALSE)
  })
}

# ============================================================================
# RUN
# ============================================================================
shinyApp(ui = ui, server = server)

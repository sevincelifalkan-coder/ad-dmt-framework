# Default Parameters - Ireland
# HSE Payer Perspective | HIQA Guidelines
# Sevinc Elif Sen

get_ireland_params <- function() {
  list(
    country = "Ireland",
    wtp = 45000,
    discount_rate = 0.04,
    horizon_years = 30,
    cycles_per_year = 12,

    # Population (CSO Ireland 2024)
    total_pop = 5150000,
    prev_65plus = 0.165,
    ad_prevalence = 0.07,
    early_stage_pct = 0.35,
    diagnosed_pct = 0.55,
    amyloid_pos_pct = 0.70,
    eligible_pct = 0.80,

    # Transition probabilities (annual, natural history)
    tp_mci_mild = 0.15,
    tp_mild_mod = 0.20,
    tp_mod_sev  = 0.25,
    tp_mci_death = 0.02,
    tp_mild_death = 0.04,
    tp_mod_death = 0.08,
    tp_sev_death = 0.25,

    prop_mci = 0.70,

    # Utilities (Landeiro et al. 2020)
    u_mci = 0.73, u_mild = 0.68, u_mod = 0.54, u_sev = 0.38,
    u_death = 0.00, u_aria_decrement = 0.05,

    # Health state costs EUR 2024
    c_mci = 4200, c_mild = 12800, c_mod = 28500, c_sev = 52000,

    # Caregiver costs (societal)
    cg_mci = 2000, cg_mild = 8500, cg_mod = 18000, cg_sev = 12000,

    # Lecanemab (CLARITY-AD)
    lec_hr = 0.69, lec_tx_months = 18, lec_drug_annual = 24766,
    lec_infusion = 385, lec_n_infusions = 26, lec_diag = 3200,
    lec_monitoring = 2100, lec_aria_e = 0.126, lec_aria_h = 0.173,
    lec_aria_cost = 1500,

    # Donanemab (TRAILBLAZER-ALZ 2)
    don_hr = 0.73, don_tx_months = 12, don_drug_total = 24500,
    don_infusion = 385, don_n_infusions = 13, don_diag = 3200,
    don_monitoring = 2450, don_aria_e = 0.240, don_aria_h = 0.313,
    don_aria_cost = 1500
  )
}

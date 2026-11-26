# =============================================================================
# 4EK453 Principy ekonometrie – Přednáška 9: Endogenita a instrumentální proměnné
# Jan Zouhar, Katedra ekonometrie FIS VŠE v Praze
#
# Skript reprodukuje všechny výsledky z prezentace. Spouštějte ho z kořene
# složky kurzu (data se čtou z podsložky data/).
# Řádky „# ---- název ----“ označují bloky, které se zobrazují v prezentaci.
#
# Potřebné balíčky: AER (funkce ivreg), sandwich, lmtest
#   install.packages(c("AER", "sandwich", "lmtest"))
# =============================================================================

# ---- nastaveni ----
options(scipen = 10)     # raději víc číslic než vědecký zápis
suppressPackageStartupMessages(library(AER))  # ivreg(): IV a 2SLS
library(sandwich)        # robustní směrodatné chyby: vcovHC()
library(lmtest)          # coeftest()


# =============================================================================
# Data: 935 mužů (Wooldridge, wage2), mzda, vzdělání, IQ a rodinné zázemí
# =============================================================================

# ---- data ----
mzdy <- read.csv("data/wage2.csv")
# lwage = log měsíční mzdy, educ = roky vzdělání, IQ = výsledek IQ testu,
# sibs = počet sourozenců, brthord = pořadí narození mezi sourozenci,
# meduc, feduc = roky vzdělání matky a otce (NA = neznámé)
nrow(mzdy)


# =============================================================================
# Odhady: OLS a IV se čtyřmi kandidáty na instrument
# =============================================================================

# ---- modely ----
ols      <- lm(lwage ~ educ, data = mzdy)
iv_sibs  <- ivreg(lwage ~ educ | sibs, data = mzdy)
iv_brth  <- ivreg(lwage ~ educ | brthord, data = mzdy)
iv_meduc <- ivreg(lwage ~ educ | meduc, data = mzdy)
iv_feduc <- ivreg(lwage ~ educ | feduc, data = mzdy)

# ---- ols ----
# OLS: vzdělání je „nemocná“ proměnná (nadání je v náhodné složce)
coeftest(ols, vcov = vcovHC(ols, type = "HC1"))

# ---- iv-summary ----
# Výstup ivreg() vypadá skoro jako u lm(); za svislou čarou jsou instrumenty.
# Robustní se; diagnostiku přidá volba diagnostics = TRUE
summary(iv_sibs, vcov. = vcovHC, type = "HC1")


# =============================================================================
# IV jako podíl: první stupeň a redukovaný tvar
# =============================================================================

# ---- iv-sibs ----
# první stupeň: vliv instrumentu na nemocnou proměnnou
prvni <- lm(educ ~ sibs, data = mzdy)
# redukovaný tvar: vliv instrumentu na y
redukovany <- lm(lwage ~ sibs, data = mzdy)
# IV odhad = podíl obou sklonů (stejné číslo jako z ivreg())
coef(redukovany)[["sibs"]] / coef(prvni)[["sibs"]]
coef(iv_sibs)[["educ"]]

# ---- iv-vsechny ----
# Srovnání čtyř kandidátů. ivreg() vynechá řádky, kde instrument je NA.
modely <- list(sibs = iv_sibs, brthord = iv_brth,
               meduc = iv_meduc, feduc = iv_feduc)
srovnani <- t(sapply(names(modely), function(z) {
  m  <- modely[[z]]
  d  <- mzdy[!is.na(mzdy[[z]]), ]
  s  <- summary(m, vcov. = vcovHC, type = "HC1", diagnostics = TRUE)
  c(n        = nobs(m),
    b_iv     = coef(m)[["educ"]],
    se_iv    = s$coefficients["educ", "Std. Error"],
    # síla prvního stupně: robustní F (test slabých instrumentů)
    F_prvni  = s$diagnostics["Weak instruments", "statistic"],
    r_z_educ = cor(d[[z]], d$educ),
    r_z_IQ   = cor(d[[z]], d$IQ))
}))
round(srovnani, 3)


# =============================================================================
# Je instrument zdravý? IQ jako viditelná část náhodné složky
# =============================================================================

# ---- iq-pomer ----
# OLS: vychýlení úměrné corr(x, u); IV: úměrné corr(z, u) / corr(z, x).
# IQ je nedokonalé měřítko nadání (část u), ale ukazuje, kterým směrem se věci mají.
c(OLS = cor(mzdy$educ, mzdy$IQ),
  round(abs(srovnani[, "r_z_IQ"] / srovnani[, "r_z_educ"]), 3))


# =============================================================================
# Slabý instrument: s IQ mezi kontrolními proměnnými
# =============================================================================

# ---- iv-iq ----
# Kontrolní proměnná (IQ) se uvádí před svislou čarou i za ní: je zdravá,
# a proto je sama sobě instrumentem.
iv_sibs_iq <- ivreg(lwage ~ educ + IQ | sibs + IQ, data = mzdy)
# Weak instruments: sourozenci teď vysvětlují vzdělání mnohem méně
summary_iq <- summary(iv_sibs_iq, vcov. = vcovHC, type = "HC1",
                      diagnostics = TRUE)
summary_iq
ols_iq <- lm(lwage ~ educ + IQ, data = mzdy)
coeftest(ols_iq, vcov = vcovHC(ols_iq, type = "HC1"))


# =============================================================================
# 2SLS „ručně“ ve dvou krocích (jen pro pochopení, v praxi ivreg())
# =============================================================================

# ---- rucne ----
# 1. stupeň: zdravá část vzdělání = vyrovnané hodnoty
mzdy$educ_hat <- fitted(lm(educ ~ sibs, data = mzdy))
# 2. stupeň: mzdu vysvětlujeme zdravou částí vzdělání
rucne <- lm(lwage ~ educ_hat, data = mzdy)
# stejný odhad, jiná se (klasické, pro srovnání) a jiný R^2
rbind(rucne = c(coef(summary(rucne))["educ_hat", 1:2],
                R2 = summary(rucne)$r.squared),
      ivreg = c(coef(summary(iv_sibs))["educ", 1:2],
                R2 = summary(iv_sibs)$r.squared))

# ---- rezidua ----
# Rozdíl je v reziduích: ivreg() počítá y - b0 - b1*educ, ruční 2. stupeň
# y - b0 - b1*educ_hat
b <- coef(iv_sibs)
c(rucne = sd(residuals(rucne)),
  ivreg = sd(mzdy$lwage - b[1] - b[2] * mzdy$educ))


# =============================================================================
# Test přeidentifikujících omezení (Sarganův): dva instrumenty, jedna nemocná
# =============================================================================

# ---- sargan ----
rodice <- mzdy[complete.cases(mzdy[, c("meduc", "feduc")]), ]
# každý instrument zvlášť
sapply(c("meduc", "feduc"), function(z)
  coef(ivreg(as.formula(paste("lwage ~ educ |", z)), data = rodice))[["educ"]])
# oba najednou (2SLS) a diagnostika: Sargan = shodují se instrumenty?
iv_rodice <- ivreg(lwage ~ educ | meduc + feduc, data = rodice)
summary(iv_rodice, vcov = sandwich, diagnostics = TRUE)$diagnostics

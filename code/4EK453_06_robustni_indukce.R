# =============================================================================
# 4EK453 Principy ekonometrie – Přednáška 6: Robustní statistická indukce
# Jan Zouhar, Katedra ekonometrie FIS VŠE v Praze
#
# Skript reprodukuje všechny výsledky z prezentace. Spouštějte ho z kořene
# složky kurzu (data se čtou z podsložky data/).
# Řádky „# ---- název ----“ označují bloky, které se zobrazují v prezentaci.
#
# Potřebné balíčky: sandwich, lmtest, car
#   install.packages(c("sandwich", "lmtest", "car"))
# =============================================================================

# ---- nastaveni ----
options(scipen = 10)  # raději víc číslic než vědecký zápis (0.109 místo 1.090e-01)
library(sandwich)     # sendvičové odhady varianční matice: vcovHC(), vcovCL(), NeweyWest()
library(lmtest)       # coeftest(), bptest(), lrtest()
library(car)          # linearHypothesis(), deltaMethod()
set.seed(453)         # bootstrap a simulace jsou náhodné; seed zajistí stejné výsledky
se <- function(V) sqrt(diag(V))  # směrodatné chyby = odmocniny z diagonály varianční matice


# =============================================================================
# Ojeté škodovky: model z přednášky 5
# =============================================================================

# ---- skoda-model ----
sk <- read.csv("data/skoda.csv", sep = ";", fileEncoding = "CP1250")
sk$palivo <- factor(sk$palivo)  # faktor: úrovně zůstanou i ve výběrech bez LPG (8 aut)
m <- lm(log(cena) ~ I(km/1000) + rok + model + kombi + palivo, data = sk)

# ---- skoda-se ----
# Tři odhady směrodatných chyb pro tentýž model
se_klasicke <- se(vcov(m))                  # z sigma^2 (X'X)^{-1}, předpokládá homoskedasticitu
se_HC3      <- se(vcovHC(m, type = "HC3"))  # sendvič
round(cbind(odhad = coef(m), se_klasicke, se_HC3), 5)

# ---- skoda-rezidua ----
# Heteroskedasticita: rezidua vs. vyrovnané hodnoty; Breuschův-Paganův test
plot(fitted(m), resid(m), xlab = "vyrovnané hodnoty", ylab = "rezidua")
bptest(m)

# ---- skoda-HC ----
# Varianty sendviče: HC0 (White), HC1, HC3
V <- list(klasicke = vcov(m),
          HC0 = vcovHC(m, type = "HC0"),
          HC1 = vcovHC(m, type = "HC1"),
          HC3 = vcovHC(m, type = "HC3"))
round(sapply(V, se), 5)

# ---- skoda-coeftest ----
# Výpis koeficientů s robustními směrodatnými chybami, t-statistikami a p-hodnotami
coeftest(m, vcov = vcovHC(m, type = "HC3"))


# ---- mc-uvod ----
# =============================================================================
# Monte Carlo: pokrytí 95% intervalů spolehlivosti při heteroskedasticitě
# =============================================================================
# Populační model: y = 1 + 1*x + u, x ~ Exp(1), u = x * e, e ~ N(0, 1),
# tj. var(u | x) = x^2 (rozptyl roste s x). Zajímá nás, v kolika procentech
# výběrů obsahuje 95% interval spolehlivosti skutečnou hodnotu beta_1 = 1.

# ---- mc-pokryti ----
pokryti <- function(n, R = 1000) {
  vysledky <- replicate(R, {
    x <- rexp(n)
    y <- 1 + 1*x + x * rnorm(n)
    mm <- lm(y ~ x)
    c_krit <- qt(0.975, df = n - 2)
    chyby <- c(klasicke = sqrt(vcov(mm)[2, 2]),
               HC0 = sqrt(vcovHC(mm, type = "HC0")[2, 2]),
               HC3 = sqrt(vcovHC(mm, type = "HC3")[2, 2]))
    abs(coef(mm)[2] - 1) <= c_krit * chyby  # TRUE = interval pokrývá beta_1
  })
  rowMeans(vysledky)
}
pokryti_tab <- sapply(c(n_30 = 30, n_100 = 100, n_1000 = 1000), pokryti)
pokryti_tab  # podíl výběrů, ve kterých 95% interval pokryl skutečnou hodnotu


# ---- wages-uvod ----
# =============================================================================
# Mzdy v panelu: 595 lidí, každý pozorován 7 let (4,165 řádků)
# =============================================================================
# lwage … log mzdy, exper … roky praxe, wks … odpracované týdny v roce,
# ed … roky vzdělání, union, sex, black, married, smsa (bydlí ve městě),
# id … identifikátor člověka, time … rok (1–7)

# ---- wages-model ----
w <- read.csv("data/Wages.csv")
mw <- lm(lwage ~ exper + I(exper^2) + wks + ed + union + sex + black + married + smsa,
         data = w)

# ---- wages-cluster ----
# shluk = člověk (id)
V_CL <- vcovCL(mw, cluster = ~ id, type = "HC1")
round(cbind(odhad = coef(mw), klasicke = se(vcov(mw)),
            HC1 = se(vcovHC(mw, type = "HC1")),
            shlukove = se(V_CL)), 5)

# ---- wages-coeftest ----
coeftest(mw, vcov = V_CL)


# ---- okun-uvod ----
# =============================================================================
# Časové řady: Okunův zákon (čtvrtletní data)
# =============================================================================
# unemp … míra nezaměstnanosti (%), growth … čtvrtletní růst HDP (%)

# ---- okun-hac ----
ok <- read.csv("data/okun2.csv", sep = ";")
# změna nezaměstnanosti
ok$d_unemp <- c(NA, diff(ok$unemp))
mo <- lm(d_unemp ~ growth, data = ok)
V_NW <- NeweyWest(mo, lag = 4, prewhite = FALSE)
round(cbind(odhad = coef(mo), klasicke = se(vcov(mo)),
            HC1 = se(vcovHC(mo, type = "HC1")),
            NeweyWest = se(V_NW)), 5)


# ---- bootstrap-uvod ----
# =============================================================================
# Bootstrap
# =============================================================================

# ---- bootstrap-smycka ----
# Párový bootstrap: losujeme s vracením celá pozorování
n <- nrow(sk)
B <- 2000  # počet replikací
bs <- matrix(NA, nrow = B, ncol = length(coef(m)),
             dimnames = list(NULL, names(coef(m))))
for (b in 1:B) {
  vyber <- sample(n, n, replace = TRUE)  # s vracením
  bs[b, ] <- coef(lm(formula(m), data = sk[vyber, ]))
}
sum(is.na(bs[, "palivoLPG"]))  # replikace bez LPG
# sd napříč replikacemi
se_boot <- apply(bs, 2, sd, na.rm = TRUE)
round(cbind(se_klasicke, se_HC3, se_boot), 5)

# ---- bootstrap-hist ----
hist(bs[, "rok"], breaks = 50, main = "Bootstrapové odhady koeficientu u rok",
     xlab = "koeficient u rok")

# ---- bootstrap-shluky ----
# Panel: losujeme celé lidi (všech 7 let jednoho člověka najednou)
radky_cloveka <- split(seq_len(nrow(w)), w$id)   # pro každé id čísla jeho řádků
N_lidi <- length(radky_cloveka)
B_w <- 1000
bs_poz <- bs_shl <- matrix(NA, nrow = B_w, ncol = length(coef(mw)),
                           dimnames = list(NULL, names(coef(mw))))
for (b in 1:B_w) {
  # (a) chybně: losujeme jednotlivé řádky
  vyber_poz <- sample(nrow(w), nrow(w), replace = TRUE)
  bs_poz[b, ] <- coef(lm(formula(mw), data = w[vyber_poz, ]))
  # (b) správně: losujeme lidi a bereme všechny jejich řádky
  lide <- sample(N_lidi, N_lidi, replace = TRUE)
  vyber_shl <- unlist(radky_cloveka[lide])
  bs_shl[b, ] <- coef(lm(formula(mw), data = w[vyber_shl, ]))
}
round(cbind(shlukove    = se(V_CL),
            boot_radky  = apply(bs_poz, 2, sd),
            boot_shluky = apply(bs_shl, 2, sd)), 5)

# ---- vrchol ----
# Vrchol mzdového profilu: -b1 / (2 b2)
vrchol <- function(b)
  -b[["exper"]] / (2 * b[["I(exper^2)"]])
vrchol(coef(mw))  # bodový odhad
# vrchol v každé replikaci (bootstrap po lidech)
vrcholy_bs <- apply(bs_shl, 1, vrchol)
sd(vrcholy_bs)    # bootstrapová se
quantile(vrcholy_bs, c(0.025, 0.975))  # 95% CI

# ---- vrchol-delta ----
# Pro srovnání: delta metoda (lineární aproximace funkce odhadů)
deltaMethod(mw, "-b1 / (2*b2)", parameterNames = paste0("b", 0:9), vcov. = V_CL)


# ---- jackknife-uvod ----
# =============================================================================
# Jackknife
# =============================================================================

# ---- jackknife ----
# Vynecháme postupně každé pozorování a model odhadneme znovu
jk <- t(sapply(1:n, function(i) coef(lm(formula(m), data = sk[-i, ]))))
jk_prumer <- colMeans(jk)
V_jk <- (n - 1) / n * crossprod(sweep(jk, 2, jk_prumer))  # (n-1)/n * součet čtverců odchylek
round(cbind(se_HC3, se_jackknife = se(V_jk)), 5)


# ---- fgls-uvod ----
# =============================================================================
# Vážené nejmenší čtverce (FGLS)
# =============================================================================

# ---- fgls ----
# 1. model rozptylu: log(û^2) vysvětlujeme pomocí regresorů
sk$log_u2 <- log(resid(m)^2)
pomocna <- lm(update(formula(m), log_u2 ~ .), data = sk)
vahy <- 1 / exp(fitted(pomocna))  # 1 / odhad rozptylu
# 2. vážené nejmenší čtverce + robustní se
m_wls <- lm(formula(m), data = sk, weights = vahy)
round(cbind(OLS = coef(m), se_OLS = se_HC3,
            WLS = coef(m_wls),
            se_WLS = se(vcovHC(m_wls, type = "HC3"))), 5)


# ---- testy-uvod ----
# =============================================================================
# Porovnávání modelů: záleží cena na palivu?
# =============================================================================

# ---- robust-F ----
H0 <- c("palivodiesel = 0", "palivoLPG = 0")
linearHypothesis(m, H0)              # klasický F-test
linearHypothesis(m, H0,              # robustní (Waldův) test
                 vcov. = vcovHC(m, type = "HC3"))

# ---- wlrlm ----
m0 <- update(m, . ~ . - palivo)  # bez paliva
SSR_U <- deviance(m)    # neomezený model
SSR_R <- deviance(m0)   # omezený model
c(W  = n * (SSR_R - SSR_U) / SSR_U,
  LR = n * log(SSR_R / SSR_U),
  LM = n * (SSR_R - SSR_U) / SSR_R,
  mF = 2 * linearHypothesis(m, H0)$F[2],  # m = 2
  krit = qchisq(0.95, df = 2))

# ---- lm-pomocna ----
# LM = n * R^2 z pomocné regrese: rezidua omezeného
# modelu vysvětlujeme pomocí všech regresorů
pom <- lm(update(formula(m), resid(m0) ~ .), data = sk)
n * summary(pom)$r.squared

# ---- lrtest ----
lrtest(m0, m)   # LR test přímo z log-věrohodností obou modelů


# ---- trojice-hladina ----
# Skutečná hladina W, LR a LM (s kritickou hodnotou chi^2) a F-testu, když H0 platí:
# 4 regresory, testujeme 2 restrikce (beta_3 = beta_4 = 0), nominální hladina 5 %
set.seed(4535)
hladina <- function(n, R = 5000) {
  zamitnuti <- replicate(R, {
    X <- cbind(1, matrix(rnorm(n * 4), n, 4))
    y <- drop(X[, 1:3] %*% c(1, 1, 1)) + rnorm(n)  # H0 platí
    SSR_U <- sum(lm.fit(X, y)$residuals^2)
    SSR_R <- sum(lm.fit(X[, 1:3], y)$residuals^2)
    c(W  = n * (SSR_R - SSR_U) / SSR_U > qchisq(0.95, 2),
      LR = n * log(SSR_R / SSR_U) > qchisq(0.95, 2),
      LM = n * (SSR_R - SSR_U) / SSR_R > qchisq(0.95, 2),
      F  = (SSR_R - SSR_U) / 2 / (SSR_U / (n - 5)) > qf(0.95, 2, n - 5))
  })
  rowMeans(zamitnuti)  # podíl zamítnutí platné H0
}
hladiny <- sapply(c(15, 30, 100, 1000), hladina)
colnames(hladiny) <- paste0("n=", c(15, 30, 100, 1000))
round(hladiny, 3)

# ---- lrm3-uvod ----
# =============================================================================
# Doplněk: LRM.3 a numerická přesnost (rok a jeho mocniny)
# =============================================================================
# Matice X s regresory year, year^2, year^3 má matematicky plnou hodnost,
# počítač ji ale kvůli konečné přesnosti za regulární nepovažuje.

# ---- lrm3-numerika ----
set.seed(2005)
roky <- data.frame(year = rep(1990:2020, each = 20))
roky$t <- roky$year - 2005  # vycentrovaný rok
roky$y <- 0.02*roky$t - 0.00003*roky$t^3 +
  rnorm(nrow(roky), sd = 0.3)
coef(lm(y ~ year + I(year^2) + I(year^3), data = roky))
coef(lm(y ~ t + I(t^2) + I(t^3), data = roky))

# ---- lrm3-kappa ----
# kappa: kolikrát se může zvětšit relativní
# chyba; při kappa kolem 10^d ztratíme zhruba
# d ze 16 platných číslic
X2 <- model.matrix(~ year + I(year^2), roky)
kappa(X2, exact = TRUE)  # 10^11: lm() zvládne
# vzorec (X'X)^(-1): kappa(X'X) = kappa(X)^2
chyba <- tryCatch(solve(crossprod(X2)),
                  error = conditionMessage)
cat(strwrap(chyba, 40), sep = "\n")

# ---- lrm3-kappa-kubicky ----
# Kubický model: kappa před a po vycentrování roku
kappa(model.matrix(~ year + I(year^2) + I(year^3), data = roky), exact = TRUE)
kappa(model.matrix(~ t + I(t^2) + I(t^3), data = roky), exact = TRUE)

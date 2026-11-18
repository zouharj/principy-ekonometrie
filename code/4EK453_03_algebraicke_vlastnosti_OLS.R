# =============================================================================
# 4EK453 Principy ekonometrie – Přednáška 3: Algebraické vlastnosti OLS
# Jan Zouhar, Katedra ekonometrie FIS VŠE v Praze
#
# Skript reprodukuje všechny výsledky z prezentace. Spouštějte ho z kořene
# složky kurzu (data se čtou z podsložky data/).
# Řádky „# ---- název ----“ označují bloky, které se zobrazují v prezentaci.
# =============================================================================

# ---- nastaveni ----
options(scipen = 10)  # raději víc číslic než vědecký zápis (0.109 místo 1.090e-01)


# =============================================================================
# Námět 1: změny jednotek (ojeté škodovky)
# =============================================================================

# ---- n1-data ----
sk <- read.csv("data/skoda.csv", sep = ";", fileEncoding = "CP1250")
head(sk)

# ---- n1-m1 ----
# (1) nájezd v km, rok výroby
m1 <- lm(log(cena) ~ km + rok, data = sk)
summary(m1)

# ---- n1-m2 ----
# (2) nájezd v desítkách tisíc km
m2 <- lm(log(cena) ~ I(km/10000) + rok, data = sk)
summary(m2)

# ---- n1-m3 ----
sk$stari <- 2004 - sk$rok  # stáří vozu (data jsou z roku 2004)
m3 <- lm(log(cena) ~ km + stari, data = sk)
summary(m3)

# ---- n1-m4 ----
sk$x2 <- 2*sk$km - 3*sk$rok  # syntetická proměnná: lineární kombinace km a rok
m4 <- lm(log(cena) ~ km + x2, data = sk)
summary(m4)

# ---- n1-rozdily ----
# Vyrovnané hodnoty a rezidua jsou ve všech modelech stejné (až na zaokrouhlení)
max(abs(fitted(m2) - fitted(m1)))
max(abs(fitted(m3) - fitted(m1)))
max(abs(fitted(m4) - fitted(m1)))
max(abs(resid(m4) - resid(m1)))

# ---- n1-matice-A ----
# Řešení: regresory modelu (4) jsou XA pro regulární A, koeficienty jsou A^{-1} beta
X <- model.matrix(m1)
A <- rbind(c(1, 0,  0),
           c(0, 1,  2),
           c(0, 0, -3))
all.equal(X %*% A, model.matrix(m4), check.attributes = FALSE)  # TRUE
solve(A) %*% coef(m1)  # = coef(m4)
coef(m4)


# =============================================================================
# Námět 2: změny referenční kategorie
# =============================================================================

# ---- n2-fel ----
sk$model <- relevel(factor(sk$model), ref = "Felicia")
r_fel <- lm(log(cena) ~ I(km/10000) + rok + model, data = sk)
summary(r_fel)

# ---- n2-oct ----
sk$model <- relevel(factor(sk$model), ref = "Octavia")
r_oct <- lm(log(cena) ~ I(km/10000) + rok + model, data = sk)
summary(r_oct)

# ---- n2-sup ----
sk$model <- relevel(factor(sk$model), ref = "Superb")
r_sup <- lm(log(cena) ~ I(km/10000) + rok + model, data = sk)
summary(r_sup)

# ---- n2-vztahy ----
# Koeficienty z modelu s referencí Octavia dopočítané z modelu s referencí Felicia
b <- coef(r_fel)
c(intercept      = b[["(Intercept)"]] + b[["modelOctavia"]],
  modelFelicia   = -b[["modelOctavia"]],
  modelSuperb    = b[["modelSuperb"]] - b[["modelOctavia"]])
coef(r_oct)


# =============================================================================
# Námět 3a: trendy (HDP USA a lososi z Norska)
# =============================================================================

# ---- n3-data ----
sa <- read.csv("data/salmon.csv")
names(sa)[1] <- "year"
naivni <- lm(gdp ~ salmon, data = sa)  # zdánlivá regrese
summary(naivni)

# ---- n3-modely ----
# Řady očištěné od trendu (přípona _det, detrended) = rezidua z regrese,
# která řadu vysvětluje pomocí interceptu a year
sa$gdp_det    <- resid(lm(gdp ~ year, data = sa))
sa$salmon_det <- resid(lm(salmon ~ year, data = sa))

t1 <- lm(gdp ~ year + salmon, data = sa)          # (1) trend v regresi
t2 <- lm(gdp_det ~ 0 + salmon_det, data = sa)     # (2) očištěné obě strany
t3 <- lm(gdp_det ~ year + salmon, data = sa)      # (3) očištěné y, trend v regresi
t4 <- lm(gdp ~ 0 + salmon_det, data = sa)         # (4) očištěné jen x
t5 <- lm(gdp_det ~ salmon, data = sa)             # (5) očištěné jen y, bez trendu

# ---- n3-srovnani ----
sapply(list(t1, t2, t3, t4, t5), function(m) tail(coef(m), 1))  # koeficient u salmon
sapply(list(t1, t2, t3, t4, t5), deviance)                       # SSR
sapply(list(t1, t2, t3, t4, t5), function(m) summary(m)$r.squared)  # R^2
# Pozor: u modelů bez interceptu (2) a (4) vypisuje R necentrovaný R^2.


# =============================================================================
# Námět 3b: sezónnost (prodej ledniček v USA, čtvrtletně 1978–1985)
# =============================================================================
# frig  … prodej ledniček
# dur   … výdaje na zboží dlouhodobé spotřeby
# houst … zahájené stavby bytových jednotek v USA za čtvrtletí (tis.),
#         FRED HOUSTNSA (sezónně neočištěno)

# ---- n3s-data ----
dg <- read.csv("data/durgoods2.csv")
dg$ctvrtleti <- factor(substr(dg$obs, 6, 6))  # 1, 2, 3, 4

# ---- n3s-modely ----
s0 <- lm(frig ~ dur + houst, data = dg)  # naivní model bez sezónnosti a trendu

# Sezónně očištěné a odtrendované řady (přípona _sa, seasonally adjusted)
# = rezidua z regrese, která řadu vysvětluje pomocí čtvrtletních dummy a trendu
ocisti <- function(v) resid(lm(v ~ ctvrtleti + time, data = dg))
dg$frig_sa  <- ocisti(dg$frig)
dg$dur_sa   <- ocisti(dg$dur)
dg$houst_sa <- ocisti(dg$houst)

s1 <- lm(frig ~ ctvrtleti + time + dur + houst, data = dg)     # (1) dummy a trend v regresi
s2 <- lm(frig_sa ~ 0 + dur_sa + houst_sa, data = dg)           # (2) očištěné všechny řady
s3 <- lm(frig_sa ~ ctvrtleti + time + dur + houst, data = dg)  # (3) očištěné y, dummy a trend v regresi
s4 <- lm(frig ~ 0 + dur_sa + houst_sa, data = dg)              # (4) očištěné jen x
s5 <- lm(frig_sa ~ dur + houst, data = dg)                     # (5) očištěné jen y, bez dummy a trendu

# ---- n3s-srovnani ----
modely_s <- list(s0, s1, s2, s3, s4, s5)
sapply(modely_s, function(m) tail(coef(m), 2))                # koeficienty u dur a houst
sapply(modely_s, deviance)                                    # SSR
sapply(modely_s, function(m) summary(m)$r.squared)            # R^2 (u (2) a (4) necentrovaný)


# =============================================================================
# Námět 4: diagonála (X'X)^{-1} a 1 / [SST_j (1 - R_j^2)]
# =============================================================================

# ---- n4 ----
X <- model.matrix(r_fel)  # regresory (reference Felicia)
diag_inv <- diag(solve(t(X) %*% X))
vzorec <- sapply(2:ncol(X), function(j) {
  xj   <- X[, j]
  SSTj <- sum((xj - mean(xj))^2)
  # R2j: xj vysvětlujeme pomocí ostatních regresorů
  R2j  <- summary(lm(xj ~ X[, -c(1, j)]))$r.squared
  1 / (SSTj * (1 - R2j))
})
cbind(diag_inv = diag_inv[-1], vzorec)

# ---- r2 ----
# =============================================================================
# Další: centrovaný a necentrovaný R^2
# =============================================================================

y <- log(sk$cena); yhat <- fitted(m1); uhat <- resid(m1)
SST <- sum((y - mean(y))^2); SSE <- sum((yhat - mean(y))^2); SSR <- sum(uhat^2)
c(SST = SST, SSE_plus_SSR = SSE + SSR)
c(R2 = SSE / SST, R2_u = sum(yhat^2) / sum(y^2))

# ---- r2-parcialni ----
# Parciální R^2 proměnné x_j (FWL): y i x_j očistíme od ostatních regresorů
X <- model.matrix(r_fel); y <- log(sk$cena)
j_vse <- setNames(2:ncol(X), colnames(X)[-1])
parc <- sapply(j_vse, function(j) {
  y_o <- resid(lm(y ~ X[, -c(1, j)]))      # y očištěné od ostatních
  x_o <- resid(lm(X[, j] ~ X[, -c(1, j)])) # x_j očištěné od ostatních
  summary(lm(y_o ~ x_o))$r.squared
})
# přírůstek R^2, když x_j přidáme do modelu bez ní
prir <- sapply(j_vse, function(j)
  summary(r_fel)$r.squared - summary(lm(y ~ X[, -c(1, j)]))$r.squared)
rbind(parcialni = parc, prirustek = prir)
c(soucet_parc = sum(parc), soucet_prir = sum(prir), R2 = summary(r_fel)$r.squared)

# =============================================================================
# 4EK453 Principy ekonometrie – Přednáška 7: Výběr modelu a specifikace
# Jan Zouhar, Katedra ekonometrie FIS VŠE v Praze
#
# Skript reprodukuje všechny výsledky z prezentace. Spouštějte ho z kořene
# složky kurzu (data se čtou z podsložky data/).
# Řádky „# ---- název ----“ označují bloky, které se zobrazují v prezentaci.
#
# Potřebné balíčky: sandwich, lmtest, glmnet, hdm (splines je součástí R)
#   install.packages(c("sandwich", "lmtest", "glmnet", "hdm"))
# Doplňkový materiál o splinech, GAM a LOESS: 4EK453_07_doplnek_spliny.R
# =============================================================================

# ---- nastaveni ----
options(scipen = 10)  # raději víc číslic než vědecký zápis
library(splines)      # spliny: bs(), ns()
library(sandwich)     # shlukové směrodatné chyby: vcovCL()
library(lmtest)       # coeftest(), resettest()
library(glmnet)       # ridge a lasso: glmnet(), cv.glmnet()
library(hdm)          # dvojí výběr pomocí lassa: rlassoEffect()
set.seed(453)         # simulace a křížová validace jsou náhodné
se <- function(V) sqrt(diag(V))


# =============================================================================
# Námět 1: kterou kontrolní proměnnou vybere AIC?
# =============================================================================

# ---- n1-sim ----
# Populační model: y = x + 0.5 z1 + 0.5 z2 + u, skutečný efekt x je 1.
# z1 je zavádějící proměnná (ovlivňuje x i y), z2 ovlivňuje jen y.
n  <- 2000
z1 <- rnorm(n); z2 <- rnorm(n)
x  <- 0.6 * z1 + rnorm(n)
y  <- x + 0.5 * z1 + 0.5 * z2 + rnorm(n)
n1 <- list("bez kontrolních proměnných" = lm(y ~ x), "+ z1" = lm(y ~ x + z1), "+ z2" = lm(y ~ x + z2))
sapply(n1, function(m) c(beta_x = coef(m)[["x"]], AIC = AIC(m)))


# =============================================================================
# Námět 2: lasso vyřadí zavádějící proměnnou (Monte Carlo)
# =============================================================================

# ---- n2-lasso ----
# 50 kandidátů na kontrolní proměnné z1, ..., z50; skutečný efekt x je 0.5.
# z1 silně ovlivňuje x a slabě y (zavádějící), z2 ovlivňuje jen y, ostatní nic.
# Naivní postup: lasso vysvětluje y pomocí x (bez penalizace) a kontrolních proměnných, lambda
# křížovou validací, pak OLS s vybranými kontrolními proměnnými.
# Dvojí výběr: rlassoEffect(), viz snímek o dvojím výběru.
jeden_vyber <- function(n = 200, p = 50) {
  Z <- matrix(rnorm(n * p), n)
  x <- 0.9 * Z[, 1] + 0.45 * rnorm(n)
  y <- 0.5 * x + 0.25 * Z[, 1] + Z[, 2] + rnorm(n)
  cv  <- cv.glmnet(cbind(x, Z), y, penalty.factor = c(0, rep(1, p)))
  vybrane <- which(as.numeric(coef(cv, s = "lambda.1se"))[-(1:2)] != 0)
  naivni <- summary(lm(y ~ x + Z[, vybrane]))$coefficients["x", 1:2]
  dvoji  <- rlassoEffect(Z, y, x, method = "double selection")
  c(naivni = naivni[[1]], naivni_se = naivni[[2]], dvoji = dvoji$alpha[[1]],
    dvoji_se = dvoji$se[[1]], z1_vybrana = 1 %in% vybrane)
}
R <- 300
lasso_mc <- t(replicate(R, jeden_vyber()))
pokryti <- function(b, s) mean(abs(b - 0.5) < 1.96 * s)
c(prumer_naivni = mean(lasso_mc[, "naivni"]), prumer_dvoji = mean(lasso_mc[, "dvoji"]),
  sd_naivni = sd(lasso_mc[, "naivni"]), sd_dvoji = sd(lasso_mc[, "dvoji"]),
  pokryti_naivni = pokryti(lasso_mc[, "naivni"], lasso_mc[, "naivni_se"]),
  pokryti_dvoji = pokryti(lasso_mc[, "dvoji"], lasso_mc[, "dvoji_se"]),
  z1_vybrana = mean(lasso_mc[, "z1_vybrana"]))


# =============================================================================
# Námět 3: test po vyřazení nevýznamné kontrolní proměnné (Monte Carlo)
# =============================================================================

# ---- n3-pretest ----
# y = x + gama * z + u, corr(x, z) = 0.7, n = 100. Testujeme H0: beta = 1 (platí).
# Postup „předběžného testu“: když je z nevýznamná (|t| < 1.96), vyřadíme ji
# a beta testujeme v modelu bez z. Jak často zamítneme platnou H0?
predbezny_test <- function(gama, n = 100, R = 4000) {
  zamitnuti <- replicate(R, {
    z <- rnorm(n)
    x <- 0.7 * z + sqrt(1 - 0.7^2) * rnorm(n)
    y <- x + gama * z + rnorm(n)
    s_plny <- summary(lm(y ~ x + z))$coefficients
    s_test <- if (abs(s_plny["z", 3]) > 1.96) s_plny else summary(lm(y ~ x))$coefficients
    t_plny <- (s_plny["x", 1] - 1) / s_plny["x", 2]
    t_test <- (s_test["x", 1] - 1) / s_test["x", 2]
    c(s_z = abs(t_plny) > 1.96, po_testu = abs(t_test) > 1.96)
  })
  rowMeans(zamitnuti)
}
gamy <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7)
hladiny <- sapply(gamy, predbezny_test)
colnames(hladiny) <- gamy
round(100 * hladiny, 1)  # skutečná hladina testu v %, nominálně 5 %


# =============================================================================
# Obce: efekt místního referenda na volební účast
# =============================================================================

# ---- obce-data ----
# 6,245 obcí ve sněmovních volbách 2002, 2006 a 2010 (Dvořák, Zouhar
# a Novák, 2017). referendum = 1, pokud se v obci v daném volebním období konalo
# místní referendum; ucast v %. Směrodatné chyby shlukujeme podle obce.
obce <- read.csv("data/obce_v_CR.csv", fileEncoding = "CP1250")
obce$rok <- factor(obce$rok)
se_obce <- function(m) se(vcovCL(m, cluster = ~ ID))[["referendum"]]

# ---- obce-referendum ----
o1 <- lm(ucast ~ referendum + rok, data = obce)
o2 <- lm(ucast ~ referendum + rok + log_velikost, data = obce)
o3 <- lm(ucast ~ referendum + rok + log_velikost + nezamestnanost + vek_0_44 + vek_75plus,
         data = obce)
sapply(list(o1, o2, o3), function(m) c(beta = coef(m)[["referendum"]], se = se_obce(m)))

# ---- obce-decily ----
# Velikost obce souvisí s referendem i s účastí: zavádějící proměnná
decil <- cut(obce$log_velikost, quantile(obce$log_velikost, 0:10 / 10), include.lowest = TRUE)
round(cbind(referendum_v_procentech = 100 * tapply(obce$referendum, decil, mean),
            ucast = tapply(obce$ucast, decil, mean)), 1)


# =============================================================================
# Vychýlení vs. rozptyl (Monte Carlo)
# =============================================================================

# ---- mse-n ----
# y = x + z + u, x = 0.5 z + e: vynechání z vychýlí odhad efektu x.
# MSE odhadu beta (= vychýlení^2 + rozptyl) podle rozsahu výběru
mse_n <- function(n, R = 2000) {
  b <- replicate(R, {
    z <- rnorm(n); x <- 0.5 * z + rnorm(n); y <- x + z + rnorm(n)
    c(bez_z = coef(lm(y ~ x))[["x"]], se_z = coef(lm(y ~ x + z))[["x"]])
  })
  rowMeans((b - 1)^2)
}
rozsahy <- c(5, 10, 20, 50, 100, 200, 500)
mse <- sapply(rozsahy, mse_n)
colnames(mse) <- rozsahy
round(mse, 4)

# ---- terc ----
# Tři modely pro efekt x (skutečně 1), n = 100:
#   málo kontrolních proměnných:  vynechá zavádějící z (vychýlení),
#   správně:       zahrne z,
#   hodně kontrolních proměnných: zahrne z a navíc pět proměnných w, které ovlivňují jen x
#                  (žádné vychýlení neodstraní, jen zvětší rozptyl).
terc <- t(replicate(2000, {
  n <- 100
  z <- rnorm(n); W <- matrix(rnorm(5 * n), n)
  x <- 0.5 * z + W %*% rep(0.5, 5) + 0.5 * rnorm(n)
  y <- x + z + rnorm(n)
  c(malo = coef(lm(y ~ x))[["x"]], spravne = coef(lm(y ~ x + z))[["x"]],
    hodne = coef(lm(y ~ x + z + W))[["x"]])
}))
round(rbind(prumer = colMeans(terc), sd = apply(terc, 2, sd)), 3)


# =============================================================================
# Prediktivní kritéria: AIC, BIC, křížová validace
# =============================================================================

# ---- obce-tvary ----
# Tvar velikosti obce jako kontrolní proměnné: lineárně, kvadraticky,
# diskretizace do 20 intervalů (kvantilových), přirozený spline
hranice_int <- quantile(obce$log_velikost, 0:20 / 20)
obce$intervaly <- cut(obce$log_velikost, hranice_int, include.lowest = TRUE)
tvary <- list(
  lineárně     = lm(ucast ~ referendum + rok + log_velikost + nezamestnanost + vek_0_44 + vek_75plus, data = obce),
  kvadraticky  = lm(ucast ~ referendum + rok + poly(log_velikost, 2) + nezamestnanost + vek_0_44 + vek_75plus, data = obce),
  "20 intervalů" = lm(ucast ~ referendum + rok + intervaly + nezamestnanost + vek_0_44 + vek_75plus, data = obce),
  spline       = lm(ucast ~ referendum + rok + ns(log_velikost, 6) + nezamestnanost + vek_0_44 + vek_75plus, data = obce)
)
round(sapply(tvary, function(m) c(beta = coef(m)[["referendum"]], se = se_obce(m),
                                  R2 = summary(m)$r.squared, AIC = AIC(m), BIC = BIC(m))), 3)

# ---- obce-cv ----
# 10násobná křížová validace; celé obce (všechny tři roky) jsou vždy ve stejné části,
# jinak by se model učil na téže obci, kterou pak predikuje
casti <- sample(rep(1:10, length.out = length(unique(obce$ID))))[match(obce$ID, unique(obce$ID))]
cv_mse <- function(vzorec) {
  chyby <- numeric(nrow(obce))
  for (k in 1:10) {
    m <- lm(vzorec, data = obce[casti != k, ])
    chyby[casti == k] <- obce$ucast[casti == k] - predict(m, obce[casti == k, ])
  }
  mean(chyby^2)
}
# LOOCV bez opakovaného odhadu: rezidua dělená (1 - páka)
loocv <- function(m) mean((resid(m) / (1 - hatvalues(m)))^2)
round(sapply(tvary, function(m) c(CV10 = cv_mse(formula(m)), LOOCV = loocv(m))), 2)

# ---- obce-reset ----
# RESET test: přidá do modelu druhou a třetí mocninu vyrovnaných hodnot
resettest(tvary$lineárně, power = 2:3)


# =============================================================================
# Regularizace: ridge a lasso
# =============================================================================

# ---- obce-cesty ----
# Účast vysvětlujeme pomocí několika proměnných a jejich čtverců;
# glmnet proměnné sám standardizuje (koeficienty vrací v původních jednotkách)
X_male <- model.matrix(~ referendum + rok + log_velikost + I(log_velikost^2) + nezamestnanost +
                         I(nezamestnanost^2) + vek_0_44 + vek_75plus, data = obce)[, -1]
# Cesty koeficientů pro graf: standardizované proměnné, ať jsou koeficienty srovnatelné
lasso <- glmnet(scale(X_male), obce$ucast, alpha = 1)  # alpha = 1 lasso
ridge <- glmnet(scale(X_male), obce$ucast, alpha = 0)  # alpha = 0 ridge
# Volba lambda 10násobnou křížovou validací (stejné části jako výše)
cv_lasso <- cv.glmnet(X_male, obce$ucast, alpha = 1, foldid = casti)
c(lambda_min = cv_lasso$lambda.min, lambda_1se = cv_lasso$lambda.1se)
round(cbind(OLS = coef(lm(obce$ucast ~ X_male)),
            lasso_1se = as.numeric(coef(cv_lasso, s = "lambda.1se"))), 4)

# ---- obce-dvoji-vyber ----
# Mnoho technických kontrolních proměnných:
# spliny spojitých proměnných, zvlášť po rocích
X_tech <- model.matrix(
  ~ rok * (ns(log_velikost, 8) +
             ns(nezamestnanost, 4) +
             ns(vek_0_44, 4) + ns(vek_75plus, 4)),
  data = obce)[, -1]
ncol(X_tech)
dvoji <- rlassoEffect(x = X_tech, y = obce$ucast,
                      d = obce$referendum,
                      method = "double selection")
summary(dvoji)

# ---- obce-lasso-naivni ----
# Pro srovnání jednoduchý výběr: lasso vybírá kontrolní proměnné jen podle predikce
# účasti, pak OLS. Velikost obce silně ovlivňuje účast, takže ji lasso ponechá.
vyber_y <- which(coef(rlasso(X_tech, obce$ucast))[-1] != 0)
length(vyber_y)
naivni_obce <- lm(obce$ucast ~ obce$referendum + X_tech[, vyber_y])
coef(summary(naivni_obce))[2, 1:2]

# ---- kvadrat-uvod ----
# =============================================================================
# Kvadratický člen: o kolik vzroste mzda mezi 10 a 11 lety praxe?
# =============================================================================

# ---- kvadrat-psid ----
# PSID, stejný model jako v přednášce 6: log mzdy s praxí a jejím čtvercem
w <- read.csv("data/Wages.csv")
mw <- lm(lwage ~ exper + I(exper^2) + wks + ed + union + sex + black + married + smsa,
         data = w)
b1 <- coef(mw)[["exper"]]
b2 <- coef(mw)[["I(exper^2)"]]
# sklon (mezní efekt) b1 + 2 b2 exper v bodech 10, 10.5 a 11
sklon <- b1 + 2 * b2 * c(10, 10.5, 11)
# rozdíl dosazených hodnot pro 11 a 10 let praxe: stejný jako sklon v 10.5
rozdil <- (b1 * 11 + b2 * 11^2) - (b1 * 10 + b2 * 10^2)
c(sklon, rozdil = rozdil)
100 * (exp(rozdil) - 1)  # přesně v procentech

# ---- spline-uvod ----
# =============================================================================
# Spliny: polynomy po částech s omezeními v uzlech (podrobně v doplňku)
# =============================================================================

# ---- spline-data ----
# Simulovaná data (stejná jako v doplňkovém materiálu), uzly 2.5, 5 a 7.5
set.seed(1)
xs <- sort(runif(500, 0, 10))
ys <- 10 + sin(xs) + rnorm(500, 0, 0.4) + sort(cumsum(rnorm(500, 0, 0.24))) - (1:500)/200
uzly <- c(2.5, 5, 7.5)
set.seed(453)

# ---- spline-omezeni ----
# Useknuté mocniny (x - k)_+^p; vypuštění mocniny p u všech uzlů = spojitost p-té derivace
usek <- function(x, k, p) if (p == 0) as.numeric(x > k) else pmax(x - k, 0)^p
baze <- function(x, vypustit) {
  B <- sapply(1:3, function(p) x^p)
  for (k in uzly) for (p in setdiff(0:3, vypustit)) B <- cbind(B, usek(x, k, p))
  B
}
mriz <- seq(0, 10, length.out = 401)
krivka <- function(f) {
  B <- f(xs)
  m <- lm(ys ~ B)
  list(y = drop(cbind(1, f(mriz)) %*% coef(m)), parametru = length(coef(m)))
}
omezeni <- list(
  "bez omezení"             = krivka(function(x) baze(x, integer(0))),
  "spojitost"               = krivka(function(x) baze(x, 0)),
  "+ spojitá 1. derivace"   = krivka(function(x) baze(x, 0:1)),
  "+ spojitá 2. derivace"   = krivka(function(x) baze(x, 0:2))
)
sapply(omezeni, `[[`, "parametru")

# ---- spline-stupne ----
stupne <- list(
  "lineární spline"           = krivka(function(x) bs(x, knots = uzly, degree = 1, Boundary.knots = c(0, 10))),
  "kvadratický spline"        = krivka(function(x) bs(x, knots = uzly, degree = 2, Boundary.knots = c(0, 10))),
  "kubický spline"            = krivka(function(x) bs(x, knots = uzly, degree = 3, Boundary.knots = c(0, 10))),
  "přirozený kubický spline"  = krivka(function(x) ns(x, knots = uzly, Boundary.knots = range(xs)))
)
sapply(stupne, `[[`, "parametru")

# ---- obce-polynomy ----
# Účast podle velikosti obce: polynomy stupně 1, 2 a 5 (graf v prezentaci
# ukazuje průměry ve 40 intervalech velikosti, aby nebylo třeba kreslit 18,735 bodů)
pol <- lapply(c(1, 2, 5), function(st) lm(ucast ~ poly(log_velikost, st), data = obce))
sapply(pol, function(m) summary(m)$r.squared)


# =============================================================================
# Dynamika: Okunův zákon s rozloženým zpožděním
# =============================================================================

# ---- okun-data ----
# Čtvrtletní data ČR: míra nezaměstnanosti (unemp, %) a růst HDP (growth, %)
ok <- read.csv("data/okun2.csv", sep = ";")
ok$d_unemp <- c(NA, diff(ok$unemp))  # změna nezaměstnanosti
zpozdi <- function(v, k) c(rep(NA, k), head(v, length(v) - k))
for (k in 1:4) ok[[paste0("growth_", k)]] <- zpozdi(ok$growth, k)

# ---- okun-zpozdeni ----
# Modely s 0 až 4 zpožděními růstu, všechny na stejném vzorku (kvůli AIC a BIC)
vzorek <- complete.cases(ok)
zpozdeni <- lapply(0:4, function(q) {
  vzorec <- reformulate(c("growth", if (q > 0) paste0("growth_", 1:q)), "d_unemp")
  lm(vzorec, data = ok[vzorek, ])
})
round(sapply(zpozdeni, function(m) c(okamzity = coef(m)[["growth"]],
                                     dlouhodoby = sum(coef(m)[-1]),
                                     AIC = AIC(m), BIC = BIC(m))), 3)

# ---- okun-adl ----
# ADL model: místo dalších zpoždění růstu zpožděná změna nezaměstnanosti;
# dlouhodobý multiplikátor = okamžitý efekt / (1 - koeficient u zpožděné y)
ok$d_unemp_1 <- zpozdi(ok$d_unemp, 1)
adl <- lm(d_unemp ~ growth + d_unemp_1, data = ok[vzorek, ])
coef(adl)
coef(adl)[["growth"]] / (1 - coef(adl)[["d_unemp_1"]])

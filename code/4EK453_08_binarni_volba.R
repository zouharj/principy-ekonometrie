# =============================================================================
# 4EK453 Principy ekonometrie – Přednáška 8: Binární volba a numerika
# Jan Zouhar, Katedra ekonometrie FIS VŠE v Praze
#
# Skript reprodukuje všechny výsledky z prezentace. Spouštějte ho z kořene
# složky kurzu (data se čtou z podsložky data/).
# Řádky „# ---- název ----“ označují bloky, které se zobrazují v prezentaci.
#
# Potřebné balíčky: sandwich, lmtest, marginaleffects, brglm2
#   install.packages(c("sandwich", "lmtest", "marginaleffects", "brglm2"))
# =============================================================================

# ---- nastaveni ----
options(scipen = 10)     # raději víc číslic než vědecký zápis
library(sandwich)        # robustní směrodatné chyby pro LPM: vcovHC()
library(lmtest)          # coeftest(), lrtest()
library(marginaleffects) # mezní efekty: avg_slopes(), slopes()
library(brglm2)          # Firthova korekce při perfektní predikci


# =============================================================================
# Data: Global Entrepreneurship Monitor (GEM), Česko 2006 a 2011
# =============================================================================

# ---- gem-data ----
gem <- read.csv("data/GEM_CZE.csv", stringsAsFactors = TRUE)
# podnikani: v původních datech typ podnikání; 1 = podniká (začíná, nově
# nebo dlouhodobě), 0 = nepodniká ("zadne")
gem$podnikani <- as.numeric(gem$podnikani != "zadne")
# referenční kategorie: nejčastější vzdělání, bez sebedůvěry
gem$vzdelani   <- relevel(gem$vzdelani, "SS s mat.")
gem$sebeduvera <- relevel(gem$sebeduvera, "Ne")
gem$rok <- factor(gem$rok)
# pro odhady používáme jen řádky bez chybějících hodnot
f <- podnikani ~ sebeduvera + pohlavi + vzdelani + vek + I(vek^2) + rok
gem <- gem[complete.cases(gem[, all.vars(f)]), ]
c(n = nrow(gem), podil_podnikajicich = mean(gem$podnikani))


# =============================================================================
# Lineární pravděpodobnostní model (LPM)
# =============================================================================

# ---- lpm ----
lpm <- lm(f, data = gem)
# heteroskedasticita je v LPM zaručená: var(y | x) = p(x)[1 - p(x)]
coeftest(lpm, vcov = vcovHC(lpm, type = "HC1"))
# predikované pravděpodobnosti mimo interval [0, 1]
range(fitted(lpm))
sum(fitted(lpm) < 0)


# =============================================================================
# Metoda maximální věrohodnosti: Bernoulliho rozdělení bez vysvětlujících
# proměnných
# =============================================================================

# ---- ml-bernoulli ----
# L(p) = pravděpodobnost, že bychom znovu vylosovali výběr s přesně stejnými
# hodnotami y; maximum je v p = podíl jedniček
y <- gem$podnikani
logL_bern <- function(p) sum(y) * log(p) + sum(1 - y) * log(1 - p)
optimize(logL_bern, c(0.001, 0.999), maximum = TRUE)
mean(y)


# =============================================================================
# Logit a probit
# =============================================================================

# ---- logit ----
logit <- glm(f, data = gem,
             family = binomial(link = "logit"))
round(coef(summary(logit)), 4)

# ---- probit ----
probit <- glm(f, data = gem, family = binomial(link = "probit"))
# koeficienty logitu jsou zhruba 1.6 až 1.9krát větší než u probitu
round(cbind(logit = coef(logit), probit = coef(probit),
            podil = coef(logit) / coef(probit)), 3)

# ---- logL ----
# logaritmická věrohodnost: logit, probit a model jen s konstantou
logit0 <- glm(podnikani ~ 1, data = gem, family = binomial)
c(logit = logLik(logit), probit = logLik(probit), konstanta = logLik(logit0))
# McFaddenovo pseudo-R2
1 - logLik(logit) / logLik(logit0)


# =============================================================================
# Interpretace: poměr šancí a mezní efekty
# =============================================================================

# ---- pomer-sanci ----
# exp(beta) = poměr šancí
round(exp(coef(logit)), 3)

# ---- me-jedinec ----
# efekt sebedůvěry pro dva lidi: muž 45 let s VŠ, žena 60 let se ZŠ
osoby <- gem[c(1, 1), all.vars(f)]  # dva řádky dat, hodnoty přepíšeme
osoby$pohlavi[]  <- c("Muz", "Zena")
osoby$vzdelani[] <- c("VS", "pouze ZS")
osoby$vek <- c(45, 60)
osoby$rok[] <- "2011"
comparisons(logit, variables = "sebeduvera", newdata = osoby)

# ---- ame ----
# průměrné mezní efekty (AME), funkce z balíčku marginaleffects
avg_slopes(logit, variables = c("sebeduvera", "pohlavi"))

# ---- ame-probit ----
avg_slopes(probit, variables = c("sebeduvera", "pohlavi"))

# ---- mem ----
# mezní efekt pro "průměrné" pozorování (MEM): jiný pojem než AME
slopes(logit, variables = c("sebeduvera", "pohlavi"), newdata = "mean")


# =============================================================================
# Indukce: test věrohodnostním poměrem (LR)
# =============================================================================

# ---- lr-test ----
# sdružená významnost vzdělání: LR test
logit_bez_vzd <- update(logit, . ~ . - vzdelani)
lrtest(logit_bez_vzd, logit)

# ---- sep-uvod ----
# =============================================================================
# Perfektní predikce (separace)
# =============================================================================

# ---- sep-tabulka ----
# věk v intervalech po 10 letech
hranice <- c(18, 30, 40, 50, 60, 70, 80, 90)
nazvy <- c("18-29", "30-39", "40-49", "50-59",
           "60-69", "70-79", "80+")
gem$vek10 <- cut(gem$vek, hranice,
                 right = FALSE, labels = nazvy)
table(gem$vek10, gem$podnikani)

# ---- sep-logit ----
f10 <- podnikani ~ sebeduvera + pohlavi + vzdelani + vek10 + rok
logit10 <- glm(f10, data = gem, family = binomial)
# koeficienty u věkových intervalů (referenční 18-29)
round(coef(summary(logit10))[grep("vek10", names(coef(logit10))), ], 3)

# ---- sep-presnost ----
# Přísnější kritérium zastavení algoritmu: odhad u 80+ se změní, ostatní ne
logit10_p <- update(logit10, control = list(epsilon = 1e-14, maxit = 100))
round(cbind(standardni = coef(logit10), prisnejsi = coef(logit10_p)), 3)

# ---- sep-iterace ----
# Odhad u 80+ a logaritmická věrohodnost po jednotlivých iteracích
po_k_iteracich <- function(k) suppressWarnings(
  update(logit10, control = list(epsilon = 1e-30, maxit = k)))
iter_b80  <- sapply(1:25, function(k) coef(po_k_iteracich(k))[["vek1080+"]])
iter_logL <- sapply(1:25, function(k) as.numeric(logLik(po_k_iteracich(k))))
round(cbind(iterace = 1:25, b80 = iter_b80, logL = iter_logL), 3)

# ---- sep-testy ----
# Waldův test u 80+ nic neříká; LR test porovná model s intervalem 80+
# a model, kde je 80+ sloučeno s referenční kategorií 18-29
gem$vek10_bez80 <- gem$vek10
levels(gem$vek10_bez80)[levels(gem$vek10_bez80) == "80+"] <- "18-29"
logit10_h0 <- glm(podnikani ~ sebeduvera + pohlavi + vzdelani + vek10_bez80 + rok,
                  data = gem, family = binomial)
lrtest(logit10_h0, logit10)

# ---- sep-reseni ----
# Tři řešení: (1) sloučit intervaly 70-79 a 80+, (2) vyřadit 80+ (to dělá
# např. Stata), (3) Firthova korekce (penalizovaná věrohodnost)
gem$vek10_70 <- gem$vek10
levels(gem$vek10_70)[levels(gem$vek10_70) %in% c("70-79", "80+")] <- "70+"
r_sloucit <- glm(podnikani ~ sebeduvera + pohlavi + vzdelani + vek10_70 + rok,
                 data = gem, family = binomial)
r_vyradit <- glm(f10, data = droplevels(subset(gem, vek10 != "80+")),
                 family = binomial)
r_firth   <- glm(f10, data = gem, family = binomial, method = "brglmFit")
koef <- c("sebeduveraAno", "pohlaviZena")
round(rbind(sloucit = coef(r_sloucit)[koef], vyradit = coef(r_vyradit)[koef],
            firth = coef(r_firth)[koef], puvodni = coef(logit10)[koef]), 3)
round(coef(r_firth)["vek1080+"], 3)


# =============================================================================
# Numerika: jak software hledá maximum
# =============================================================================

# ---- irls ----
# glm(): Newtonova metoda (IRLS)
# trace vypíše deviaci = -2 log L
logit_t <- update(logit,
                  control = list(trace = TRUE))

# ---- optim-data ----
# Logaritmickou věrohodnost logitu napíšeme sami a maximalizujeme
# obecným optimalizátorem
X <- model.matrix(logit)
yv <- logit$y
minus_logL <- function(b, X) {
  eta <- X %*% b
  # plogis(eta, log.p = TRUE) = log(p) spočtený přesně i pro velká |eta|
  -sum(yv * plogis(eta, log.p = TRUE) + (1 - yv) * plogis(-eta, log.p = TRUE))
}

# ---- optim-skalovani ----
start <- rep(0, ncol(X))
o_puvodni <- optim(start, minus_logL, X = X, method = "BFGS")
# stejný model, věk v desítkách let (vek/10, vek^2/100)
X10 <- X
X10[, "vek"] <- X[, "vek"] / 10
X10[, "I(vek^2)"] <- X[, "I(vek^2)"] / 100
o_skal <- optim(start, minus_logL, X = X10, method = "BFGS")
c(konvergence_puvodni = o_puvodni$convergence, logL_puvodni = -o_puvodni$value,
  konvergence_skal = o_skal$convergence, logL_skal = -o_skal$value,
  logL_glm = as.numeric(logLik(logit)))
# odhad u sebedůvěry: optim bez škálování, se škálováním, glm()
c(o_puvodni$par[2], o_skal$par[2], coef(logit)[["sebeduveraAno"]])
# číslo podmíněnosti matice X (přednáška 6) před a po škálování
c(puvodni = kappa(X, exact = TRUE), skal = kappa(X10, exact = TRUE))


# =============================================================================
# Počítače nepočítají přesně
# =============================================================================

# ---- presnost ----
0.1 + 0.2 == 0.3
# p = 1 - 4e-18 se uloží jako přesně 1
plogis(40) == 1
log(1 - plogis(40))
# totéž přesně: log(1 - p) = -40
plogis(-40, log.p = TRUE)
# věrohodnost GEM: součin pravděpodobností
p_i <- ifelse(yv == 1, fitted(logit),
              1 - fitted(logit))
prod(p_i)
# dvojnásobek dat: podtečení na 0
prod(c(p_i, p_i))
# logaritmická věrohodnost: součet logaritmů
sum(log(p_i))

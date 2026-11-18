#' ---
#' title: "Vyhlazování pomocí regrese: polynomy, spliny, GAM a LOESS"
#' subtitle: "Doplňkový materiál k přednášce 7 (Výběr modelu a specifikace)"
#' author: "Jan Zouhar | Katedra ekonometrie FIS VŠE v Praze"
#' output:
#'   html_document:
#'     number_sections: yes
#'     toc: yes
#'     toc_float: yes
#' ---

#' Materiál ukazuje, jak lze pomocí obyčejné lineární regrese (OLS) odhadnout
#' nelineární vztah mezi dvěma proměnnými: polynomy, spliny, GAM a pro srovnání
#' LOESS. Vychází z přednášky ke kurzu 4EK621. Skript spouštějte z kořenové složky
#' kurzu (data se čtou z podsložky `data/`).

#+ include = FALSE
knitr::opts_chunk$set(fig.width = 9, fig.height = 5.5, dpi = 110, message = FALSE, warning = FALSE)

#' # Příprava

#' Balíčky:
library(tidyverse) # práce s daty a grafy (ggplot2)
library(splines)   # B-spliny bs() a přirozené spliny ns()
library(mgcv)      # zobecněné aditivní modely (GAM)

#' Nastavení: čísla bez vědeckého zápisu, barvy FIS VŠE a jednotný vzhled grafů.
options(scipen = 10)
fis_zelena <- "#00957d"
fis_seda   <- "#54585a"
# barvy pro čtyři úseky mezi uzly
barvy_useku <- c("#00957d", "#0072b0", "#54585a", "#7f4f9f")
theme_set(theme_minimal(base_size = 14))

#' # Data

#' ## Simulovaná data s nelineárním vztahem

#' Vygenerujeme 500 bodů: *x* rovnoměrně mezi 0 a 10 (seřazené), *y* jako součet
#' sinusovky, náhodného šumu, pomalu se měnící náhodné složky a lineárního trendu.
#' Skutečný tvar vztahu je tedy složitý a žádný jednoduchý funkční tvar ho přesně
#' nepopíše.
set.seed(1)
x <- runif(500, 0, 10) |> sort()
y <- 10 + sin(x) + rnorm(500, 0, 0.4) + sort(cumsum(rnorm(500, 0, 0.24))) - (1:500)/200
df <- data.frame(x, y)

ggplot(df, aes(x, y)) +
    geom_point(colour = "grey55", size = 1) +
    ggtitle("Simulovaná data")

#' # Polynomy odhadnuté pomocí OLS

#' ## Pomocná funkce pro graf

#' Funkce `graf_vzorce()` vykreslí data a regresní funkci odhadnutou pomocí `lm()`
#' se zadaným vzorcem `f` (včetně 95% pásu spolehlivosti pro střední hodnotu).
graf_vzorce <- function(f, nazev) {
    ggplot(df, aes(x, y)) +
        geom_point(colour = "grey55", size = 1) +
        geom_smooth(method = "lm", formula = f, colour = fis_zelena, linewidth = 1.5) +
        ggtitle(paste0(nazev, ": ", format(f)))
}

#' ## Polynomy různých stupňů

#' Každý vzorec dostane název, který se použije v nadpisu grafu. Funkce `poly()`
#' vytvoří ortogonální polynomy: vyrovnané hodnoty jsou stejné jako s `x + I(x^2)`,
#' jen jsou sloupce matice **X** navzájem nekorelované (numericky stabilnější,
#' viz téma *Počítače nepočítají přesně*).
polynomy <- list(
    Lineární          = y ~ x,
    Kvadratický       = y ~ x + I(x^2),
    Kvadratický_poly  = y ~ poly(x, 2),  # stejné vyrovnané hodnoty
    Kubický           = y ~ poly(x, 3),
    Polynom_4_stupně  = y ~ poly(x, 4),
    Polynom_5_stupně  = y ~ poly(x, 5)
)

#' Funkce `imap()` z balíčku purrr projde seznam a kromě vzorce předá i jeho název.
polynomy |> imap(graf_vzorce)

#' ## Co se odehrálo uvnitř `lm()`

#' Odhady samotné nás tu příliš nezajímají; koeficienty polynomů nemají názornou
#' interpretaci, důležitá je výsledná křivka.
polynomy |> map(\(f) lm(f, data = df) |> summary())

#' # Spliny odhadnuté pomocí OLS

#' **Spline** je funkce složená z polynomů, které na sebe navazují v bodech zvaných
#' **uzly** (*knots*). Zvolíme tři uzly, které rozdělí obor *x* na čtyři úseky.
uzly <- c(2.5, 5, 7.5)

#' ## Spline jako polynomy po částech s omezeními v uzlech

#' Nejdřív odhadneme v každém úseku zvlášť kubický polynom bez jakýchkoli omezení
#' a pak postupně požadujeme, aby funkce byla v uzlech spojitá, aby měla spojitou
#' první derivaci (bez zlomů) a nakonec spojitou druhou derivaci (zakřivení).
#' Poslední varianta je **kubický spline**.
#'
#' Technicky to jde pomocí tzv. useknutých mocnin $(x - k)_+^p$, které jsou nulové
#' vlevo od uzlu $k$ a rovny $(x - k)^p$ vpravo od něj (pro $p = 0$ jde o indikátor
#' $x > k$). K polynomu $x, x^2, x^3$ přidáme u každého uzlu useknuté mocniny
#' $p = 0, 1, 2, 3$, takže v každém úseku dostaneme vlastní polynom. Každé
#' vypuštění jedné mocniny u všech uzlů přidá jedno omezení: bez $p = 0$ nejsou
#' v uzlech skoky, bez $p = 1$ zlomy, bez $p = 2$ skoky ve druhé derivaci.
#' Každé omezení tak ubere jeden parametr na uzel.
usek_mocnina <- function(x, k, p) if (p == 0) as.numeric(x > k) else pmax(x - k, 0)^p

baze_po_castech <- function(x, vypustit) {
    B <- sapply(1:3, function(p) x^p)
    for (k in uzly) for (p in setdiff(0:3, vypustit)) B <- cbind(B, usek_mocnina(x, k, p))
    B
}

#' Pomocná funkce: odhadne OLS s danou maticí regresorů a vrátí vyrovnané hodnoty
#' na jemné mřížce bodů *x* (kvůli hladkému grafu), včetně čísla úseku mezi uzly.
mriz <- seq(0, 10, length.out = 801)
krivka <- function(baze, nazev) {
    B <- baze(x)
    m <- lm(y ~ B)
    tibble(x = mriz,
           yhat = drop(cbind(1, baze(mriz)) %*% coef(m)),
           usek = factor(findInterval(mriz, uzly) + 1),
           model = paste0(nazev, " (", length(coef(m)), " parametrů)"))
}

omezeni <- bind_rows(
    krivka(\(z) baze_po_castech(z, integer(0)), "1. bez omezení"),
    krivka(\(z) baze_po_castech(z, 0),          "2. spojitost"),
    krivka(\(z) baze_po_castech(z, 0:1),        "3. + spojitá 1. derivace"),
    krivka(\(z) baze_po_castech(z, 0:2),        "4. + spojitá 2. derivace = kubický spline")
)

graf_useku <- function(d) {
    ggplot(d, aes(x, yhat)) +
        geom_point(data = df, aes(x, y), colour = "grey75", size = 0.6) +
        geom_vline(xintercept = uzly, linetype = "dashed", colour = "grey50") +
        geom_line(aes(colour = usek, group = usek), linewidth = 1.3) +
        scale_colour_manual(values = barvy_useku, guide = "none") +
        facet_wrap(~ model, ncol = 2) +
        labs(y = "y")
}

#+ fig.height = 7
graf_useku(omezeni)

#' Poslední panel dává přesně totéž jako funkce `bs()` s kubickou bází, jen jsou
#' sloupce matice **X** jinak zvolené (B-spliny jsou numericky stabilnější než
#' useknuté mocniny):
m_ruc <- lm(y ~ baze_po_castech(x, 0:2))
m_bs  <- lm(y ~ bs(x, knots = uzly, degree = 3))
all.equal(fitted(m_ruc), fitted(m_bs))

#' ## Spliny různých stupňů

#' Stejný princip funguje pro libovolný stupeň: **lineární spline** je spojitá
#' lomená čára (jako indikátorové proměnné pro intervaly, jen se sklonem),
#' **kvadratický spline** má spojitou první derivaci. **Přirozený kubický spline**
#' (`ns()`) je navíc za krajními uzly lineární, takže se na okrajích dat nechová
#' tak divoce jako polynomy.
stupne <- bind_rows(
    krivka(\(z) bs(z, knots = uzly, degree = 1, Boundary.knots = c(0, 10)), "1. lineární spline"),
    krivka(\(z) bs(z, knots = uzly, degree = 2, Boundary.knots = c(0, 10)), "2. kvadratický spline"),
    krivka(\(z) bs(z, knots = uzly, degree = 3, Boundary.knots = c(0, 10)), "3. kubický spline"),
    krivka(\(z) ns(z, knots = uzly, Boundary.knots = range(x)),            "4. přirozený kubický spline")
)

#+ fig.height = 7
graf_useku(stupne)

#' Pro srovnání tytéž spliny i s pásy spolehlivosti přes `geom_smooth()`:
spliny <- list(
    BS_stupen_1 = y ~ bs(x, knots = uzly, degree = 1),
    BS_stupen_2 = y ~ bs(x, knots = uzly, degree = 2),
    BS_stupen_3 = y ~ bs(x, knots = uzly, degree = 3),
    NS_stupen_3 = y ~ ns(x, knots = uzly)
)
spliny |> imap(graf_vzorce)

#' ## Co se odehrálo uvnitř `lm()`

#' Koeficienty patří bázovým funkcím, ne původní proměnné *x*; jednotlivě je
#' neinterpretujeme.
spliny |> map(\(f) lm(f, data = df) |> summary())

#' ## Jak vypadají bázové funkce

#' Funkce `bs()` a `ns()` vytvoří z jedné proměnné *x* několik sloupců matice **X**
#' (bázových funkcí). Regresní funkce je jejich lineární kombinací s odhadnutými
#' koeficienty, takže jde stále o lineární regresní model. Vytáhneme sloupce
#' z matice modelu (bez interceptu) a vykreslíme je proti *x*.
baze_z_modelu <- function(m) {
    model.matrix(m) |>
        as.data.frame() |>
        select(-matches("Intercept")) |>
        set_names(\(n) paste0("báze ", seq_along(n))) |>
        mutate(x = df$x) |>
        pivot_longer(-x, values_to = "hodnota", names_to = "bazova_funkce")
}

graf_baze <- function(d, nazev) {
    ggplot(d, aes(x, hodnota, colour = bazova_funkce)) +
        geom_line(linewidth = 1) +
        geom_vline(xintercept = uzly, linetype = "dashed", colour = "grey50") +
        labs(colour = NULL) +
        ggtitle(nazev)
}

spliny |>
    map(\(f) lm(f, data = df) |> baze_z_modelu()) |>
    imap(graf_baze)

#' ## GAM: uzly a hladkost volí software

#' Zobecněný aditivní model (GAM) z balíčku mgcv nepotřebuje zadat uzly ani stupeň.
#' Použije bázi s mnoha uzly a „vlnitost“ funkce tlumí penalizací (podobně jako
#' ridge regrese v přednášce); sílu penalizace zvolí podle dat. Zde kubické
#' regresní spliny (`bs = "cr"`).
vzorec_gam <- y ~ s(x, bs = "cr")

ggplot(df, aes(x, y)) +
    geom_point(colour = "grey55", size = 1) +
    geom_smooth(method = "gam", formula = vzorec_gam, colour = fis_zelena, linewidth = 1.5) +
    ggtitle("Spline odhadnutý jako GAM (mgcv)")

#' Jak vypadá odhadnutý model? `ggplot()` volal funkci `gam()` z balíčku mgcv.
#' Hodnota `edf` (*effective degrees of freedom*) říká, kolika parametrům odpovídá
#' vyhlazená funkce po penalizaci.
gam_odhad <- gam(vzorec_gam, data = df)
summary(gam_odhad)  # kde jsou koeficienty u bázových funkcí?
coef(gam_odhad)     # tady

#' # LOESS

#' LOESS není regrese s pevným funkčním tvarem: v okolí každého bodu odhadne
#' váženou regresi (bližší body mají větší váhu). Parametr `span` určuje, jak
#' velká část dat tvoří okolí: velký span = hladká křivka (velké vychýlení, malý
#' rozptyl), malý span = křivka sleduje šum. Opět kompromis mezi vychýlením
#' a rozptylem.
graf_loess <- function(span) {
    ggplot(df, aes(x, y)) +
        geom_point(colour = "grey55", size = 1) +
        geom_smooth(method = "loess", formula = y ~ x, span = span,
                    colour = fis_zelena, linewidth = 1.5) +
        ggtitle(paste0("LOESS: span = ", span))
}

c(0.9, 0.7, 0.5, 0.3, 0.1) |> map(graf_loess)

#' # Binární a podílová vysvětlovaná proměnná

#' Data o obcích (Dvořák, Zouhar a Novák, 2017, *Political Research Quarterly*):
#' účast ve sněmovních volbách a místní referendum. Použijeme rok 2010.
#' Účast je v datech v procentech, převedeme ji na podíl.
obce <- read.csv("data/obce_v_CR.csv", fileEncoding = "CP1250") |>
    filter(rok == 2010) |>
    mutate(ucast = ucast / 100)

summary(obce)

#' U binární vysvětlované proměnné (konalo se referendum?) by přímka ani polynom
#' nemusely zůstat mezi 0 a 1. Proto vyhladíme pomocí logitu (podrobněji později
#' v kurzu). Body jsou svisle trochu rozházené (`geom_jitter()`), aby se nepřekrývaly.
binarni_vyhlazeni <- function(...) {
    geom_smooth(method = "glm", formula = y ~ x,
                method.args = list(family = "binomial"), ...)
}

ggplot(obce, aes(log_velikost, referendum)) +
    geom_jitter(height = 0.01, alpha = 0.25, colour = fis_seda) +
    binarni_vyhlazeni(colour = fis_zelena, linewidth = 1.5) +
    labs(x = "log(počet obyvatel)", y = "referendum") +
    ggtitle("Vyhlazení binární proměnné (logit)")

#' Podíl (účast) leží také mezi 0 a 1. Použijeme logit s kvazi-věrohodností
#' (*fractional logit*, Papke a Wooldridge, 1996), který podíl vysvětluje stejnou
#' funkcí jako pravděpodobnost u binární proměnné.
podilove_vyhlazeni <- function(...) {
    geom_smooth(method = "glm", formula = y ~ x,
                method.args = list(family = quasibinomial(link = "logit")), ...)
}

ggplot(obce, aes(log_velikost, ucast)) +
    geom_point(alpha = 0.25, colour = fis_seda) +
    podilove_vyhlazeni(colour = fis_zelena, linewidth = 1.5) +
    labs(x = "log(počet obyvatel)", y = "účast") +
    ggtitle("Vyhlazení podílu (frakční logit, kvazi-ML)")

#' Pro srovnání: přirozený spline odhadnutý pomocí OLS. Pro účast mezi zhruba 40
#' a 90 % zůstává v rozumných mezích a vztah zachytí pružněji než logit s lineárním
#' členem.
ggplot(obce, aes(log_velikost, ucast)) +
    geom_point(alpha = 0.25, colour = fis_seda) +
    geom_smooth(method = "lm", formula = y ~ ns(x, 5), colour = fis_zelena, linewidth = 1.5) +
    labs(x = "log(počet obyvatel)", y = "účast") +
    ggtitle("Přirozený spline (5 stupňů volnosti)")

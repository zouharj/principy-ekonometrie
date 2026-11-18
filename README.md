# 4EK453 Principy ekonometrie

Materiály k předmětu **4EK453 Principy ekonometrie** (Katedra ekonometrie FIS VŠE v Praze, ZS 2026/27, Jan Zouhar).

Tady najdete vždy **aktuální verze** prezentací, R skriptů, dat a dalších dokumentů. Každé PDF má na titulní straně datum verze. Oznámení, odevzdávání a hodnocení zůstávají v Moodlu a InSISu.

## Co tu je

| Složka | Obsah |
|---|---|
| [`prezentace/`](prezentace) | prezentace k přednáškám (PDF) |
| [`code/`](code) | R skripty k přednáškám: reprodukují všechna čísla, tabulky a grafy z prezentací |
| [`data/`](data) | data k přednáškám a cvičením (CSV) |
| [`cviceni/`](cviceni) | sbírka příkladů ke cvičení |
| [`zkouska/`](zkouska) | otázky a příklady k ústní zkoušce |
| [`prubezny_test/`](prubezny_test) | informace k průběžnému testu a ukázkový test s řešením |
| [`semestralni_prace/`](semestralni_prace) | pokyny k semestrální práci (PDF i zdrojový kód v LaTeXu) |
| [`matematicke_pozadi/`](matematicke_pozadi) | *Matematické pozadí*: náměty k dobrovolným konzultacím |
| [`vizualizace/`](vizualizace), [`doplnky/`](doplnky) | interaktivní materiály a doplňky (HTML, otevírejte přes odkazy níže) |

## Interaktivní materiály

Otevřou se v prohlížeči (fungují i na mobilu):

- [OLS jako ortogonální projekce](https://zouharj.github.io/principy-ekonometrie/vizualizace/viz_ols_projection.html)
- [Frischova-Waughova-Lovellova věta](https://zouharj.github.io/principy-ekonometrie/vizualizace/viz_fwl_theorem.html)
- [Výběrová rozdělení](https://zouharj.github.io/principy-ekonometrie/vizualizace/viz_sampling_distributions.html)
- [Intervaly spolehlivosti v opakovaných výběrech](https://zouharj.github.io/principy-ekonometrie/vizualizace/viz_ci_coverage.html)
- [Spliny a GAM (doplněk k přednášce 7)](https://zouharj.github.io/principy-ekonometrie/doplnky/4EK453_07_doplnek_spliny.html)

## Jak materiály stáhnout

- **Jednotlivý soubor:** klikněte na něj a pak na tlačítko *Download raw file* (ikona šipky vpravo nahoře).
- **Všechno najednou:** zelené tlačítko *Code* → *Download ZIP*.
- **S gitem (doporučeno):** `git clone https://github.com/zouharj/principy-ekonometrie`, aktualizace pak jen `git pull`. Git umí i RStudio a Positron (nový projekt ze správy verzí).
- **Upozornění na novinky:** tlačítko *Watch* nahoře (vyžaduje účet na GitHubu).

## Jak spouštět R skripty

Skripty čtou data relativní cestou, např. `read.csv("data/wage2.csv")`. Otevřete proto celou složku jako projekt: v RStudiu soubor `principy-ekonometrie.Rproj`, v Positronu *Open Folder*. Pak skripty fungují bez úprav. Některé soubory (např. `obce_v_CR.csv`, `skoda.csv`) jsou v kódování Windows-1250, skripty s tím počítají (`fileEncoding = "CP1250"`).

Na začátek vlastních skriptů doporučuji `options(scipen = 10)`, aby R nevypisovalo čísla ve vědeckém zápisu.

Data z Českého panelového šetření domácností o populismu (`populismus_cisty.csv`) tu nejsou: podmínky Českého sociálněvědního datového archivu nedovolují jejich další šíření. Najdete je v Moodlu (jen pro účastníky kurzu).

## Chyby a dotazy

Najdete-li v materiálech chybu, napište mi e-mail nebo založte *Issue* zde na GitHubu.

Materiály jsou určené studentům předmětu; jejich další šíření konzultujte s autorem.

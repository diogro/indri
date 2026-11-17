# LG/SM F3 intercross: trait classes

The columns of
[lgsm_f3](https://diogro.github.io/indri/reference/lgsm_f3.md) in each
trait class.

## Usage

``` r
lgsm_traits
```

## Format

A named list of character vectors: `growth` (`growth12` to `growth910`,
the nine weekly gains; the original QTL analyses used the first seven),
`weight` (`WEEK1` to `WEEK10`, raw weekly weights), `area` (`area1` to
`area7`) and `necropsy` (`FATPAD`, `HEART`, `KIDNEY`, `SPLEEN`,
`LIVER`).

## Source

Built by `data-raw/lgsm.R`.

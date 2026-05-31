# Read genotype dosage matrix

Read a genotype matrix in dosage format (0/1/2, NA for missing). Rows
are individuals, columns are SNP markers.

## Usage

``` r
read_genotype(file, sep = ",", header = TRUE)
```

## Arguments

- file:

  Path to a CSV/TSV file. First column is individual ID, remaining
  columns are SNP dosage values (0, 1, 2, or NA).

- sep:

  Field separator (default: comma).

- header:

  Logical; does the file contain a header row? Default TRUE.

## Value

A numeric matrix with rownames = individual IDs, colnames = marker
names.

## See also

Other data-io:
[`read_anchors()`](https://luansheng.github.io/AAPA/reference/read_anchors.md),
[`read_parents()`](https://luansheng.github.io/AAPA/reference/read_parents.md),
[`simulate_aapa_data()`](https://luansheng.github.io/AAPA/reference/simulate_aapa_data.md)

## Examples

``` r
# Create a temporary genotype file
tmp <- tempfile(fileext = ".csv")
geno_data <- data.frame(
  id = c("IND1", "IND2"),
  SNP1 = c(0, 1),
  SNP2 = c(2, 1)
)
write.csv(geno_data, tmp, row.names = FALSE)
geno <- read_genotype(tmp)
#> ℹ Reading genotype file: /tmp/Rtmp5VLxv7/file196d7be4443.csv
#> ✔ Read 2 individuals x 2 markers
unlink(tmp)
```

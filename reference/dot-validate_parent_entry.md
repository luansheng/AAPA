# Validate a single parent entry

Check that one family entry in an `aapa_parents` object satisfies the
required field, naming, and marker constraints.

## Usage

``` r
.validate_parent_entry(family, family_name)
```

## Arguments

- family:

  A single parent entry.

- family_name:

  Expected family name from `names(parents)`.

## Value

The input family entry, invisibly.

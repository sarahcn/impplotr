---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# impplotr

The goal of impplotr is to generate plots and summaries from imputation metrics provided by the University of Michigan Imputation Server. It is currently optimized to work with the Haplotype Reference Consortium (HRC r1.1) reference panel. 

## Example

While you can run each plotting and summarizing function separately, the convenience function `makeAllPlots` will run them all in one step, e.g.


```r

# imp_dir should contain metrics files ("chr#.info.gz")
# out_dir is desired location of plots and summaries files
# start_chr is desired start chromosome (i.e., first metrics file to read)
# end_chr is desired end chromosome (i.e., last metrics file to read)

makeAllPlots(imp_dir="metrics_directory_name",
			out_dir="output_directory_name",
			start_chr=1, end_chr=22)
```

# Planned updates
I created this package from a set of functions I've used internally for many years to process and report imputation metrics, initially from local imputation with IMPUTE2 and 1000 Genomes Phase 3, and more recently with the UM server. My goal is to help others who want to parse the large amounts of metrics data generated by contemporary imputation. However, I recognize several things could be improved with the current release. My to do list includes:

1. remove `start_chr` and `end_chr` arguments from the individual plotting and summarizing functions, as they currently only go in the output file name(s)
2. give the user more granular control over output file names (i.e., rather than specifying output directory with `out_dir`, provide the base file path AND name

I welcome other suggestions via GitHub issues tracking!
~Sarah

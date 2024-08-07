# biome_gs

Code repository to run analyses for Nishiguchi et al. "Ecological features affecting patterns of speciation among terrestrial mammals."

The code is in [R](https://cran.r-project.org/). The [targets package](https://wlandau.github.io/targets/index.html) is used to manage the workflow.

To run all analyses:

1. [Clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
2. [Prepare data](#Data)
3. [Run the analysis pipeline](#Analysis-pipeline)

## Data

Before running the analysis pipeline, you need to prepare some data and code in the following procedures:

- [Download spatial data of mammals](https://www.iucnredlist.org/resources/spatial-data-download). Then, place the files in the `data/MAMMALS` folder.
- [Download Köppen-Geiger climate classification maps](https://doi.org/10.6084/m9.figshare.6396959). Then, place the files in the `data/Beck_KG_V1` folder.
- [Download DR statistics calculated from 10,000 node-dated trees](https://doi.org/10.5061/dryad.tb03d03). Then, place the file in the `data/tree` folder.
- [Download the code for FiSSE](https://github.com/macroevolution/fisse/blob/master/run_fisse/traitDependent_functions.R). Then, place the file in the `R` folder. Please note that the code for calculating p-values was modified from `sum(delta_true > delta) / (reps + 1)` to `min(sum(delta_true > delta) / reps, 1 - sum(delta_true > delta) / reps) * 2` to derive p-values presented in our manuscript.

[The mammalian tree](https://doi.org/10.1371/journal.pbio.3000494) will be automatically downloaded through [the analysis pipeline](#Analysis-pipeline). For more information, please see our manuscript.

## Analysis pipeline

[A docker image is provided](https://hub.docker.com/r/aurelia01/biome_gs) to run the pipeline. You can [install docker from here](https://docs.docker.com/install/).

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine):

```
cd /path/to/repo
```

Run the pipeline:
```
docker run --rm -v ${PWD}:/tmpdir -w /tmpdir aurelia01/biome_gs:1.3 Rscript -e 'targets::tar_make()'
```

Or run the pipeline in parallel (e.g. with 8 workers):
```
docker run --rm -v ${PWD}:/tmpdir -w /tmpdir aurelia01/biome_gs:1.3 Rscript -e 'targets::tar_make_clustermq(workers=8)'
```

The output files will be written to the `data/output` folder.

## License

The code in this repository is licensed under the [MIT license](LICENSE)
# biome_gs

Code repository to run analyses for Nishiguchi et al. "Two contrasting evolutionary drivers shape biome niche breadths of terrestrial mammals."

The code is in [R](https://cran.r-project.org/). The [targets package](https://wlandau.github.io/targets/index.html) is used to manage the workflow.

To run all analyses:

1. [Clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
2. [Prepare the input data](#Data)
3. [Run the analysis pipeline](#Analysis-pipeline)

## Data

Before running the analysis pipeline, you need to prepare some data in the following procedures:

- [Download spatial data of mammals](https://www.iucnredlist.org/resources/spatial-data-download). Then, place the files in the `data/MAMMALS` folder.
- [Download spatial data of birds](http://datazone.birdlife.org/species/requestdis). Then, place the files in the `data/Aves` folder.
- [Download the Köppen-Geiger climate classification maps](https://doi.org/10.6084/m9.figshare.6396959). Then, place the files in the `data/Beck_KG_V1` folder.
- [Download 10000 avian trees of Hackett backbone at stage 1](https://data.vertlife.org) and concatenate these files into a single file. Then, convert the Newick format to the Nexus format. Using the tree file as an input, run [TreeAnnotator](https://beast.community/treeannotator) and generate the maximum clade credibility (MCC) tree of birds. Finally, place the MCC tree file in the `data/tree` folder.

Note that [the mammalian tree](https://doi.org/10.1371/journal.pbio.3000494) will be automatically downloaded through [the analysis pipeline](#Analysis pipeline). For more information, please see the manuscript.

## Analysis pipeline

[A docker image is provided](https://hub.docker.com/r/aurelia01/biome_gs) to run the pipeline. You can [install docker from here](https://docs.docker.com/install/).

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine):

```
cd /path/to/repo
```

Run the pipeline:
```
docker run --rm -v ${PWD}:/tmpdir -w /tmpdir aurelia01/biome_gs:1.2 Rscript -e 'targets::tar_make()'
```

Or run the pipeline in parallel (e.g. with 4 workers):
```
docker run --rm -v ${PWD}:/tmpdir -w /tmpdir aurelia01/biome_gs:1.2 Rscript -e 'targets::tar_make_clustermq(workers=4)'
```

The output files will be written to the `data/output` folder.

## License

The code in this repository is licensed under the [MIT license](LICENSE)
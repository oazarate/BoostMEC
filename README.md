# BoostMEC: Boosting and Markov for Efficient CRISPR
  
A pipeline for predicting CRISPR-Cas9 cleavage efficiency using feature engineering and [LightGBM](https://lightgbm.readthedocs.io/).

## Requirements
  
The software requirements for BoostMEC are listed below. To simplify the management of software packages, we also provide a [Docker image](#using-the-docker-image) with the required packages and code.

R 4.1.0  
Python 3.7.1  
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) 2.4.14

The Python and R package requirements are available [here](docker_files/requirements.txt) and [here](docker_files/requirements.r), respectively.

## Instructions for use

BoostMEC requires that the files under the [model](model) directory always reside in a shared location.

### Data format
To use BoostMEC, sgRNA target regions 30 nt in length are required, consisting of the 4 nt context in the 5' end + 20 nt sgRNA + 3 nt PAM + 3 nt context in the 3' end. Data should be formatted as a CSV file with 2 required columns: dataset and x30mer, which specify the sgRNA source dataset name and the target region 30-mer, respectively. We provide a small example dataset comprised of sequences from [Doench et al. 2014](https://www.nature.com/articles/nbt.3026) and [Moreno-Mateos et al. 2015](https://www.nature.com/articles/nmeth.3543) (data sourced from [this repo](https://github.com/maximilianh/crisporPaper)) in [data/example-test-sequences.csv](data/example-test-sequences.csv)

### Predicting efficiency
To obtain predictions with BoostMEC, run the Python script [model/boostmec_pipeline.py](model/boostmec-pipeline.py) with the required `--file` argument, providing the location and name of the CSV file with the sgRNA 30-mers of interest. These must be formatted as specified above. From within the BoostMEC repo, we can obtain predictions for the provided example dataset as below:

```console
python model/boostmec-pipeline.py --file data/example-test-sequences.csv
```
This command will produce three output files in the same location as the provided CSV file:
  
* __\[Data file name\]-with-features.csv__: This appends the original sgRNA target region CSV with all the extracted features. Note that the column names for position-specific k-mer features are numbered according to their absolute position on the 30-mer here, as opposed to the position relative to the sgRNA 20-mer target used in graphs.

* __\[Data file name\]-with-features-encoded.csv__: This contains the encoded version of the features used in BoostMEC, for use with LightGBM in R, and is useful in the creation of [interpretation plots](#interpretation-plots).

* __\[Data file name\]-boostmec-predictions.csv__: This contains the BoostMEC editing efficiency predictions for each sgRNA in the original source file.

### Interpretation plots
  
BoostMEC provides a wrapper for LightGBM's `lgb.interprete()` function that simplifies the creation of interpretation plots for sgRNA efficiency predictions. These interpretation plots provide the overall contribution of each feature (positive or negative) for a specific prediction. By default, only the top 10 features are shown (this can be modified in code).
  
To create an interpretation plot, the BoostMEC pipeline must have already run and produced the encoded feature CSV file. The required arguments are `--path` which provides the path for the BoostMEC model objects directory, `--dataset` which provides the location and filename for the encoded feature CSV file, and `--row` which specifies which row contains the prediction of interest. The command below will produce an interpretation plot for the first sgRNA in the example dataset under the [data](data) directory. Rows are 1-indexed (R convention).

```console
Rscript model/interpretation.R --path model --dataset data/example-test-sequences-with-features-encoded.csv --row 1
```

### Visualizing model trees
  
BoostMEC also provides a wrapper for LightGBM's `create_tree_digraph()` function that simplifies the creation of tree visualizations for each tree in the BoostMEC model. In the code below, we produce a visualization for the first tree in the model. Note that trees are 0-indexed (Python convention).

```console
python model/plot_trees.py --index 0
```

### Using the Docker image
  
The most convenient and reliable way to use BoostMEC is through our Docker image, which contains all the required software to run the pipeline. All the required files from this repo are also available in the BoostMEC folder in the image. To use use BoostMEC through Docker, first ensure that you have [Docker](https://www.docker.com/) installed. Then, pull the image and run it in interactive mode (creating your container at the same time):

```console
docker pull oazarate/boostmec:latest
docker run -ti oazarate/boostmec:latest
``` 

To exit the container, you can use the `exit` command. You can re-access the container using its randomly generated container name, which we will denote via [container name]. To find all your containers and their names, run:

```console
docker container ls -a
```
Re-start your container with:

```console
docker start -i [container name or id]
```

To move files to your container:

```console
docker cp your-file.csv [container_name]:/container-dir/
```

To move files from your container:

```console
docker cp [container name]:/container-dir/your-file your-host-dir
```

The dockerfile and package requirements are also available in this repository under [docker_files](docker_files).

## Model training data
  
BoostMEC was trained on high-throughput CRISPR efficiency data from from [Kim et al. 2019](https://www.science.org/doi/10.1126/sciadv.aax9249) and [Xiang et al. 2021](https://www.nature.com/articles/s41467-021-23576-0).

## Citations

If you use BoostMEC in your work, please cite our publication:

Zarate, O.A., Yang, Y., Wang, X. et al. BoostMEC: predicting CRISPR-Cas9 cleavage efficiency through boosting models. BMC Bioinformatics 23, 446 (2022). https://doi.org/10.1186/s12859-022-04998-z

In addition, please cite the following publications, which created the datasets we used for training:

Kim HK, Kim Y, Lee S, Min S, Bae JY, Choi JW, et al. SpCas9 activity prediction by DeepSpCas9, a deep learning-based model with high generalization performance. Sci Adv. 2019;5(11):9249.

Xiang X, Corsi GI, Anthon C, Qu K, Pan X, Liang X, et al. Enhancing CRISPR-Cas9 gRNA efficiency prediction by data integration and deep learning. Nat Commun. 2021;12(1):3238.
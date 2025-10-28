# Installation
Welcome to the installation guide! Follow the steps below to set up ZipStrain on your system

## Install using Conda
To install ZipStrain using Conda, follow these steps:

1- Make sure you have Conda installed on your system. If you don't have it yet, you can download and install Miniconda from [here](https://docs.conda.io/en/latest/miniconda.html).

2- Open your terminal and create a new Conda environment for ZipStrain:

```
conda create -n zipstrain_env python=3.13
```

3- Activate the newly created environment:

```
conda activate zipstrain_env
```

4- Install ZipStrain using Conda:

```
conda install -c bioconda zipstrain
```

This will install all the necessary dependencies along with ZipStrain.

5- Make sure the installation was successful by checking the version of ZipStrain:

```
zipstrain test
```

## Install using pip

To install ZipStrain using pip, follow these steps:

1- Make sure you have Python 3.13 or higher installed in your Python environment and activate it.

2- Open your terminal and run the following command to install ZipStrain:

```
pip install zipstrain
```

This will install all the necessary dependencies along with ZipStrain except samtools. Refer to samtools website for installation instructions: http://www.htslib.org/download/

3- Make sure the installation was successful by testing ZipStrain:

```
zipstrain test
```

## Use Docker

To use ZipStrain with Docker, follow these steps:

1- Make sure you have Docker installed on your system. 

2- Run ZipStrain using the Docker image:

```
docker run -it parsaghadermazi/zipstrain:latest zipstrain test
```

***Note***: Apptainer works in a very similar way to Docker:

```
apptainer run docker://parsaghadermazi/zipstrain:latest zipstrain test
```

## Setting up Nextflow (Optional)

You can use Conda to install Nextflow:

```
conda install -c bioconda nextflow
```

For more installation options, refer to the Nextflow installation guide: https://www.nextflow.io/docs/latest/getstarted.html#installation

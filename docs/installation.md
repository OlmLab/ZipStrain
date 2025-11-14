# Installation
Welcome to the installation guide! Follow the steps below to set up ZipStrain on your system


## Install using pip

To install ZipStrain using pip, follow these steps:

1- Make sure you have Python 3.12 or higher installed in your Python environment and activate it.

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


# snipar

*snipar* (single nucleotide imputation of parents) is a Python package for inferring identity-by-descent (IBD) segments shared between siblings, imputing missing parental genotypes, and for performing
family based genome-wide association and polygenic score analyses using observed and/or imputed parental genotypes. 

The imputation method and the family-based GWAS and polygenic score models are described in [Young et al. 2022](https://www.nature.com/articles/s41588-022-01085-0).

THe robust estimator, unified estimator and sib-difference estimator are described in [Guan et al. 2022](https://www.biorxiv.org/content/10.1101/2022.10.24.513611v1).

# Main features:

Infer identity-by-descent segments shared between siblings (ibd.py). 

Impute missing parental genotypes given the observed genotypes in a nuclear family (impute.py).

Perform family based GWAS using various estimators (gwas.py). 

Compute polygenic scores for probands, siblings, and parents from SNP weights using observed/imputed parental genotypes, and perform family
 based analysis of polygenic scores (pgs.py script). 
 
 Compute genome-wide correlations between different effects estimated by gwas.py (correlate.py). 

 Simulate genotype and phenotypes under different scenarios: direct and indirect genetic effects, vertical transmission, assortative mating (simulate.py). 

# Documentation

Documentation: https://snipar.rtfd.io/

It is recommended to read the guide: https://snipar.rtfd.io/en/latest/guide.html

And to work through the tutorial: https://snipar.rtfd.io/en/latest/tutorial.html

# Installing Using pip

*snipar* currently supports Python 3.7-3.9 on Linux, Windows, and Mac OSX. We recommend using a python distribution such as Anaconda 3 (https://store.continuum.io/cshop/anaconda/). 

# Virtual Environment

You may encounter problems with the installation due to Python version incompatability or package conflicts with your existing Python environment. To overcome this, you can try installing in a virtual environment. In a bash shell, this could be done by using the following commands in your directory of choice:
    
    python -m venv path-to-where-you-want-the-virtual-environment-to-be

You can activate and use the environment using

    source path-to-where-you-want-the-virtual-environment-to-be/bin/activate

# Installing From Source
To install from source, clone the git repository, and in the directory
containing the *snipar* source code, at the shell type:

    pip install .

# Python version incompatibility 

*snipar* does not currently support Python 3.10 or higher due to version incompatibilities of dependencies. 
To overcome this, create a Python3.9 environment using conda and install using pip in the conda environment:
	
    conda create -n myenv python=3.9
	conda activate myenv
	pip install snipar

# Apple ARM processor machines

There can be difficulties install *snipar* on Apple ARM processor machines due
to lack of available versions of scientific computing software made for these
processors' architectures. A workaround for this is to use *Snipar* in a docker
container.

The following steps will guide you on how to create a suitable container for
your M1/M2 MacBook and seamlessly integrate it into a VSCode environment. These
steps assume you have little knowledge about Docker, so we will start from
scratch.

## 1. Installing Docker Engine

Ensure that Docker Engine is installed on your machine. If it's already
installed, you can skip this step. Otherwise, you can install it by following
the instructions provided at the following link:

[Install Docker Engine on a MacBook](https://docs.docker.com/desktop/install/mac-install/)

## 2. Creating the Docker Container

To install the Snipar package, you need to create a Docker container that
emulates a virtual Ubuntu machine where Snipar can be installed.

To create the appropriate Docker container, follow these steps:

- Start the Docker engine. Just open the Docker application from your
  Applications folder.

- While the engine is running, open your terminal and execute the following
  command to create a container named "snipar_container" (you can choose a
  different name if you prefer):

  ```bash
  docker run --name snipar_container -it amd64/python:3.9.9-slim-buster /bin/bash
  ```

After running this command, you should see the "snipar_container" listed in the
Docker GUI under the "Containers" tab.

- You can close the terminal at this point.

## 3. Running the Container

To use the environment of the created container within VSCode or other IDEs,
ensure that the container is running. You can do this easily through the Docker
GUI:

- Open the Docker GUI.
- In the "Containers" section of the dashboard, locate "snipar_container" (or
  the name you chose for your container).
- Click the play button to start the container.

Once the container is running, you can access its terminal and files through the
Docker GUI. Keep in mind that the files within the container are isolated from
your macOS files. You won't be able to access them using Finder, but you can
manage them through the Docker GUI, including uploading and downloading.

## 4. Attaching the Container to VSCode

If you prefer a smoother development experience with VSCode, follow these steps
to attach the container to VSCode:

- Install the "Dev Containers" extension by visiting the following link in the
  VS Marketplace:

  [Dev Containers Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

- Once installed, a new icon will appear on the left sidebar of VSCode labeled "
  Remote Explorer."

- Click on "Remote Explorer" and, under the "Dev Containers" section, locate "
  snipar_container."

- Next to the container name, you'll find a button that says "Attach in a new
  window." Click on this button, and VSCode will open a new window attached to
  the container.

In the newly attached VSCode window, the terminal will be connected to the
container, similar to the "Exec" tab in the Docker GUI. You can also open
folders from the container environment (not your Mac itself) using the "Open
Folder" button in VSCode. This makes it more convenient to manage files while
writing code, and you can run your modules using the terminal directly within
VSCode.

## 5. Installing *snipar* Package in the Container

Up to this point, you've set up the environment for snipar but haven't installed
it in the container. To install snipar, follow these steps:

- Open the terminal in the VSCode window attached to "snipar_container" (or use
  the "Exec" tab in the Docker GUI).
- Run the following commands:
  ```bash
  pip install --upgrade pip
  pip install snipar
  ```
   
# Running tests
To check that the code is working properly and that the C modules have been compiled, you can run the tests using this command:

    python -m unittest snipar.tests
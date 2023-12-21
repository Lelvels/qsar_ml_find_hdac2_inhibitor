# Harmonizing QSAR machine learning-based models and docking approaches for the identification of novel HDAC2 inhibitors

by Dao Quang Tung, Do Thi Mai Dung, Nguyen Thanh Cong, Dao Ngoc Nam Hai, Daniel Baecker, Phan Thi Phuong Dung, Nguyen Hai Nam, Nguyen Ngoc An*

*Correspondence: [ngocan@vnu.edu.vn](mailto:ngocan@vnu.edu.vn) (N.N.A)

This repository is the prove of our work for our paper, which has been submitted for publication in Some Journal, doi link.

### Dependencies and implementation

You will need a working Python environment to run the code. The recommended way to set up your environment is through the [Anaconda Python distribution](https://www.anaconda.com/download/) which provides the `conda` package manager. Anaconda can be installed in your user directory and does not interfere with the system Python installation. The required dependencies are specified in the file `environment.yml`.

Run the following command in the repository folder (where `environment.yml` is located) to create a separate environment and install all required dependencies in it:

```
conda env create -f environment.yml
```

Then verify that the new environment was installed correctly:

```
conda env list
```

Our screening dataset was stored using PostgreSQL database, installation is availible in [PostgreSQL official website](https://www.postgresql.org/). Any request for the screening dataset, please contact with the correspondence [ngocan@vnu.edu.vn](mailto:ngocan@vnu.edu.vn).

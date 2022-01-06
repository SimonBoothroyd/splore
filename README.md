## Scroll through and exPLORE molecule sets

[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/SimonBoothroyd/splore.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/SimonBoothroyd/splore/context:python)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The `splore` framework aims to offer a simple graphical interface for scrolling through and exploring data sets of 
molecules.

![The GUI provided by `splore`](docs/gui.png)

### Getting Started

The GUI can easily be launched from the command line using the `splore` command:

```shell
# Load molecules from a local file
splore --file path-to-molecules.sdf

# Load molecules from a public QCArchive dataset
splore --qcf-dataset "OpenFF BCC Refit Study COH v1.0" --qcf-datatype basic
splore --qcf-dataset "OpenFF Rowley Biaryl v1.0" --qcf-datatype td
```

A full list of options can be printed using the `--help` flag:

```shell
splore --help                                                                   
Usage: splore [OPTIONS]

Options:
  --file FILE     The path to the file of molecules (.smi, .sdf, .sdf.gz) to
                  display.  [required]
  --port INTEGER  The port to run the GUI on.  [default: 8000; required]
  --help          Show this message and exit.

```

### Installation

The required dependencies for this framework can be installed using `conda`:

```shell
conda env create --name splore --file devtools/conda-envs/test-env.yaml
python setup.py develop
```

after which the GUI can be built by running:

```shell
cd frontend
npm install
npm run build -- -c production --output-path ../splore/_static --resources-output-path --deploy-url static/
cd ..
```

### License

The main package is release under the [MIT license](LICENSE). 

### Copyright

Copyright (c) 2021, Simon Boothroyd

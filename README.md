# ExtraCellTalk
`ExtraCellTalk` is a Python-based tool for the prediction of long-distance cell communication mediated by extracellular proteins.

![](https://github.com/xupeng-zang/ExtraCellTalk/blob/main/ExtraCellTalk.jpg)

## Contents

- [Overview](#Overview)
- [Repo Contents](#Repo-Contents)
- [System Requirements](#System-Requirements)
- [Installation Guide](#Installation-Guide)
- [Demo](#Demo)
- [License](#License)
- [Citation](#Citation)

## Overview

`ExtraCellTalk`, a tool for predicting long-distance cell communication mediated by extracellular proteins as intermediate bridges. The tool is based on the prediction of protein ligands, receptors and their interactions, similar to other previously released prediction tools for cell communication, but extracellular proteins are introduced as intermediate bridges for prediction. 

## Repo Contents

- [data](https://github.com/xupeng-zang/ExtraCellTalk/tree/main/data): `ExtraCellTalk` built-in database files and sample data.
- [res](https://github.com/xupeng-zang/ExtraCellTalk/tree/main/res): The results obtained after running `ExtraCellTalk` program on the sample data.
- [ExtraCellTalk.ipynb](https://github.com/xupeng-zang/ExtraCellTalk/blob/main/ExtraCellTalk.ipynb): Example of the `ExtraCellTalk` program in action.

## System Requirements

#### Hardware requirements

`ExtraCellTalk` requires only a standard computer with enough RAM to support the in-memory operations.

#### Software requirements

##### OS Requirements

This tool is supported for Windows and Linux. The tool has been tested on the following systems:

- Windows: 11
- Linux: 

##### Python Dependencies

`ExtraCellTalk` mainly depends on the Python scientific stack.

```python
numpy
pandas
pyarrow
matplotlib
statsmodels
rpy2
```

##### R Dependencies

The implementation of `SankeyPlot` function needs to call R in Python and depends on the following R packages.

```
ggplot2
ggalluvial
RColorBrewer
```

## Installation Guide

#### Install from Github

```python
git clone https://github.com/xupeng-zang/ExtraCellTalk.git
```

## Demo

For interactive demos of the functions, please check out the built-in example file [ExtraCellTalk.ipynb](https://github.com/xupeng-zang/ExtraCellTalk/blob/main/ExtraCellTalk.ipynb).

## License

This project is covered under the **MIT License**.

## Citation

For usage of the tool and associated manuscript, please await the publication of manuscript.

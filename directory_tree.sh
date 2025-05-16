#!bash.
├── Dockerfile
├── README.md
├── api
├── data
├── docs
│   ├── Makefile
│   ├── build
│   └── source
│       ├── README.md -> ../../README.md
│       ├── _static
│       ├── _templates
│       ├── conf.py
│       ├── index.rst
│       ├── modules.rst
│       ├── process.rst
│       ├── scsilhouette.cli.rst
│       ├── scsilhouette.compute.cli.rst
│       ├── scsilhouette.compute.rst
│       ├── scsilhouette.download.cli.rst
│       ├── scsilhouette.download.rst
│       ├── scsilhouette.nsforest.cli.rst
│       ├── scsilhouette.nsforest.rst
│       ├── scsilhouette.rst
│       ├── scsilhouette.viz.cli.rst
│       └── scsilhouette.viz.rst
├── environment.yml
├── main.nf
├── nextflow.config
├── pyproject.toml
├── results
├── setup.cfg
├── src
│   ├── scsilhouette
│   │   ├── __init__.py
│   │   ├── cli.py
│   │   ├── compute.py
│   │   ├── download.py
│   │   ├── nsforest.py
│   │   ├── utils.py
│   │   └── viz.py
│   └── scsilhouette.egg-info
│       ├── PKG-INFO
│       ├── SOURCES.txt
│       ├── dependency_links.txt
│       ├── entry_points.txt
│       └── top_level.txt
└── tests
    ├── test_cli.py
    ├── test_compute.py
    ├── test_download.py
    └── test_viz.py

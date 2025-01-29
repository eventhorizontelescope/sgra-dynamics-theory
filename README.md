# Sgr A* Dynamic Imaging Theory Tools

This repository collects software tools used for theoretical
interpretations in the **Sgr A\* Dynamic Imaging Project**.
Rather than a single cohesive `python` package, it serves as a
collection of independent scripts written in various programming
languages, leveraging different libraries and packages as needed.

## Usage

This repository may be cloned and run locally or in any unix-like
systems.
To ensure a standardized workflow:

+ **Simulation data and analytical models** should be placed or
  symlinked in the `models/` directory.

+ **Observational data and related inputs** should be placed or
  symlinked in the `data/` directory.

+ **Intermediate data products**, such as position angles measured
  from simulations or other processed outputs, may be placed in the
  `cache/` directory.

Each script should follow this directory structure for consistency and
ease of automation.

## Structure

Each **observable** analyzed in the project corresponds to a dedicated
directory.
Within each directory:

+ Include a `README.md` file detailing the purpose, dependencies, and
  usage of the scripts.

+ Provide a `run.sh` script that, when executed **without arguments**,
  produces the **default analysis** included in the **Sgr A\***
  papers.
  This ensures automation and reproducibility.

+ Place scripts and supporting files necessary for computation and
  analysis.

## Documentation

+ **Implementation details, equations, and algorithmic explanations**
  should be documented in the **"Notes for Theoretical Interpretations
  in the Sgr A\* Dynamic Imaging Project"** LaTeX document, which will
  become part of the dynamic imaging collaboration and/or individual
  official papers.

+ The repository's `README.md` files should focus on practical usage
  rather than detailed theoretical background.

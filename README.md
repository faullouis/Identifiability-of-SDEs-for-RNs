# Identifiability of SDEs for reaction networks

This repository contains the source code and a tutorial associated to our paper :

**"Identifiability of SDEs for reaction networks"**  
_Louis Faul, Linard Hoessly and Panqiu Xia (2025)_

##  Project Description

We implement methods to decide identifiability of ODEs and SDEs for reaction networks. The repository includes Julia scripts to:
- Verify rate-identifiability of a reaction network both w.r.t. its ODE and w.r.t. its SDE, and if so give the corresponding rate constants.
- Decide if two reaction networks are confoundable both w.r.t. their ODEs and w.r.t. their SDE, and if so give the corresponding rate constants.
- Decide if two reaction networks are linearly conjugated both w.r.t. their ODEs and w.r.t. their SDEs, and if so give the corresponding rate constants.

  See the file Tutorial.ipynb for some useful examples.
---

## Installation

Clone the repository:

```bash
git clone https://github.com/faullouis/Identifiability-of-SDEs-for-RNs.git

This project uses the commercial Gurobi solver. You must:

1. Install Gurobi from [https://www.gurobi.com](https://www.gurobi.com)
2. Set your `GUROBI_HOME` environment variable:
   ```bash
   export GUROBI_HOME="/path/to/gurobi"

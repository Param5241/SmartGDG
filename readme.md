**SmartCGD – Automated DXF Generation System**


**Overview:**
SmartCGD is a Python-based automation project designed to generate DXF drawings dynamically based on structured JSON input data and predefined DXF templates.
Each template has a dedicated processing script that reads input data, applies it to the corresponding DXF template, and produces a finalized DXF output.


**This architecture ensures:**
Clean separation of data, templates, scripts, and assets
Scalability for adding new templates
Maintainable and predictable execution flow



**Project Directory Structure:**

SmartCGD/
│
├── Data/
│   └── *.json
│       └── Input data files (template-specific JSON)
│
├── Template/
│   ├── template1.dxf
│   ├── template2.dxf
│   ├── template3.dxf
│   └── template4.dxf
│
├── icons/
│   └── *.png
│       └── PNG symbol/icon files used inside DXF drawings
│
├── scripts/
│   ├── one.py
│   ├── two.py
│   ├── three.py
│   └── four.py
│
└── README.md



**Folder Description**

**1. Data/**
Contains JSON input files
Each JSON file specifies:
  Template identifier (e.g., template1)
  Drawing data (lines, symbols, text, tables, metadata, etc.)
  Acts as the single source of truth for DXF generation

**2. Template/**
Contains DXF template files
Each template defines:
  Fixed layout
  Title blocks
  Legends
  Tables
  Drawing boundaries
Naming convention:
  template1.dxf
  template2.dxf
  template3.dxf
  ...

**3. icons/**
Contains PNG icon files
Used for:
  Symbols
  Fittings
  Legends
  Utility markers
Icons are inserted into DXF drawings by scripts

**4. scripts/**
Contains Python scripts, one per template
Each script:
  Loads the corresponding DXF template
  Reads JSON data
  Applies geometry, symbols, tables, and annotations
  Exports the final DXF file

  

**Workflow Explanation**

**Prepare Input Data**
  Place the JSON file inside the Data/ folder
  JSON must include which template it belongs to (e.g., template1)
**Select the Matching Script**
  If JSON is for template1 → run one.py
  If JSON is for template2 → run two.py
  And so on…
**Script Execution**
  Script reads:
    JSON data from Data/
    DXF template from Template/
    Icons from icons/
    
**Output Generation**
Final DXF file is generated
Output DXF is ready for use in AutoCAD or compatible CAD software

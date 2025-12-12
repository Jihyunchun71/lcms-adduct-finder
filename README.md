# LCMS Adduct Finder

**Automated Targeted Feature Extraction & Adduct Verification Tool for LC-MS Data.**

This Python tool is designed for the **targeted analysis** of LC-MS data. By providing a list of **Chemical Formulas**, it automatically performs a comprehensive scan for various adduct forms (e.g., `[M+H]+`, `[M+Na]+`). It extracts Extracted Ion Chromatograms (EIC) and rigorously evaluates **peak quality** using Gaussian fitting to determine the reliability of the detected signals.

---

## Key Features

* **Targeted Extraction**: Instantly converts chemical formulas (e.g., `C6H12O6`) into target m/z values, enabling precise extraction of specific metabolites or compounds.
* **Multi-Adduct Verification**:
    * Automatically scans for **16+ different adduct types** (Monomers, Dimers, Na/NH4 adducts, etc.) simultaneously.
    * Helps confirm the identity of a substance by checking if multiple adducts elute at the same Retention Time (RT).
* **Peak Quality & Existence Check**:
    * **Gaussian Scoring**: Fits a Gaussian curve to the raw peak data and calculates an $R^2$ score.
    * Distinguishes high-quality peaks ("Excellent/Good") from noise or irregular shapes ("Poor/Noise").
* **Precision Mass Calculation**: Uses high-precision logic considering electron mass:
    $$m/z = \frac{(M \times n + \Delta) - (Charge \times m_e)}{|Charge|}$$
* **Visual Inspection**: Automatically saves EIC plots (Raw vs. Gaussian Fit) as PNG images to visually verify the existence of the compound.

---

## Supported Adducts

The tool automatically detects the ionization mode (Positive/Negative) and scans for the following adducts to maximize detection coverage:

| Mode | Adduct Types |
| :--- | :--- |
| **Positive (+)** | `[M+H]+`, `[M+Na]+`, `[M+NH4]+`, `[M+ACN+H]+`, `[2M+H]+`, etc. |
| **Negative (-)** | `[M-H]-`, `[M+Cl]-`, `[M+HCOO]-`, `[M+FA-H]-`, `[2M-H]-`, etc. |

---

## Prerequisites

### 1. Data Conversion
This tool reads **`.mzML`** files.
If you have vendor formats (e.g., Thermo `.raw`), please convert them using **MSConvert (ProteoWizard)** first.

### 2. Python Dependencies
Install the required libraries:

```bash
pip install pyopenms pandas openpyxl molmass scipy matplotlib
```

---

## Usage

### Step 1. Prepare Data Folder
Place your `.mzML` files in a specific folder (e.g., `C:/Data/MyProject`).

### Step 2. Prepare Input Excel
Create a file named `file_list.xlsx` with a sheet named **`Final`**.
The structure should be as follows:

| RawFile | Mode | Formula 1 | Formula 2 | ... |
| :--- | :--- | :--- | :--- | :--- |
| `sample_01.raw` | POS | C6H12O6 | C10H16N5O13P3 | ... |
| `sample_02.raw` | NEG | C6H12O6 | | ... |

* **RawFile**: Filename (extension can be `.raw` or `.mzML`).
* **Mode**: `POS` or `NEG`.
* **Columns 3~**: List the chemical formulas you want to target in that file.

### Step 3. Configure & Run
Open the script (`main.py`) and set your data folder path:

```python
# In main.py
MZML_DATA_FOLDER = r"C:\Path\To\Your\mzML_Files"
```

Then run the script:

```bash
python main.py
```

---

## Output Examples

The tool generates reports that help you answer: *"Does my target compound exist in this sample, and is the signal reliable?"*

### 1. Excel Report (`Final_Result_With_Plots.xlsx`)
* **All_Features**: Raw data for every detected adduct.
* **[Formula_Name] Sheets**:
    * **Area Table**: Pivot table of peak areas for each adduct.
    * **Retention Time Table**: Check if RTs are consistent across different adducts.

### 2. EIC Plots (`EIC_Plots_Export/`)
Visual validation of the detected peaks.

* **Blue Line**: Raw EIC Data
* **Red Dashed Line**: Gaussian Fit Curve (High overlap indicates a high-quality peak)


---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

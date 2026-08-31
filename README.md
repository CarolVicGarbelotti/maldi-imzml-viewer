# Standalone MALDI imzML Viewer

A Windows-friendly desktop viewer for SCiLS-exported, RMS-normalized imzML datasets.

## Installation

Open PowerShell in this folder and create a dedicated environment:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If PowerShell blocks activation, the environment can be used without activation:

```powershell
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe maldi_viewer.py
```

After installation, start with `run_viewer.bat`. It automatically uses the local `.venv` when present.

## Workflow

1. Click **Open imzML** and choose the `.imzML` file. Keep its `.ibd` file beside it.
2. Choose an average-spectrum bin width. `0.01 Da` is a practical starting point.
3. Hover over the spectrum to snap to the nearest local peak and show its m/z and intensity.
4. Click a peak to reconstruct an ion image.
5. Hover over the ion image for the exact array value at a pixel.
6. Click an acquired pixel to show its spectrum.
7. Move/resize the cyan rectangle and click **Average rectangle ROI** for an ROI spectrum.
8. Click **Whole-image mean** to restore the total mean spectrum.

## Normalization

The primary display is labeled **Stored intensity (SCiLS RMS-normalized)** because these imzML files were exported from SCiLS after RMS normalization. It is not original detector-count intensity.

**Additional TIC normalization** divides the already RMS-normalized selected-ion value by the sum of all stored values in that pixel. This is a second normalization and should only be used intentionally.

## Publication export

- Ion images: TIFF, PNG, PDF, or SVG.
- Spectra: SVG, PDF, TIFF, or PNG.
- Default raster resolution: 600 dpi.
- Enter the dataset calibration in µm/pixel to add a physical scale bar. Leave it at zero to omit the scale bar.

The export is redrawn with Matplotlib; it is not a screenshot of the interactive window.

## Performance notes

- The mean spectrum and each new ion image require one pass through the imzML spectra.
- Previously calculated `(m/z, tolerance)` images are cached in RAM.
- Changing between stored RMS intensity and additional TIC normalization is immediate because both images are calculated together.
- The cache is cleared whenever a new dataset is opened and disappears when the application closes.
- While a spectrum or image calculation is running, additional plot clicks are temporarily ignored. The status line changes back to **Ready** when another peak can be selected.

# Standalone MALDI imzML Viewer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24447064.svg)](https://doi.org/10.5281/zenodo.24447064)
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

After installation, start with `run_viewer.bat` if `python` points to this environment, or use the explicit `.venv` command above.

## Workflow

1. Click **Open imzML** and choose the `.imzML` file. Keep its `.ibd` file beside it.
2. Choose an average-spectrum bin width. `0.01 Da` is a practical starting point.
3. Hover over the spectrum to snap to the nearest local peak and show its m/z and intensity.
4. Click a peak to reconstruct an ion image.
5. Hover over the ion image for the exact array value at a pixel.
6. Click an acquired pixel to show its spectrum.
7. Move/resize the cyan rectangle and click **Average rectangle ROI** for an ROI spectrum.
8. Click **Whole-image mean spectrum** to restore the spectrum averaged across all acquired pixels.
9. Click **Total ion image** to display the sum of the complete stored spectrum at each pixel.

The total ion image uses the intensities stored in the imzML. For SCiLS RMS-normalized exports, it is therefore the sum of the stored RMS-normalized spectrum, not original detector-count TIC.

## Multi-ion overlays

1. Calculate an ion image by clicking its spectral peak.
2. Choose an overlay color and click **Add current ion**.
3. Repeat for each ion, then choose a scaling mode and click **Build overlay**.
4. Use the interactive ion table to show/hide ions, change colors, or remove individual ions.
5. Use **Clear overlay** only when you want to begin a completely new combination.

When an overlay is already displayed, table edits update it immediately. Hidden ions remain in the list and can be restored without recalculating their ion images.

Available colors are magenta, green, cyan, red, blue, yellow, orange, and white.

- **Separate per ion** scales every ion to its own selected upper percentile. This is useful for comparing spatial distributions and colocalization, but color brightness does not represent abundance differences between ions.
- **Shared across ions** uses one percentile-derived denominator for every ion. This preserves relative intensity differences when the ion values are scientifically comparable.

Overlay exports include an ion-color legend and support the same full-image/ROI, pixel-axis, scale-bar color, and scale-bar position controls as single-ion exports.

## Normalization

The primary display is labeled **Stored intensity (SCiLS RMS-normalized)** because these imzML files were exported from SCiLS after RMS normalization. It is not original detector-count intensity.

**Additional TIC normalization** divides the already RMS-normalized selected-ion value by the sum of all stored values in that pixel. This is a second normalization and should only be used intentionally.

## Publication export

- Ion images: TIFF, PNG, PDF, or SVG.
- Total ion images and multi-ion overlays use the same export formats and controls.
- Spectra: SVG, PDF, TIFF, or PNG.
- Default raster resolution: 600 dpi.
- Export either the full image or the current rectangular ROI.
- Pixel axes are hidden by default and can be restored with **Show pixel axes**.
- Choose a white or black scale bar and place it in any corner.
- Enter the dataset calibration in µm/pixel to add a physical scale bar. Leave it at zero to omit the scale bar.

The export is redrawn with Matplotlib; it is not a screenshot of the interactive window.

## Performance notes

- The mean spectrum and total ion image are calculated together during the initial pass through the imzML spectra.
- Each new ion image requires one additional pass through the spectra.
- The cached total ion image is reused as the denominator for optional TIC normalization.
- Previously calculated `(m/z, tolerance)` images are cached in RAM.
- Changing between stored RMS intensity and additional TIC normalization is immediate because both images are calculated together.
- The cache is cleared whenever a new dataset is opened and disappears when the application closes.

## Author and affiliation

**Carolina Garbelotti**  
Universidade de São Paulo  
ORCID: [0000-0003-2631-4473](https://orcid.org/0000-0003-2631-4473)

## Citation

If you use this software, please cite:

> Garbelotti, C. (2026). *MALDI imzML Viewer* (Version 0.1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24447064

The repository includes a [CITATION.cff](CITATION.cff) file so GitHub can generate APA and BibTeX citations. The archived version is available through [Zenodo](https://doi.org/10.5281/zenodo.24447064).

## Scope and limitations

- The viewer is intended for exploratory visualization and figure preparation; it does not replace validation in vendor software or a validated quantitative pipeline.
- Intensities are interpreted exactly as stored in the imzML file. Input exported with RMS or another normalization remains normalized in the viewer.
- The total-ion image is the per-pixel sum of stored intensities and is not original detector-count TIC when the source data were already normalized.
- **Additional TIC normalization** is a second normalization applied to the stored data and should be used intentionally.
- Ion images depend on the selected m/z and tolerance; the viewer does not perform automatic mass recalibration, peak identification or deconvolution.
- Physical scale bars require the user to provide the correct spatial calibration in µm/pixel.
- Performance depends on dataset size, storage speed and available system memory.
- The viewer reads imzML/ibd pairs and does not convert proprietary vendor formats.

## Development transparency

The software was developed through an AI-assisted programming workflow using OpenAI ChatGPT/Codex. Scientific requirements, domain-specific decisions, testing with experimental data and validation of outputs were conducted by the researcher.

## License

This project is distributed under the [MIT License](LICENSE).

# Changelog

All notable changes to MALDI imzML Viewer are documented in this file.

## [0.1.0] - 2026-11-22

### Added

- Loading of paired imzML and ibd datasets.
- Memory-conscious calculation of the whole-image mean spectrum.
- Interactive spectrum hover and peak selection.
- Reconstruction and caching of selected-ion images.
- Total-ion image calculated from the intensities stored in the imzML file.
- Pixel-spectrum inspection and rectangular ROI mean spectra.
- Stored-intensity and optional additional TIC-normalized displays.
- Custom multi-ion overlays with editable colors, visibility, scaling and removal.
- Full-image and ROI-only publication export.
- Configurable DPI, pixel axes and physical scale bar.
- TIFF, PNG, PDF and SVG image export.
- Windows launcher using a local Python virtual environment.

### Scientific note

For SCiLS-exported RMS-normalized imzML datasets, displayed stored intensities
and the total-ion image are derived from the RMS-normalized values in the file,
not from original detector counts.

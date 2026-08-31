from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pyqtgraph as pg
from pyimzml.ImzMLParser import ImzMLParser
from PySide6 import QtCore, QtGui, QtWidgets


pg.setConfigOption("imageAxisOrder", "row-major")
pg.setConfigOptions(antialias=False)


OVERLAY_COLORS = {
    "Magenta": (1.0, 0.0, 1.0),
    "Green": (0.0, 1.0, 0.0),
    "Cyan": (0.0, 1.0, 1.0),
    "Red": (1.0, 0.0, 0.0),
    "Blue": (0.0, 0.35, 1.0),
    "Yellow": (1.0, 1.0, 0.0),
    "Orange": (1.0, 0.45, 0.0),
    "White": (1.0, 1.0, 1.0),
}


@dataclass
class Dataset:
    path: Path
    parser: ImzMLParser
    coordinates: np.ndarray
    coordinate_to_index: dict[tuple[int, int], int]
    x_min: int
    y_min: int
    shape: tuple[int, int]
    average_mz: np.ndarray
    average_intensity: np.ndarray
    bin_width: float


class WorkerSignals(QtCore.QObject):
    result = QtCore.Signal(object)
    error = QtCore.Signal(str)
    finished = QtCore.Signal()


class Worker(QtCore.QRunnable):
    def __init__(self, function, *args):
        super().__init__()
        self.function = function
        self.args = args
        self.signals = WorkerSignals()

    @QtCore.Slot()
    def run(self):
        try:
            self.signals.result.emit(self.function(*self.args))
        except Exception as exc:  # surfaced in the GUI
            self.signals.error.emit(f"{type(exc).__name__}: {exc}")
        finally:
            self.signals.finished.emit()


class MaldiViewer(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MALDI imzML Viewer — SCiLS RMS-aware")
        self.resize(1450, 950)
        self.dataset: Dataset | None = None
        self.thread_pool = QtCore.QThreadPool.globalInstance()
        self.thread_pool.setMaxThreadCount(1)
        self._busy = False
        # Keep QRunnable objects alive until their queued Qt signals finish.
        self._active_workers: set[Worker] = set()
        self.ion_cache: dict[tuple[float, float], tuple[np.ndarray, np.ndarray]] = {}
        self.current_raw_image: np.ndarray | None = None
        self.current_tic_image: np.ndarray | None = None
        self.current_display_image: np.ndarray | None = None
        self.current_mz: float | None = None
        self.current_spectrum_mz = np.array([])
        self.current_spectrum_intensity = np.array([])
        self.current_spectrum_title = "No dataset loaded"
        self.overlay_layers: list[dict] = []
        self.current_overlay_rgb: np.ndarray | None = None
        self.display_is_overlay = False
        self._build_ui()

    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        outer = QtWidgets.QVBoxLayout(central)

        controls = QtWidgets.QHBoxLayout()
        self.open_button = QtWidgets.QPushButton("Open imzML")
        self.total_button = QtWidgets.QPushButton("Whole-image mean")
        self.roi_button = QtWidgets.QPushButton("Average rectangle ROI")
        self.export_image_button = QtWidgets.QPushButton("Export ion image")
        self.export_spectrum_button = QtWidgets.QPushButton("Export spectrum")
        for button in (self.total_button, self.roi_button, self.export_image_button, self.export_spectrum_button):
            button.setEnabled(False)
        controls.addWidget(self.open_button)
        controls.addWidget(self.total_button)
        controls.addWidget(self.roi_button)
        controls.addStretch(1)

        controls.addWidget(QtWidgets.QLabel("± m/z:"))
        self.tolerance = QtWidgets.QDoubleSpinBox()
        self.tolerance.setDecimals(4)
        self.tolerance.setRange(0.0001, 10.0)
        self.tolerance.setValue(0.05)
        self.tolerance.setSingleStep(0.01)
        controls.addWidget(self.tolerance)

        self.normalization = QtWidgets.QComboBox()
        self.normalization.addItems([
            "Stored intensity (SCiLS RMS-normalized)",
            "Additional TIC normalization",
        ])
        controls.addWidget(self.normalization)
        controls.addWidget(self.export_image_button)
        controls.addWidget(self.export_spectrum_button)
        outer.addLayout(controls)

        export_controls = QtWidgets.QHBoxLayout()
        export_controls.addWidget(QtWidgets.QLabel("Export DPI:"))
        self.dpi = QtWidgets.QSpinBox()
        self.dpi.setRange(72, 1200)
        self.dpi.setValue(600)
        export_controls.addWidget(self.dpi)
        export_controls.addWidget(QtWidgets.QLabel("µm/pixel (0 = omit scale bar):"))
        self.microns_per_pixel = QtWidgets.QDoubleSpinBox()
        self.microns_per_pixel.setDecimals(5)
        self.microns_per_pixel.setRange(0, 10000)
        self.microns_per_pixel.setValue(0)
        export_controls.addWidget(self.microns_per_pixel)
        export_controls.addWidget(QtWidgets.QLabel("Scale bar (µm):"))
        self.scale_bar_um = QtWidgets.QDoubleSpinBox()
        self.scale_bar_um.setRange(1, 100000)
        self.scale_bar_um.setValue(50)
        export_controls.addWidget(self.scale_bar_um)

        export_controls.addWidget(QtWidgets.QLabel("Export area:"))
        self.export_area = QtWidgets.QComboBox()
        self.export_area.addItems(["Full image", "Rectangle ROI"])
        export_controls.addWidget(self.export_area)

        self.show_pixel_axes = QtWidgets.QCheckBox("Show pixel axes")
        self.show_pixel_axes.setChecked(False)
        export_controls.addWidget(self.show_pixel_axes)

        export_controls.addWidget(QtWidgets.QLabel("Scale-bar color:"))
        self.scale_bar_color = QtWidgets.QComboBox()
        self.scale_bar_color.addItems(["White", "Black"])
        export_controls.addWidget(self.scale_bar_color)

        export_controls.addWidget(QtWidgets.QLabel("Position:"))
        self.scale_bar_position = QtWidgets.QComboBox()
        self.scale_bar_position.addItems([
            "Bottom left", "Bottom right", "Top left", "Top right"
        ])
        export_controls.addWidget(self.scale_bar_position)
        export_controls.addStretch(1)
        outer.addLayout(export_controls)

        overlay_controls = QtWidgets.QHBoxLayout()
        overlay_controls.addWidget(QtWidgets.QLabel("Overlay color:"))
        self.overlay_color = QtWidgets.QComboBox()
        self.overlay_color.addItems(list(OVERLAY_COLORS))
        overlay_controls.addWidget(self.overlay_color)
        self.add_overlay_button = QtWidgets.QPushButton("Add current ion")
        self.add_overlay_button.setEnabled(False)
        overlay_controls.addWidget(self.add_overlay_button)
        overlay_controls.addWidget(QtWidgets.QLabel("Scaling:"))
        self.overlay_scaling = QtWidgets.QComboBox()
        self.overlay_scaling.addItems(["Separate per ion", "Shared across ions"])
        overlay_controls.addWidget(self.overlay_scaling)
        overlay_controls.addWidget(QtWidgets.QLabel("Upper percentile:"))
        self.overlay_percentile = QtWidgets.QDoubleSpinBox()
        self.overlay_percentile.setDecimals(1)
        self.overlay_percentile.setRange(90.0, 100.0)
        self.overlay_percentile.setValue(99.5)
        self.overlay_percentile.setSingleStep(0.5)
        overlay_controls.addWidget(self.overlay_percentile)
        self.build_overlay_button = QtWidgets.QPushButton("Build overlay")
        self.build_overlay_button.setEnabled(False)
        self.clear_overlay_button = QtWidgets.QPushButton("Clear overlay")
        self.clear_overlay_button.setEnabled(False)
        overlay_controls.addWidget(self.build_overlay_button)
        overlay_controls.addWidget(self.clear_overlay_button)
        self.overlay_summary = QtWidgets.QLabel("No overlay ions selected")
        overlay_controls.addWidget(self.overlay_summary, 1)
        outer.addLayout(overlay_controls)

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.spectrum_plot = pg.PlotWidget(background="w")
        self.spectrum_plot.setLabel("bottom", "m/z")
        self.spectrum_plot.setLabel("left", "Intensity")
        self.spectrum_curve = self.spectrum_plot.plot(pen=pg.mkPen("k", width=1))
        self.spectrum_crosshair = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen("#d62728", width=1))
        self.spectrum_plot.addItem(self.spectrum_crosshair, ignoreBounds=True)
        self.spectrum_crosshair.hide()
        self.spectrum_hover = pg.TextItem(anchor=(0, 1), fill=pg.mkBrush(255, 255, 255, 220), border=pg.mkPen("#555"))
        self.spectrum_plot.addItem(self.spectrum_hover)
        self.spectrum_hover.hide()

        image_plot = pg.PlotItem()
        image_plot.invertY(True)
        image_plot.setAspectLocked(True)
        image_plot.setLabel("bottom", "x pixel")
        image_plot.setLabel("left", "y pixel")
        self.image_view = pg.ImageView(view=image_plot)
        self.image_item = self.image_view.getImageItem()
        self.roi = pg.RectROI([5, 5], [20, 20], pen=pg.mkPen("cyan", width=2))
        self.roi.addScaleHandle([1, 1], [0, 0])
        self.roi.addScaleHandle([0, 0], [1, 1])
        image_plot.addItem(self.roi)
        self.roi.hide()
        self.image_hover = pg.TextItem(anchor=(0, 1), fill=pg.mkBrush(255, 255, 255, 220), border=pg.mkPen("#555"))
        image_plot.addItem(self.image_hover)
        self.image_hover.hide()

        splitter.addWidget(self.spectrum_plot)
        splitter.addWidget(self.image_view)
        splitter.setSizes([400, 550])
        outer.addWidget(splitter)

        self.status = QtWidgets.QLabel("Open an imzML file. The matching .ibd must be in the same folder.")
        self.status.setWordWrap(True)
        outer.addWidget(self.status)

        self.open_button.clicked.connect(self.open_dataset)
        self.total_button.clicked.connect(self.show_total_spectrum)
        self.roi_button.clicked.connect(self.average_current_roi)
        self.export_image_button.clicked.connect(self.export_image)
        self.export_spectrum_button.clicked.connect(self.export_spectrum)
        self.add_overlay_button.clicked.connect(self.add_current_ion_to_overlay)
        self.build_overlay_button.clicked.connect(self.build_overlay)
        self.clear_overlay_button.clicked.connect(self.clear_overlay)
        self.normalization.currentIndexChanged.connect(self.refresh_image)
        self.tolerance.valueChanged.connect(self.recalculate_current_mz)
        self.spectrum_plot.scene().sigMouseMoved.connect(self.on_spectrum_hover)
        self.spectrum_plot.scene().sigMouseClicked.connect(self.on_spectrum_click)
        image_plot.scene().sigMouseMoved.connect(self.on_image_hover)
        image_plot.scene().sigMouseClicked.connect(self.on_image_click)

    def open_dataset(self):
        if self._busy:
            return
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open imzML", "", "imzML (*.imzML)")
        if not filename:
            return
        path = Path(filename)
        if not path.with_suffix(".ibd").exists():
            QtWidgets.QMessageBox.critical(self, "Missing .ibd", "The matching .ibd file was not found.")
            return
        bin_width, ok = QtWidgets.QInputDialog.getDouble(
            self, "Average-spectrum binning", "Bin width (Da):", 0.01, 0.0001, 1.0, 4
        )
        if not ok:
            return
        self.set_busy(True, "Reading dataset and calculating the whole-image mean spectrum...")
        worker = Worker(self._load_dataset, path, bin_width)
        worker.signals.result.connect(self._dataset_loaded)
        worker.signals.error.connect(self.show_error)
        worker.signals.finished.connect(lambda w=worker: self._worker_finished(w))
        self._active_workers.add(worker)
        self.thread_pool.start(worker)

    @staticmethod
    def _load_dataset(path: Path, bin_width: float) -> Dataset:
        parser = ImzMLParser(str(path))
        coords = np.asarray(parser.coordinates, dtype=int)
        x_min, y_min = int(coords[:, 0].min()), int(coords[:, 1].min())
        shape = (int(coords[:, 1].max()) - y_min + 1, int(coords[:, 0].max()) - x_min + 1)
        first_mz, _ = parser.getspectrum(0)
        low = np.floor(first_mz.min() / bin_width) * bin_width
        high = np.ceil(first_mz.max() / bin_width) * bin_width
        edges = np.arange(low, high + bin_width, bin_width)
        summed = np.zeros(edges.size - 1, dtype=np.float64)
        for i in range(len(coords)):
            mz, intensity = parser.getspectrum(i)
            valid = (mz >= edges[0]) & (mz < edges[-1])
            bins = np.floor((mz[valid] - edges[0]) / bin_width).astype(np.int64)
            np.add.at(summed, bins, intensity[valid])
        return Dataset(
            path, parser, coords,
            {(int(c[0]), int(c[1])): i for i, c in enumerate(coords)},
            x_min, y_min, shape, edges[:-1] + bin_width / 2,
            summed / len(coords), bin_width,
        )

    def _dataset_loaded(self, dataset: Dataset):
        self.dataset = dataset
        self.ion_cache.clear()
        self.current_raw_image = self.current_tic_image = self.current_display_image = None
        self.clear_overlay()
        self.show_spectrum(dataset.average_mz, dataset.average_intensity, "Whole-image mean spectrum — stored SCiLS RMS-normalized values")
        self.setWindowTitle(f"MALDI Viewer — {dataset.path.name}")
        self.status.setText(f"Loaded {len(dataset.coordinates):,} pixels; image grid {dataset.shape[1]} × {dataset.shape[0]}. Click a spectral peak.")
        for button in (self.total_button, self.roi_button, self.export_spectrum_button):
            button.setEnabled(True)

    def show_spectrum(self, mz, intensity, title: str):
        self.current_spectrum_mz = np.asarray(mz)
        self.current_spectrum_intensity = np.asarray(intensity)
        self.current_spectrum_title = title
        self.spectrum_curve.setData(self.current_spectrum_mz, self.current_spectrum_intensity)
        self.spectrum_plot.setTitle(title)
        self.spectrum_plot.enableAutoRange()

    def show_total_spectrum(self):
        if self.dataset:
            self.show_spectrum(self.dataset.average_mz, self.dataset.average_intensity, "Whole-image mean spectrum — stored SCiLS RMS-normalized values")

    def nearest_local_peak(self, x: float) -> int:
        mz, intensity = self.current_spectrum_mz, self.current_spectrum_intensity
        i = int(np.clip(np.searchsorted(mz, x), 0, len(mz) - 1))
        radius = max(2, int(0.05 / max(np.median(np.diff(mz)), 1e-12)))
        left, right = max(0, i - radius), min(len(mz), i + radius + 1)
        return left + int(np.argmax(intensity[left:right]))

    def on_spectrum_hover(self, scene_pos):
        if not len(self.current_spectrum_mz) or not self.spectrum_plot.sceneBoundingRect().contains(scene_pos):
            self.spectrum_hover.hide()
            return
        point = self.spectrum_plot.plotItem.vb.mapSceneToView(scene_pos)
        i = self.nearest_local_peak(point.x())
        mz, value = self.current_spectrum_mz[i], self.current_spectrum_intensity[i]
        self.spectrum_crosshair.setPos(mz)
        self.spectrum_crosshair.show()
        self.spectrum_hover.setText(f"m/z {mz:.4f}\nIntensity {value:.5g}", color="k")
        self.spectrum_hover.setPos(mz, value)
        self.spectrum_hover.show()

    def on_spectrum_click(self, event):
        if self._busy or not self.dataset or not self.spectrum_plot.sceneBoundingRect().contains(event.scenePos()):
            return
        x = self.spectrum_plot.plotItem.vb.mapSceneToView(event.scenePos()).x()
        i = self.nearest_local_peak(x)
        self.calculate_ion_image(float(self.current_spectrum_mz[i]))

    def calculate_ion_image(self, target_mz: float):
        if self._busy or not self.dataset:
            return
        self.current_mz = target_mz
        tolerance = self.tolerance.value()
        key = (round(target_mz, 6), round(tolerance, 6))
        if key in self.ion_cache:
            self.receive_ion_images(self.ion_cache[key])
            return
        self.set_busy(True, f"Calculating m/z {target_mz:.4f} ± {tolerance:.4f} Da...")
        worker = Worker(self._reconstruct_images, self.dataset, target_mz, tolerance)
        worker.signals.result.connect(lambda result: self._cache_and_receive(key, result))
        worker.signals.error.connect(self.show_error)
        worker.signals.finished.connect(lambda w=worker: self._worker_finished(w))
        self._active_workers.add(worker)
        self.thread_pool.start(worker)

    @staticmethod
    def _reconstruct_images(dataset: Dataset, target_mz: float, tolerance: float):
        raw = np.full(dataset.shape, np.nan, dtype=np.float64)
        tic_normalized = np.full(dataset.shape, np.nan, dtype=np.float64)
        low, high = target_mz - tolerance, target_mz + tolerance
        for i, coordinate in enumerate(dataset.coordinates):
            mz, intensity = dataset.parser.getspectrum(i)
            left = np.searchsorted(mz, low, side="left")
            right = np.searchsorted(mz, high, side="right")
            value = float(np.sum(intensity[left:right], dtype=np.float64))
            tic = float(np.sum(intensity, dtype=np.float64))
            row, col = int(coordinate[1]) - dataset.y_min, int(coordinate[0]) - dataset.x_min
            raw[row, col] = value
            tic_normalized[row, col] = value / tic if tic > 0 else np.nan
        return raw, tic_normalized

    def _cache_and_receive(self, key, result):
        self.ion_cache[key] = result
        self.receive_ion_images(result)

    def receive_ion_images(self, result):
        self.current_raw_image, self.current_tic_image = result
        self.refresh_image()
        self.roi.show()
        self.export_image_button.setEnabled(True)
        self.add_overlay_button.setEnabled(True)
        self.status.setText(
            f"Ready — displaying m/z {self.current_mz:.4f} ± {self.tolerance.value():.4f} Da. "
            "Click another spectral peak to calculate a new image."
        )

    def refresh_image(self):
        if self.current_raw_image is None:
            return
        self.display_is_overlay = False
        self.current_display_image = self.current_raw_image if self.normalization.currentIndex() == 0 else self.current_tic_image
        self.image_view.setImage(self.current_display_image, autoRange=True, autoLevels=True)
        label = self.normalization.currentText()
        self.image_view.getView().setTitle(f"m/z {self.current_mz:.4f} ± {self.tolerance.value():.4f} Da — {label}")

    def recalculate_current_mz(self):
        if not self._busy and self.current_mz is not None:
            self.calculate_ion_image(self.current_mz)

    def on_image_hover(self, scene_pos):
        if self.current_display_image is None or not self.image_view.getView().sceneBoundingRect().contains(scene_pos):
            self.image_hover.hide()
            return
        point = self.image_view.getView().vb.mapSceneToView(scene_pos)
        col, row = int(np.floor(point.x())), int(np.floor(point.y()))
        if 0 <= row < self.current_display_image.shape[0] and 0 <= col < self.current_display_image.shape[1]:
            x = col + self.dataset.x_min
            y = row + self.dataset.y_min
            if self.display_is_overlay:
                lines = [f"x={x}, y={y}"]
                image_key = "raw" if self.normalization.currentIndex() == 0 else "tic"
                for layer in self.overlay_layers:
                    lines.append(
                        f"{layer['mz']:.4f} ({layer['color']}): "
                        f"{layer[image_key][row, col]:.5g}"
                    )
                hover_text = "\n".join(lines)
            else:
                value = self.current_display_image[row, col]
                hover_text = f"x={x}, y={y}\nIntensity={value:.6g}"
            self.image_hover.setText(hover_text, color="k")
            self.image_hover.setPos(col, row)
            self.image_hover.show()

    def on_image_click(self, event):
        if self._busy or not self.dataset or self.current_display_image is None:
            return
        point = self.image_view.getView().vb.mapSceneToView(event.scenePos())
        x = int(np.floor(point.x())) + self.dataset.x_min
        y = int(np.floor(point.y())) + self.dataset.y_min
        index = self.dataset.coordinate_to_index.get((x, y))
        if index is None:
            return
        mz, intensity = self.dataset.parser.getspectrum(index)
        self.show_spectrum(mz, intensity, f"Pixel spectrum: x={x}, y={y} — stored SCiLS RMS-normalized values")

    def average_current_roi(self):
        if not self.dataset or self.current_display_image is None:
            return
        x0, x1, y0, y1 = self._current_roi_float_bounds()
        indices = [
            i for i, c in enumerate(self.dataset.coordinates)
            if x0 <= int(c[0]) - self.dataset.x_min < x1 and y0 <= int(c[1]) - self.dataset.y_min < y1
        ]
        if not indices:
            self.status.setText("The ROI contains no acquired pixels.")
            return
        self.set_busy(True, f"Averaging {len(indices):,} ROI spectra...")
        worker = Worker(self._average_indices, self.dataset, indices)
        worker.signals.result.connect(lambda result: self.show_spectrum(self.dataset.average_mz, result, f"ROI mean spectrum ({len(indices):,} pixels)"))
        worker.signals.error.connect(self.show_error)
        worker.signals.finished.connect(lambda w=worker: self._worker_finished(w))
        self._active_workers.add(worker)
        self.thread_pool.start(worker)

    @staticmethod
    def _average_indices(dataset: Dataset, indices: list[int]):
        accumulated = np.zeros_like(dataset.average_intensity)
        left_edge = dataset.average_mz[0] - dataset.bin_width / 2
        for i in indices:
            mz, intensity = dataset.parser.getspectrum(i)
            bins = np.floor((mz - left_edge) / dataset.bin_width).astype(np.int64)
            valid = (bins >= 0) & (bins < accumulated.size)
            np.add.at(accumulated, bins[valid], intensity[valid])
        return accumulated / len(indices)

    def _current_roi_float_bounds(self):
        pos, size = self.roi.pos(), self.roi.size()
        x0, x1 = sorted((float(pos.x()), float(pos.x() + size.x())))
        y0, y1 = sorted((float(pos.y()), float(pos.y() + size.y())))
        return x0, x1, y0, y1

    def _current_roi_array_bounds(self, image: np.ndarray):
        x0, x1, y0, y1 = self._current_roi_float_bounds()
        col0 = max(0, int(np.floor(x0)))
        col1 = min(image.shape[1], int(np.ceil(x1)))
        row0 = max(0, int(np.floor(y0)))
        row1 = min(image.shape[0], int(np.ceil(y1)))
        if row1 <= row0 or col1 <= col0:
            return None
        return row0, row1, col0, col1

    def add_current_ion_to_overlay(self):
        if self.current_mz is None or self.current_raw_image is None:
            self.status.setText("Calculate an ion image before adding it to an overlay.")
            return
        color = self.overlay_color.currentText()
        tolerance = float(self.tolerance.value())
        for layer in self.overlay_layers:
            if np.isclose(layer["mz"], self.current_mz) and np.isclose(layer["tolerance"], tolerance):
                layer["color"] = color
                self._update_overlay_summary()
                self.status.setText(f"Updated m/z {self.current_mz:.4f} to {color} in the overlay.")
                return
        self.overlay_layers.append({
            "mz": float(self.current_mz),
            "tolerance": tolerance,
            "color": color,
            "raw": self.current_raw_image.copy(),
            "tic": self.current_tic_image.copy(),
        })
        self._update_overlay_summary()
        self.build_overlay_button.setEnabled(True)
        self.clear_overlay_button.setEnabled(True)
        self.status.setText(f"Added m/z {self.current_mz:.4f} as {color}. Add another ion or build the overlay.")

    def _update_overlay_summary(self):
        if not self.overlay_layers:
            self.overlay_summary.setText("No overlay ions selected")
            return
        self.overlay_summary.setText("; ".join(
            f"{layer['mz']:.4f}={layer['color']}" for layer in self.overlay_layers
        ))

    def clear_overlay(self):
        was_displayed = self.display_is_overlay
        self.overlay_layers.clear()
        self.current_overlay_rgb = None
        self.display_is_overlay = False
        self._update_overlay_summary()
        self.build_overlay_button.setEnabled(False)
        self.clear_overlay_button.setEnabled(False)
        if was_displayed and self.current_raw_image is not None:
            self.refresh_image()

    def build_overlay(self):
        if not self.overlay_layers:
            self.status.setText("Add at least one ion before building an overlay.")
            return
        image_key = "raw" if self.normalization.currentIndex() == 0 else "tic"
        images = [np.asarray(layer[image_key], dtype=np.float64) for layer in self.overlay_layers]
        percentile = float(self.overlay_percentile.value())
        if self.overlay_scaling.currentText() == "Separate per ion":
            denominators = []
            for image in images:
                finite = image[np.isfinite(image)]
                denominators.append(float(np.percentile(finite, percentile)) if finite.size else 0.0)
        else:
            finite_parts = [image[np.isfinite(image)] for image in images]
            finite_parts = [part for part in finite_parts if part.size]
            shared = float(np.percentile(np.concatenate(finite_parts), percentile)) if finite_parts else 0.0
            denominators = [shared] * len(images)

        rgb = np.zeros((*images[0].shape, 3), dtype=np.float64)
        for layer, image, denominator in zip(self.overlay_layers, images, denominators):
            if denominator <= 0:
                continue
            scaled = np.nan_to_num(image, nan=0.0, posinf=denominator, neginf=0.0)
            scaled = np.clip(scaled / denominator, 0.0, 1.0)
            color = np.asarray(OVERLAY_COLORS[layer["color"]], dtype=np.float64)
            rgb += scaled[..., np.newaxis] * color
        self.current_overlay_rgb = np.clip(rgb, 0.0, 1.0)
        self.display_is_overlay = True
        self.image_view.setImage(
            np.round(self.current_overlay_rgb * 255).astype(np.uint8),
            autoRange=True,
            levels=(0, 255),
        )
        self.image_view.getView().setTitle(
            f"Ion overlay — {self.overlay_scaling.currentText()}, {percentile:g}th percentile — "
            f"{self.normalization.currentText()}"
        )
        self.roi.show()
        self.export_image_button.setEnabled(True)
        self.status.setText(f"Overlay ready with {len(self.overlay_layers)} ion(s).")

    def export_image(self):
        if self.current_display_image is None:
            return
        use_roi = self.export_area.currentText() == "Rectangle ROI"
        area_suffix = "_roi" if use_roi else ""
        is_overlay = self.display_is_overlay and self.current_overlay_rgb is not None
        base_name = "ion_overlay" if is_overlay else f"ion_mz_{self.current_mz:.4f}"
        filename, selected = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export publication ion image", f"{base_name}{area_suffix}.tiff",
            "TIFF (*.tif *.tiff);;PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        if not filename:
            return
        image = self.current_overlay_rgb if is_overlay else self.current_display_image
        if use_roi:
            bounds = self._current_roi_array_bounds(image)
            if bounds is None:
                QtWidgets.QMessageBox.warning(self, "Invalid ROI", "The rectangle does not overlap the ion image.")
                return
            row0, row1, col0, col1 = bounds
            image = image[row0:row1, col0:col1]

        aspect = image.shape[1] / max(image.shape[0], 1)
        figure_width = float(np.clip(4.5 * aspect + 1.2, 4.8, 9.0))
        fig, ax = plt.subplots(figsize=(figure_width, 5), constrained_layout=True)
        if is_overlay:
            ax.imshow(image, origin="upper", interpolation="nearest")
        else:
            shown = ax.imshow(image, origin="upper", cmap="inferno", interpolation="nearest")
        if self.show_pixel_axes.isChecked():
            ax.set_xlabel("x pixel")
            ax.set_ylabel("y pixel")
        else:
            ax.set_axis_off()
        if is_overlay:
            ax.set_title(
                f"Ion overlay — {self.overlay_scaling.currentText()}, "
                f"{self.overlay_percentile.value():g}th percentile"
            )
            legend_handles = [
                Patch(color=OVERLAY_COLORS[layer["color"]], label=f"m/z {layer['mz']:.4f}")
                for layer in self.overlay_layers
            ]
            ax.legend(
                handles=legend_handles,
                title=self.normalization.currentText(),
                loc="upper center",
                bbox_to_anchor=(0.5, -0.02),
                ncol=min(3, len(legend_handles)),
                frameon=False,
            )
        else:
            ax.set_title(f"m/z {self.current_mz:.4f} ± {self.tolerance.value():.4f} Da")
            cbar = fig.colorbar(shown, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(self.normalization.currentText())
        um_per_pixel = self.microns_per_pixel.value()
        if um_per_pixel > 0:
            pixels = self.scale_bar_um.value() / um_per_pixel
            if pixels >= image.shape[1] * 0.88:
                plt.close(fig)
                QtWidgets.QMessageBox.warning(
                    self,
                    "Scale bar too long",
                    "The requested scale bar is too long for this export area. "
                    "Choose a shorter bar or export a larger area.",
                )
                return
            margin_x = image.shape[1] * 0.06
            margin_y = image.shape[0] * 0.07
            position = self.scale_bar_position.currentText()
            if "left" in position.lower():
                bar_x0, bar_x1 = margin_x, margin_x + pixels
            else:
                bar_x1 = image.shape[1] - margin_x
                bar_x0 = bar_x1 - pixels
            if "bottom" in position.lower():
                bar_y = image.shape[0] - margin_y
                text_y, vertical_alignment = bar_y - image.shape[0] * 0.025, "bottom"
            else:
                bar_y = margin_y
                text_y, vertical_alignment = bar_y + image.shape[0] * 0.025, "top"
            color = self.scale_bar_color.currentText().lower()
            ax.plot([bar_x0, bar_x1], [bar_y, bar_y], color=color, linewidth=4, solid_capstyle="butt")
            ax.text(
                (bar_x0 + bar_x1) / 2,
                text_y,
                f"{self.scale_bar_um.value():g} µm",
                color=color,
                ha="center",
                va=vertical_alignment,
            )
        fig.savefig(filename, dpi=self.dpi.value(), bbox_inches="tight")
        plt.close(fig)
        self.status.setText(f"Exported {filename}")

    def export_spectrum(self):
        if not len(self.current_spectrum_mz):
            return
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export publication spectrum", "spectrum.svg",
            "SVG (*.svg);;PDF (*.pdf);;TIFF (*.tif *.tiff);;PNG (*.png)"
        )
        if not filename:
            return
        fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
        ax.plot(self.current_spectrum_mz, self.current_spectrum_intensity, color="black", linewidth=0.7)
        ax.set(xlabel="m/z", ylabel="Intensity", title=self.current_spectrum_title)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        fig.savefig(filename, dpi=self.dpi.value(), bbox_inches="tight")
        plt.close(fig)
        self.status.setText(f"Exported {filename}")

    def set_busy(self, busy: bool, message: str | None = None):
        if busy == self._busy:
            if message:
                self.status.setText(message)
            return
        self._busy = busy
        self.open_button.setEnabled(not busy)
        self.total_button.setEnabled(not busy and self.dataset is not None)
        self.roi_button.setEnabled(
            not busy and self.dataset is not None and self.current_display_image is not None
        )
        self.add_overlay_button.setEnabled(
            not busy and self.current_raw_image is not None
        )
        self.build_overlay_button.setEnabled(
            not busy and bool(self.overlay_layers)
        )
        self.clear_overlay_button.setEnabled(
            not busy and bool(self.overlay_layers)
        )
        if busy:
            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        else:
            # Override cursors are stacked by Qt. Clear the complete stack so an
            # earlier/nested task cannot leave Windows showing a permanent wait cursor.
            while QtWidgets.QApplication.overrideCursor() is not None:
                QtWidgets.QApplication.restoreOverrideCursor()
        if message:
            self.status.setText(message)

    def _worker_finished(self, worker: Worker):
        self._active_workers.discard(worker)
        self.set_busy(False)

    def show_error(self, message: str):
        QtWidgets.QMessageBox.critical(self, "MALDI viewer error", message)
        self.status.setText(message)


def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("MALDI imzML Viewer")
    window = MaldiViewer()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

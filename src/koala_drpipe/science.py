"""
Science data processing utilities for the KOALA data‑reduction manager
=====================================================================

This module orchestrates the application of **calibration products** (see the
companion *Calibration Set* module in this project) to science observations,
optionally aligns multiple RSS files astrometrically, estimates ADR trends, and
builds a KOALA datacube.

The central class :class:`ScienceSet` stores the file layout for one or more
science targets and exposes small, composable steps:

- :meth:`calibrate_object` – read RSS files and apply a :class:`CalibrationSet`.
- :meth:`register_object` – per‑exposure centroid registration and (optional)
  ADR estimate.
- :meth:`cube_object` – build a datacube with :class:`~pykoala.cubing.CubeInterpolator`.
- :meth:`process_object` – a convenience pipeline that chains the above.

NumPy‑style docstrings are provided for public methods.

Example
-------

>>> sciset = ScienceSet.from_config_dict({
...     "files": {"NGC1234": ["rss_1.fits", "rss_2.fits"]},
...     "is_pykoala": False,
...     "workdir": ".",
...     "cubing": {"kernel": "GaussianKernel", "kernel_kwargs": {"sigma_pix": 1.0}}
... })
>>> calset = CalibrationSet.from_config_yml("calibration.yml")
>>> out = sciset.process_object("NGC1234", calset, save_intermediate=True)
>>> cube = out["cube"]

Notes
-----
- This module intentionally does **not** handle sky/continuum subtraction of
  science frames; do that upstream or add a dedicated step if needed.
- ADR polynomials are computed for QC and optional downstream use; this module
  does not apply ADR shifts directly to the spectra.

"""
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple, Type
import os
import yaml

from astropy import units as u

from pykoala.data_container import RSS, Cube
from pykoala.instruments.koala_ifu import koala_rss
from pykoala.corrections.atmospheric_corrections import get_adr
from pykoala.corrections.astrometry import AstrometryCorrection
from pykoala import cubing

from koala_drpipe import vprint

__all__ = ["ScienceSet"]


class ScienceSet:
    """Collection of science targets and helpers to process them.

    Parameters
    ----------
    files : Mapping[str, Sequence[str]]
        Mapping from *object name* to list of RSS FITS paths.
    is_pykoala : bool or Mapping[str, bool], optional
        Whether the input files for each object are *pykoala*‑compatible RSS
        files (as opposed to raw KOALA RSS needing :func:`koala_rss`). If a
        single boolean is given, it is broadcast to all objects. Default is
        ``False``.
    workdir : str, optional
        Directory to store plots and intermediate products. Default is ``'.'``.
    kernel : type, optional
        Interpolator kernel **class** from :mod:`pykoala.cubing` (e.g.,
        ``cubing.GaussianKernel``). Default is :class:`~pykoala.cubing.GaussianKernel`.
    kernel_kwargs : mutable mapping, optional
        Keyword arguments passed to the kernel constructor inside
        :class:`~pykoala.cubing.CubeInterpolator`.

    Attributes
    ----------
    files : dict[str, list[str]]
        Mutable copy of input file mapping. May be updated when saving
        calibrated outputs.
    is_pykoala : dict[str, bool]
        Per‑object flag for how to read inputs.
    files_history : dict[str, list[list[str]]]
        History of file lists per object (updated when saving calibrated
        outputs).
    workdir : str
        Working directory for artefacts.
    kernel, kernel_kwargs
        Cubing kernel configuration.
    """

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    def __init__(
        self,
        files: Mapping[str, Sequence[str]],
        is_pykoala: bool | Mapping[str, bool] = False,
        workdir: str = ".",
        kernel: Type[Any] = cubing.GaussianKernel,
        kernel_kwargs: Optional[MutableMapping[str, Any]] = None,
    ) -> None:
        self.files: Dict[str, List[str]] = {k: list(v) for k, v in files.items()}

        if isinstance(is_pykoala, Mapping):
            self.is_pykoala: Dict[str, bool] = {k: bool(is_pykoala.get(k, False)) for k in self.files}
        else:
            self.is_pykoala = {k: bool(is_pykoala) for k in self.files}

        self.files_history: Dict[str, List[List[str]]] = {k: [] for k in self.files}
        self.workdir = workdir
        self.kernel = kernel
        self.kernel_kwargs = dict(kernel_kwargs or {})

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------
    def _read_rss_list(self, names: Sequence[str], pykoala_fmt: bool) -> List[RSS]:
        """Read a list of RSS files into memory.

        Parameters
        ----------
        names : sequence of str
            Paths to the RSS files to read.
        pykoala_fmt : bool
            Whether to read with :meth:`RSS.from_fits` (True) or
            :func:`pykoala.instruments.koala_ifu.koala_rss` (False).

        Returns
        -------
        list of RSS
            Loaded RSS objects.
        """
        reader = RSS.from_fits if pykoala_fmt else koala_rss
        return [reader(p) for p in names]

    # ------------------------------------------------------------------
    # steps
    # ------------------------------------------------------------------
    def calibrate_object(
        self,
        name: str,
        calibration_set: Any,
        calibration_corr_order: Optional[List[str]] = None,
        save: bool = False,
        overwrite: bool = False,
    ) -> Optional[List[RSS]]:
        """Apply a :class:`CalibrationSet` to all RSS exposures of *name*.

        Parameters
        ----------
        name : str
            Object name present in :attr:`files`.
        calibration_set : CalibrationSet
            Instance that bundles throughput, telluric, extinction and/or flux
            calibration corrections. Only available corrections are applied.
        calibration_corr_order : list of str, optional
            Order to apply within :meth:`CalibrationSet.apply`. If ``None``,
            the calibration set's default is used.
        save : bool, optional
            If ``True``, write calibrated RSS files to disk and update
            :attr:`files` to point to the new files. Default ``False``.
        overwrite : bool, optional
            Overwrite original files in place instead of creating ``*_cal.fits``.

        Returns
        -------
        list of RSS or None
            The calibrated RSS list if ``save=False``; otherwise ``None``.
        """
        vprint(f"Processing object: {name}")
        if name not in self.files:
            raise KeyError(f"Input object '{name}' not found")

        rss_set = self._read_rss_list(self.files[name], self.is_pykoala[name])
        rss_set = calibration_set.apply(rss_set, calibration_corr_order)

        if not save:
            return rss_set

        # Save to disk and update bookkeeping
        self.files_history[name].append(self.files[name].copy())
        self.is_pykoala[name] = True
        for i, (path, rss) in enumerate(zip(self.files[name], rss_set)):
            out_path = path if overwrite else os.path.splitext(path)[0] + "_cal.fits"
            rss.to_fits(out_path, overwrite=True)
            if not overwrite:
                self.files[name][i] = out_path
        return None

    def register_object(
        self,
        rss: Sequence[RSS],
        *,
        object_name: Optional[str] = None,
        centroider: str = "gauss",
        qc_plot: bool = True,
    ) -> Tuple[List[RSS], List[Tuple[float, float]]]:
        """Astrometrically register a list of RSS exposures.

        For each exposure, the brightest source centroid is measured and the
        exposure is shifted accordingly. An ADR polynomial trend is also
        computed per exposure for QC.

        Parameters
        ----------
        rss : sequence of RSS
            Input exposures to register.
        object_name : str, optional
            Only used for plot filenames. If ``None``, a generic name is used.
        centroider : {"gauss", "com"}
            Centroiding method passed to :class:`AstrometryCorrection`.
        qc_plot : bool, optional
            Save QC plot(s) to :attr:`workdir`.

        Returns
        -------
        registered : list of RSS
            The registered RSS exposures.
        offsets : list of tuple(float, float)
            Applied (dRA, dDec) offsets per exposure in arcsec.
        """
        name = object_name or "object"

        astrom_corr = AstrometryCorrection()
        offsets, fig = astrom_corr.register_centroids(rss, object_name=name, qc_plot=qc_plot, centroider=centroider)
        if qc_plot and fig is not None:
            fig.savefig(os.path.join(self.workdir, f"{name}_register.png"), dpi=200, bbox_inches="tight")

        out: List[RSS] = []
        for i, (ex, off) in enumerate(zip(rss, offsets)):
            astrom_corr.apply(ex, offset=off)
            # Estimate ADR for diagnostics
            _, _, adr_fig = get_adr(ex, max_adr=0.8, pol_deg=2, plot=qc_plot)
            if qc_plot and adr_fig is not None:
                adr_fig.savefig(os.path.join(self.workdir, f"{name}_adr_{i}.png"), dpi=200, bbox_inches="tight")
            out.append(ex)
        return out, list(offsets)

    def cube_object(
        self,
        rss: Sequence[RSS],
        *,
        spatial_pix_size: u.Quantity = 0.5 * u.arcsec,
        spectra_pix_size: Optional[float] = None,
        qc_plots: bool = False,
        cube_info: Optional[Dict[str, Any]] = None,
    ) -> Cube:
        """Interpolate a set of registered RSS exposures into a datacube.

        Parameters
        ----------
        rss : sequence of RSS
            Registered exposures to combine.
        spatial_pix_size : `~astropy.units.Quantity`, optional
            Cube spatial pixel size. Default ``0.5 arcsec``.
        spectra_pix_size : float, optional
            Spectral sampling in the output cube. If ``None`` (default), use the
            first RSS wavelength step.
        qc_plots : bool, optional
            Save cube QC plots.
        cube_info : dict, optional
            Extra metadata attached to the cube (e.g., ``{"name": "NGC1234"}``).

        Returns
        -------
        cube : :class:`~pykoala.data_container.Cube`
            The built datacube.
        """
        # This information should be provided by the configuration file or
        # computed if not provided.
        if spectra_pix_size is None:
            spectra_pix_size = rss[0].wavelength[1] - rss[0].wavelength[0]

        wcs = cubing.build_wcs_from_rss(
            rss, spatial_pix_size=spatial_pix_size, spectra_pix_size=spectra_pix_size
        )
        interpolator = cubing.CubeInterpolator(
            rss_set=rss,
            wcs=wcs,
            kernel=self.kernel,
            qc_plots=qc_plots,
            **self.kernel_kwargs,
        )
        cube = interpolator.build_cube(cube_info=cube_info or {})

        if qc_plots:
            if "stack_cube" in interpolator.cube_plots:
                interpolator.cube_plots["stack_cube"].savefig(
                    os.path.join(self.workdir, f"{cube.info.get('name','object')}_cube_qc.png"),
                    dpi=200,
                    bbox_inches="tight",
                )
            if "weights" in interpolator.cube_plots:
                interpolator.cube_plots["weights"].savefig(
                    os.path.join(self.workdir, f"{cube.info.get('name','object')}_cube_weights.png"),
                    dpi=200,
                    bbox_inches="tight",
                )
        return cube

    def process_object(
        self,
        name: str,
        calibration_set: Any,
        *,
        calibration_corr_order: Optional[List[str]] = None,
        save_intermediate: bool = True,
        overwrite: bool = False,
        centroider: str = "gauss",
        qc_plots: bool = True,
        spatial_pix_size: u.Quantity = 0.5 * u.arcsec,
        spectra_pix_size: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Full pipeline for a single object: calibrate → register → cube.

        Parameters
        ----------
        name : str
            Object name present in :attr:`files`.
        calibration_set : CalibrationSet
            See :meth:`calibrate_object`.
        calibration_corr_order : list of str, optional
            See :meth:`calibrate_object`.
        save_intermediate : bool, optional
            Save calibrated RSS to disk and update internal state. Default ``True``.
        overwrite : bool, optional
            Overwrite original RSS when saving intermediates.
        centroider : {"gauss", "com"}
            See :meth:`register_object`.
        qc_plots : bool, optional
            Save QC plots for registration/ADR and cube steps.
        spatial_pix_size, spectra_pix_size :
            See :meth:`cube_object`.

        Returns
        -------
        out : dict
            ``{"rss_cal": [...], "rss_reg": [...], "offsets": [...], "cube": Cube}``.
            If ``save_intermediate=True``, ``rss_cal`` will be ``None`` and
            calibrated filenames are stored in :attr:`files[name]`.
        """
        # 1) Calibrate
        rss_cal = self.calibrate_object(
            name,
            calibration_set,
            calibration_corr_order=calibration_corr_order,
            save=save_intermediate,
            overwrite=overwrite,
        )
        if save_intermediate:
            rss_cal = self._read_rss_list(self.files[name], True)

        # 2) Register
        rss_reg, offsets = self.register_object(
            rss_cal, object_name=name, centroider=centroider, qc_plot=qc_plots
        )

        # 3) Cube
        cube = self.cube_object(
            rss_reg,
            spatial_pix_size=spatial_pix_size,
            spectra_pix_size=spectra_pix_size,
            qc_plots=qc_plots,
            cube_info={"name": name},
        )

        return {
            "rss_cal": rss_cal if not save_intermediate else None,
            "rss_reg": rss_reg if not save_intermediate else None,
            "offsets": offsets,
            "cube": cube,
        }

    # ------------------------------------------------------------------
    # configuration constructors
    # ------------------------------------------------------------------
    @classmethod
    def from_config_yml(cls, yaml_file: str) -> "ScienceSet":
        """Create a :class:`ScienceSet` from a YAML configuration.

        The expected YAML section is::

            ScienceSet:
              files:
                NGC1234: ["rss_1.fits", "rss_2.fits"]
              is_pykoala: false  # or a mapping per object
              workdir: "."
              cubing:
                kernel: GaussianKernel
                kernel_kwargs: {sigma_pix: 1.0}

        Parameters
        ----------
        yaml_file : str
            Path to the configuration file.

        Returns
        -------
        ScienceSet
            A new instance.
        """
        with open(yaml_file, "r") as fh:
            config = yaml.safe_load(fh)
        if "ScienceSet" not in config:
            raise KeyError("ScienceSet key not found in YAML file")
        return cls.from_config_dict(config["ScienceSet"])

    @classmethod
    def from_config_dict(cls, config: Mapping[str, Any]) -> "ScienceSet":
        """Create a :class:`ScienceSet` from a Python mapping.

        Parameters
        ----------
        config : mapping
            A mapping with keys ``files`` (required), and optionally
            ``is_pykoala``, ``workdir``, and a ``cubing`` block containing
            ``kernel`` and ``kernel_kwargs``.

        Returns
        -------
        ScienceSet
            A new instance.
        """
        vprint("Initialising ScienceSet from config")
        files = config["files"]
        is_pykoala = config.get("is_pykoala", False)
        workdir = config.get("workdir", ".")

        # Resolve kernel class from name string (defaults to GaussianKernel)
        cubing_cfg = config.get("cubing", {})
        kernel_name = str(cubing_cfg.get("kernel", "GaussianKernel"))
        kernel_cls = getattr(cubing, kernel_name, cubing.GaussianKernel)
        kernel_kwargs = cubing_cfg.get("kernel_kwargs", {})

        return cls(
            files=files,
            is_pykoala=is_pykoala,
            workdir=workdir,
            kernel=kernel_cls,
            kernel_kwargs=kernel_kwargs,
        )

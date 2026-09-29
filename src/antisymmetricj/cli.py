"""Command-line interface for AntisymmetricJ."""

from __future__ import annotations

import argparse
from pathlib import Path

from antisymmetricj.distribution import weighted_distribution
from antisymmetricj.geometry import f_zy
from antisymmetricj.orientations import (
    OrientationSet,
    generate_zcw,
    load_orientations,
    write_orientations,
)
from antisymmetricj.plotting import (
    plot_fzy_overlay,
    plot_spectrum,
    write_fzy_distribution_table,
    write_spectrum_table,
)
from antisymmetricj.spectrum import calculate_spectrum

DEFAULT_DISTRIBUTION_ORIENTATIONS = 28_656
DEFAULT_SPECTRUM_BINS = 4096
DEFAULT_SPECTRUM_FWHM_HZ = 1.0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="antisymmetricj",
        description="Calculate antisymmetric J-coupling f_zy distributions.",
        allow_abbrev=False,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    make_zcw_set = subparsers.add_parser(
        "make_zcw_set",
        usage="antisymmetricj make_zcw_set ORIENTATIONS_OUTFILE [options]",
        help="write a legacy ZCW orientation file",
        allow_abbrev=False,
    )
    make_zcw_set.add_argument(
        "orientations_outfile",
        type=Path,
        metavar="ORIENTATIONS_OUTFILE",
    )
    make_zcw_set.add_argument("--min-orientations", type=int, default=832_039)
    make_zcw_set.set_defaults(func=_cmd_make_zcw_set)

    calc_fzy = subparsers.add_parser(
        "calc_fzy",
        usage=(
            "antisymmetricj calc_fzy [ORIENTATIONS_INFILE] "
            "DISTRIBUTION_OUTFILE [options]"
        ),
        help=(
            "calculate a weighted f_zy or f_zy^2 distribution; pass either "
            "DISTRIBUTION_OUTFILE or ORIENTATIONS_INFILE DISTRIBUTION_OUTFILE"
        ),
        allow_abbrev=False,
    )
    calc_fzy.add_argument(
        "paths",
        nargs="+",
        type=Path,
        metavar="DISTRIBUTION_OUTFILE",
        help=(
            "either DISTRIBUTION_OUTFILE for default in-memory ZCW orientations, "
            "or ORIENTATIONS_INFILE DISTRIBUTION_OUTFILE to read orientations "
            "from a file"
        ),
    )
    calc_fzy.add_argument("--plot_fzy", type=Path)
    calc_fzy.add_argument("--c-yx", type=float, default=1.0)
    calc_fzy.add_argument("--c-zx", type=float, default=1.0)
    calc_fzy.add_argument("--bins", type=int, default=200)
    calc_fzy.add_argument(
        "--min-orientations",
        type=int,
        default=DEFAULT_DISTRIBUTION_ORIENTATIONS,
        help=(
            "minimum ZCW orientations to generate when no orientation file is "
            f"supplied (default: {DEFAULT_DISTRIBUTION_ORIENTATIONS})"
        ),
    )
    calc_fzy.add_argument("--squared", action="store_true")
    calc_fzy.add_argument(
        "--format", choices=("auto", "plain", "simpson"), default="auto"
    )
    calc_fzy.add_argument(
        "--angle-unit",
        choices=("radians", "degrees"),
        default="radians",
    )
    calc_fzy.set_defaults(func=_cmd_calc_fzy)

    calc_spectrum = subparsers.add_parser(
        "calc_spectrum",
        usage=(
            "antisymmetricj calc_spectrum [ORIENTATIONS_INFILE] "
            "SPECTRUM_OUTFILE [options]"
        ),
        help=(
            "calculate a Gaussian-broadened Fig. 3-style AB NMR spectrum; pass "
            "either SPECTRUM_OUTFILE or ORIENTATIONS_INFILE SPECTRUM_OUTFILE"
        ),
        allow_abbrev=False,
    )
    calc_spectrum.add_argument(
        "paths",
        nargs="+",
        type=Path,
        metavar="SPECTRUM_OUTFILE",
        help=(
            "either SPECTRUM_OUTFILE for default in-memory ZCW orientations, "
            "or ORIENTATIONS_INFILE SPECTRUM_OUTFILE to read orientations "
            "from a file"
        ),
    )
    calc_spectrum.add_argument("--plot_spectrum", type=Path)
    calc_spectrum.add_argument("--sigma-i-hz", type=float, required=True)
    calc_spectrum.add_argument("--sigma-s-hz", type=float, required=True)
    calc_spectrum.add_argument("--j-iso-hz", type=float, required=True)
    calc_spectrum.add_argument("--j-zy-anti-hz", type=float, required=True)
    calc_spectrum.add_argument("--c-yx", type=float, default=1.0)
    calc_spectrum.add_argument("--c-zx", type=float, default=1.0)
    calc_spectrum.add_argument("--bins", type=int, default=DEFAULT_SPECTRUM_BINS)
    calc_spectrum.add_argument(
        "--fwhm-hz",
        type=float,
        default=DEFAULT_SPECTRUM_FWHM_HZ,
    )
    calc_spectrum.add_argument(
        "--min-orientations",
        type=int,
        default=DEFAULT_DISTRIBUTION_ORIENTATIONS,
        help=(
            "minimum ZCW orientations to generate when no orientation file is "
            f"supplied (default: {DEFAULT_DISTRIBUTION_ORIENTATIONS})"
        ),
    )
    calc_spectrum.add_argument(
        "--format", choices=("auto", "plain", "simpson"), default="auto"
    )
    calc_spectrum.add_argument(
        "--angle-unit",
        choices=("radians", "degrees"),
        default="radians",
    )
    calc_spectrum.set_defaults(func=_cmd_calc_spectrum)

    return parser


def _cmd_make_zcw_set(args: argparse.Namespace) -> None:
    orientations = generate_zcw(min_orientations=args.min_orientations)
    write_orientations(orientations, args.orientations_outfile)


def _cmd_calc_fzy(args: argparse.Namespace) -> None:
    orientations, distribution_outfile, orientation_source = _load_or_generate(
        args,
        command="calc_fzy",
        output_name="DISTRIBUTION_OUTFILE",
    )

    values = f_zy(
        orientations.alpha,
        orientations.beta,
        c_yx=args.c_yx,
        c_zx=args.c_zx,
    )
    distribution = weighted_distribution(
        values,
        orientations.weights,
        squared=args.squared,
        bins=args.bins,
    )
    observable = "f_zy^2" if args.squared else "f_zy"
    metadata = [
        "generated by antisymmetricj calc_fzy",
        f"observable {observable}",
        f"c_yx {args.c_yx}",
        f"c_zx {args.c_zx}",
        orientation_source,
        f"orientations {len(orientations)}",
        "normalization fraction_per_bin",
    ]
    write_fzy_distribution_table(
        distribution_outfile,
        distribution,
        metadata=metadata,
    )
    if args.plot_fzy is not None:
        plot_fzy_overlay(
            args.plot_fzy,
            [distribution],
            labels=[f"{observable}, c_yx={args.c_yx:g}, c_zx={args.c_zx:g}"],
            xlabel=f"${observable}$",
        )


def _cmd_calc_spectrum(args: argparse.Namespace) -> None:
    orientations, spectrum_outfile, orientation_source = _load_or_generate(
        args,
        command="calc_spectrum",
        output_name="SPECTRUM_OUTFILE",
    )

    values = f_zy(
        orientations.alpha,
        orientations.beta,
        c_yx=args.c_yx,
        c_zx=args.c_zx,
    )
    spectrum = calculate_spectrum(
        values,
        orientations.weights,
        sigma_i_hz=args.sigma_i_hz,
        sigma_s_hz=args.sigma_s_hz,
        j_iso_hz=args.j_iso_hz,
        j_zy_anti_hz=args.j_zy_anti_hz,
        bins=args.bins,
        fwhm_hz=args.fwhm_hz,
    )
    spectrum_limit = abs(args.j_iso_hz) + abs(args.j_zy_anti_hz)
    metadata = [
        "generated by antisymmetricj calc_spectrum",
        "observable AB spectrum",
        f"sigma_I_hz {args.sigma_i_hz}",
        f"sigma_S_hz {args.sigma_s_hz}",
        f"delta_S_hz {args.sigma_i_hz - args.sigma_s_hz}",
        f"J_iso_hz {args.j_iso_hz}",
        f"J_zy_anti_hz {args.j_zy_anti_hz}",
        f"c_yx {args.c_yx}",
        f"c_zx {args.c_zx}",
        f"fwhm_hz {args.fwhm_hz}",
        f"bins {args.bins}",
        f"spectrum_range_hz {-spectrum_limit} {spectrum_limit}",
        orientation_source,
        f"orientations {len(orientations)}",
        "normalization raw_table2_relative_intensity_times_orientation_weight",
    ]
    write_spectrum_table(
        spectrum_outfile,
        spectrum,
        metadata=metadata,
    )
    if args.plot_spectrum is not None:
        plot_spectrum(
            args.plot_spectrum,
            spectrum,
            title=(
                "AB spectrum, "
                f"J_iso={args.j_iso_hz:g} Hz, "
                f"J_zy^anti={args.j_zy_anti_hz:g} Hz"
            ),
        )


def _load_or_generate(
    args: argparse.Namespace,
    *,
    command: str,
    output_name: str,
) -> tuple[OrientationSet, Path, str]:
    if len(args.paths) == 1:
        output_file = args.paths[0]
        orientations = generate_zcw(min_orientations=args.min_orientations)
        orientation_source = f"generated_zcw_min_orientations {args.min_orientations}"
    elif len(args.paths) == 2:
        orientations_infile, output_file = args.paths
        orientations = load_orientations(
            orientations_infile,
            format=args.format,
            angle_unit=args.angle_unit,
        )
        orientation_source = f"orientations_infile {orientations_infile}"
    else:
        raise SystemExit(
            f"{command} expects either {output_name} or "
            f"ORIENTATIONS_INFILE {output_name}"
        )

    return orientations, output_file, orientation_source


if __name__ == "__main__":
    raise SystemExit(main())

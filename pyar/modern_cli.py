"""Task-oriented CLI pilot, alongside the existing PyAR entry points.

The command tree owns task discovery; each command owns its argument parser
and execution. Command modules are imported only when their task is selected.
"""

import argparse
from importlib import import_module
import sys

from pyar import __version__


# Add future tasks here without importing their scientific dependencies at startup.
_COMMANDS = {
    "energies": ("pyar.scripts.energy_table", "Print a relative-energy table from XYZ files"),
    "compare": ("pyar.scripts.compare", "Compare two structures, energies, geometry, and connectivity"),
    "clustering": ("pyar.scripts.clustering", "Cluster or filter XYZ structures"),
    "optimize": ("pyar.scripts.optimize", "Optimize each XYZ structure independently"),
    "scan-bond": ("pyar.scripts.scan_bond", "Relaxed bond scan with optional reaction-path continuation"),
}


def build_parser():
    """Build the lightweight task dispatcher."""
    parser = argparse.ArgumentParser(
        prog="pyar", description="PyAR task-oriented command interface"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")
    for name, (_, description) in _COMMANDS.items():
        # The task's shared parser handles help and every argument after its name.
        commands.add_parser(name, help=description, add_help=False)
    return parser


def main(argv=None):
    """Select a task and delegate its complete argument list unchanged."""
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv[:1])
    module_name, _ = _COMMANDS[args.command]
    command = import_module(module_name)
    entry = getattr(command, 'modern_main', command.main)
    return entry(argv[1:], prog=f"{parser.prog} {args.command}")


if __name__ == "__main__":
    main()

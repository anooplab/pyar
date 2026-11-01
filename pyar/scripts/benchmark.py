"""Developer benchmark commands, loaded only on request."""
import argparse
from importlib import import_module

TASKS = {
    'conformer': 'pyar.scripts.conformer_benchmark',
    'clustering': 'pyar.scripts.benchmark_clustering',
    'orientations': 'pyar.scripts.benchmark_orientations',
    'ts': 'pyar.scripts.ts_optimizer_benchmark',
}


def main(argv=None, *, prog=None):
    parser = argparse.ArgumentParser(add_help=False, prog=prog or 'pyar benchmark',
                                     description='Run developer benchmark campaigns.')
    parser.add_argument('task', choices=sorted(TASKS))
    if argv is None or not argv or argv[0] in {'-h', '--help'}:
        parser.print_help()
        return
    args, remaining = parser.parse_known_args(argv)
    module = import_module(TASKS[args.task])
    return module.main(remaining, prog=f'{parser.prog} {args.task}')

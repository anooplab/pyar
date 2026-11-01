"""Refresh the user CLI reference from canonical executable help."""
from pathlib import Path
import subprocess
import sys
from pyar.cli import _COMMANDS


def main():
    sections = []
    commands = [[name] for name in _COMMANDS if name not in {'gradient', 'engra'}]
    commands += [['benchmark', task] for task in ('conformer', 'clustering', 'orientations', 'ts')]
    for command in commands:
        run = subprocess.run([sys.executable, '-m', 'pyar', *command, '--help'],
                             capture_output=True, text=True, timeout=30)
        if run.returncode:
            raise RuntimeError(f"Help failed for {command}: {run.stderr}")
        sections.append('pyar ' + ' '.join(command) + '\n' + run.stdout.rstrip())
    target = Path(__file__).resolve().parents[1] / 'docs/_static/cli_reference.txt'
    target.parent.mkdir(exist_ok=True)
    target.write_text('\n\n'.join(sections) + '\n')


if __name__ == '__main__':
    main()
